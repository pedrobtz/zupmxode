# Per-subject evaluation (design §4.6, D9, D12): subjects prepared once,
# then evaluated at any theta and eta, one call each, as an estimator does.
# The same call is zupmxode_subject_eval() in C (inst/include/zupmxode.h).

#' Prepare subjects for repeated evaluation
#'
#' Reads an event table once and keeps each subject's records ready for
#' [pmx_ode_subject_eval()], which an estimator calls at every parameter
#' vector it tries. The same subjects serve the C API in
#' `inst/include/zupmxode.h`.
#'
#' The model's parameters split three ways. Those named in `eta` are the
#' subject's random effects; columns of `data` named like a parameter are
#' covariates; the rest are theta, in the order the model uses them. The
#' split is stored on the result as the attributes `theta` and `eta`.
#'
#' @inheritParams pmx_ode_solve
#' @param data An event table, as for [pmx_ode_solve()].
#' @param eta The names of the parameters that are etas. By default the
#'   model's parameters whose names start with `eta` that are not covariates.
#' @return A list of class `pmx_ode_subjects`, one prepared subject per ID
#'   in the order they first appear, each a list with `id`, `time` (its
#'   observation times) and the subject itself. The attributes `theta` and
#'   `eta` name the parameters [pmx_ode_subject_eval()] expects. A prepared
#'   subject lives for the session: one saved and loaded again must be
#'   prepared again.
#' @section Errors:
#' `pmx_ode_error_argument` for an `eta` that is not a parameter of the
#' model or is a covariate; `pmx_ode_error_events` as for [pmx_ode_solve()].
#' @seealso [pmx_ode_subject_eval()].
#' @examples
#' m <- pmx_ode_model({
#'   cl <- exp(tcl + eta.cl); v <- exp(tv); ka <- exp(tka)
#'   cp <- linCmt()
#' })
#' ev <- pmx_ode_events(amt = 100, cmt = "depot", time = c(1, 4, 12, 24))
#' subs <- pmx_ode_subjects(m, ev)
#' attr(subs, "theta")
#' pmx_ode_subject_eval(m, subs[[1]], theta = c(tcl = 0, tv = 3, tka = -0.5),
#'                      eta = c(eta.cl = 0.1), sigma = list(cp = prop(0.1)))
#' @export
pmx_ode_subjects <- function(model, data, eta = NULL) {
  call <- sys.call()
  if (!inherits(model, "pmx_ode_model")) {
    pmx_abort("pmx_ode_error_argument", "`model` must be a model.", call = call)
  }
  rec <- as_records(data, model)
  covs <- names(rec$covariates)
  if (is.null(eta)) {
    eta <- setdiff(grep("^eta", model$params, value = TRUE), covs)
  }
  if (!is.character(eta) || anyNA(eta) || anyDuplicated(eta)) {
    pmx_abort("pmx_ode_error_argument", "`eta` must be distinct parameter names.", call = call)
  }
  bad <- c(setdiff(eta, model$params), intersect(eta, covs))
  if (length(bad)) {
    pmx_abort("pmx_ode_error_argument", paste0(
      "`eta` names ", paste0("`", bad, "`", collapse = ", "),
      ", not a parameter of the model or a covariate of `data`."), call = call)
  }
  eta <- model$params[model$params %in% eta]
  theta <- setdiff(model$params, c(eta, covs))
  kind <- ifelse(model$params %in% theta, 0L, ifelse(model$params %in% eta, 1L, 2L))
  index <- vapply(seq_along(model$params), function(p) {
    name <- model$params[p]
    as.integer(switch(kind[p] + 1L, match(name, theta), match(name, eta), match(name, covs)) - 1L)
  }, 0L)
  opts <- solve_options(model)
  lc <- if (model$class == "closed-form") lincmt_spec(model) else rep(-1L, 10L)
  ptr <- model_ptr(model)
  ids <- unique(rec$id)
  subs <- lapply(ids, function(id) {
    i <- which(rec$id == id)
    off <- any(rec$evid[i] == 2L & rec$cmt[i] < 0L)
    method <- subject_method(model, opts, list(evid = rec$evid[i], cmt = rec$cmt[i]))
    allowed <- c(closed = model$class == "closed-form" && !off,
                 lsoda = TRUE, dopri5 = TRUE,
                 expm = model$class %in% c("closed-form", "linear") &&
                   length(model$states) <= linear_max_states && !off)
    cov <- matrix(as.double(unlist(rec$covariates[i, , drop = FALSE], use.names = FALSE)),
                  nrow = length(i), ncol = length(covs))
    p <- .Call(zpx_subject_new, ptr, rec$time[i], rec$evid[i], rec$amt[i], rec$cmt[i],
               rec$rate[i], rec$dur[i], rec$ii[i], rec$addl[i], rec$ss[i], cov, kind, index,
               c(length(theta), length(eta)), as.integer(c(method_codes[[method]], sum(2^(0:3)[allowed]))),
               c(opts$rtol, opts$atol, opts$ss_tol, opts$ss_cap), lc)
    structure(list(id = id, time = rec$time[i][rec$evid[i] == 0L], rows = i, ptr = p,
                   theta = theta, eta = eta, method = method),
              class = "pmx_ode_subject")
  })
  structure(subs, class = "pmx_ode_subjects", theta = theta, eta = eta)
}

#' @export
`[.pmx_ode_subjects` <- function(x, i) {
  structure(unclass(x)[i], class = "pmx_ode_subjects", theta = attr(x, "theta"),
            eta = attr(x, "eta"))
}

#' @export
print.pmx_ode_subjects <- function(x, ...) {
  cat("<pmx_ode_subjects> ", length(x), " subject", if (length(x) != 1L) "s", "\n", sep = "")
  cat("  theta: ", paste(attr(x, "theta"), collapse = ", "), "\n", sep = "")
  if (length(attr(x, "eta"))) cat("  eta:   ", paste(attr(x, "eta"), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Evaluate one prepared subject
#'
#' One subject's predictions and residual variances at its observation
#' records, at `theta` and `eta`: what an estimator needs to form the
#' likelihood. The predictions are the doubles [pmx_ode_solve()] gives for
#' the same subject and parameters. This is `zupmxode_subject_eval()` of the
#' C API (`inst/include/zupmxode.h`), one R call per subject.
#'
#' The variance follows the error model, `f` the prediction: `add(a)`
#' gives `a^2`, `prop(b)` gives `(b * f)^2`, `add(a) + prop(b)` their sum,
#' and `lnorm(s)` gives `s^2`, the variance on the log scale.
#'
#' @inheritParams pmx_ode_solve
#' @param subject One element of [pmx_ode_subjects()]'s result.
#' @param theta,eta Named numeric vectors with the parameters the subject
#'   was prepared with (`attr(subjects, "theta")` and `attr(subjects,
#'   "eta")`), in any order.
#' @param sigma The output and its error model, a list of one element, as
#'   in [pmx_ode_sim()]: `list(cp = prop(0.1))`.
#' @param tol The general solvers' tolerances, `c(rtol, atol)`; `NULL` for
#'   the model's.
#' @return A list with `time`, `pred` and `var`, one value each per
#'   observation record.
#' @section Errors:
#' `pmx_ode_error_argument` for arguments that do not match the subject
#' and model, or a subject saved and loaded again; `pmx_ode_error_events`
#' and `pmx_ode_error_solver` as for [pmx_ode_solve()].
#' @section Warnings:
#' `pmx_ode_warning_solver` when the solver fails: `pred` and `var` are
#' `NA` from there. `pmx_ode_warning_steady_state` when a steady state does
#' not converge.
#' @inherit pmx_ode_subjects examples
#' @export
pmx_ode_subject_eval <- function(model, subject, theta = double(), eta = double(), sigma,
                                 tol = NULL) {
  call <- sys.call()
  arg <- function(msg) pmx_abort("pmx_ode_error_argument", msg, call = call)
  if (!inherits(model, "pmx_ode_model")) arg("`model` must be a model.")
  if (!inherits(subject, "pmx_ode_subject")) {
    arg("`subject` must be one element of pmx_ode_subjects()'s result.")
  }
  if (!.Call(zpx_subject_valid, subject$ptr)) {
    arg("`subject` was saved and loaded again; prepare it again with pmx_ode_subjects().")
  }
  pick <- function(x, names, what) {
    if (!is.numeric(x) || length(x) != length(names) ||
          (length(x) && (is.null(names(x)) || !setequal(names(x), names) ||
                           anyDuplicated(names(x))))) {
      arg(paste0("`", what, "` must be a numeric vector named ",
                 if (length(names)) paste0("`", names, "`", collapse = ", ") else "nothing", "."))
    }
    as.double(x[names])
  }
  th <- pick(theta, subject$theta, "theta")
  et <- pick(eta, subject$eta, "eta")
  sigma <- check_sigma(sigma, model, call)
  if (length(sigma) != 1L) arg("`sigma` must name one output.")
  r <- sigma[[1]]
  err <- c(add = 1, prop = 2, combined = 3, lnorm = 4)[[r$type]]
  sd_add <- if (r$type == "lnorm") r$lnorm else if (is.null(r$add)) 0 else r$add
  sd_prop <- if (is.null(r$prop)) 0 else r$prop
  output <- match(names(sigma), c(model$states, model$captures)) - 1L
  if (is.null(tol)) tol <- c(0, 0)
  if (!is.numeric(tol) || length(tol) != 2L || anyNA(tol) || any(tol <= 0 & tol != 0)) {
    arg("`tol` must be c(rtol, atol), both positive.")
  }
  res <- .Call(zpx_subject_eval_r, model_ptr(model), subject$ptr, th, et,
               c(output, err, sd_add, sd_prop, -1, tol))
  st <- res$status
  if (st == 17L) {
    arg("`subject` was prepared for another model.")
  }
  if (st == 20L) {
    abort_events(subject$rows[res$error_row + 1L], res$message)
  }
  if (st == 21L) {
    row <- subject$rows[max(res$error_row, 0L) + 1L]
    pmx_abort("pmx_ode_error_solver", paste0("Row ", row, ": ", res$message, "."), row = row)
  }
  if (st >= 16L) {
    pmx_abort("pmx_ode_error_internal", paste0("zupmxode_subject_eval() returned ", st, ": ",
                                               .Call(zpx_status_message_r, st), "."))
  }
  if (st == 1L) {
    pmx_warn("pmx_ode_warning_solver", paste0(
      "The ", subject$method, " solver failed for subject ", subject$id, " at time ",
      format(res$failed_at), "; its outputs are NA from there."
    ), id = subject$id, time = res$failed_at)
  }
  if (st == 2L) {
    pmx_warn("pmx_ode_warning_steady_state", paste0(
      "A steady state for subject ", subject$id, " did not converge."), id = subject$id)
  }
  na <- function(x) {
    x[is.nan(x)] <- NA_real_
    x
  }
  list(time = subject$time, pred = na(res$pred), var = na(res$var))
}

# The tests' C callers (src/zpx_api.c): every subject through the
# registered table from `threads` threads; and one failing call by name.
test_eval_threads <- function(model, subjects, theta, eta, sigma, threads) {
  s <- subjects[[1]]
  r <- sigma[[1]]
  err <- c(add = 1, prop = 2, combined = 3, lnorm = 4)[[r$type]]
  output <- match(names(sigma), c(model$states, model$captures)) - 1L
  opts <- c(output, err, if (is.null(r$add)) 0 else r$add, if (is.null(r$prop)) 0 else r$prop,
            -1, 0, 0)
  .Call(zpx_test_eval_threads, model_ptr(model), lapply(subjects, `[[`, "ptr"),
        as.double(theta[s$theta]), as.double(eta[s$eta]), opts, as.integer(threads))
}

test_c_api <- function(model, subject, small, what) {
  .Call(zpx_test_c_api, model_ptr(model), subject$ptr, small$ptr, what)
}
