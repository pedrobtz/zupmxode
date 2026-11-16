# Solving (design §4.5): the deterministic core, parameters in, predictions
# out, one subject at a time: closed forms for closed-form models (D21),
# LSODA or Dormand-Prince for the rest (D5).

# The solvers, by name, as zpx_solve_general() numbers them.
solver_codes <- c(lsoda = 1L, dopri5 = 2L)

# The solvers a model can take, most specific first (design §4.3, D6).
solver_methods <- c("closed", "expm", "lsoda", "dopri5")

# The linear solver's state limit: above it, a linear model goes to LSODA
# (design §4.3, "about 20"; fixed at 20 in Stage 7).
linear_max_states <- 20L

# Options the solve takes from the model: the method (pmx_ode_method()),
# the general solvers' tolerances (design §6.3), and steady state's
# convergence tolerance and cap (iterations).
solve_options <- function(model) {
  o <- list(method = "auto", rtol = 1e-6, atol = 1e-8, ss_tol = 1e-8, ss_cap = 1000)
  o[names(model$solve_opts)] <- model$solve_opts
  # The tests' old switch for the closed form's ODE twin.
  if (isFALSE(o$closed) && identical(o$method, "auto")) o$method <- "lsoda"
  o
}

# The solver `auto` chooses for a model, and why.
auto_method <- function(model) {
  ny <- length(model$states)
  if (model$class == "closed-form") {
    return(list(method = "closed", reason = "a linCmt() model with no other derivative"))
  }
  if (model$class == "linear" && ny <= linear_max_states) {
    return(list(method = "expm", reason = paste0(
      "linear in its compartments, with coefficients constant between events (",
      ny, " of at most ", linear_max_states, " compartments)")))
  }
  if (model$class == "linear") {
    return(list(method = "lsoda", reason = paste0(
      "linear, but ", ny, " compartments is above the matrix exponential's ",
      linear_max_states)))
  }
  list(method = "lsoda", reason = "no closed form and not linear in its compartments")
}

# The solver for one subject: the model's method, or auto's choice; a
# compartment switched off takes auto from the closed form or the matrix
# exponential to LSODA, which supports it.
subject_method <- function(model, opts, rec) {
  off <- any(rec$evid == 2L & rec$cmt < 0L)
  if (!identical(opts$method, "auto")) {
    if (off && opts$method %in% c("closed", "expm")) {
      pmx_abort("pmx_ode_error_solver", paste0(
        "The ", opts$method, " solver does not support a compartment switched off ",
        "(EVID 2 with a negative CMT); use \"lsoda\" or \"dopri5\"."))
    }
    return(opts$method)
  }
  m <- auto_method(model)$method
  if (off && m %in% c("closed", "expm")) "lsoda" else m
}

#' Solve a model for one or more subjects
#'
#' Runs `model` over the event table `events` with the parameters
#' `params`, and returns the compartments and the model's outputs at every
#' observation record. Nothing is drawn at random: this is the
#' deterministic core that the population simulation builds on.
#'
#' @param model A model from [pmx_ode_model()].
#' @param events An event table: from [pmx_ode_events()], or a NONMEM-style
#'   data frame with columns `ID`, `TIME`, `EVID`, `AMT`, `CMT`, `RATE`,
#'   `DUR`, `II`, `ADDL` and `SS` (any case; all but `TIME` optional).
#'   Columns named like a model parameter are covariates, each record
#'   carrying its value.
#' @param params The model's parameters, a named numeric vector. A
#'   covariate column of `events` takes the place of its parameter.
#' @return A data frame with `id`, `time`, one column per compartment and
#'   one per model output, one row per observation record.
#' @section Errors:
#' `pmx_ode_error_events` for a record that breaks a rule of the event
#' table, naming the row; `pmx_ode_error_argument` for a parameter with no
#' value; `pmx_ode_error_solver` for a closed-form model whose rate
#' constants are not positive.
#' @section Warnings:
#' `pmx_ode_warning_solver` when the solver fails for a subject: that
#' subject's outputs are `NA` from where it failed, and the others are
#' solved (design §8). `pmx_ode_warning_steady_state` when a steady state
#' of a model with no closed form does not converge.
#' @examples
#' m <- pmx_ode_model({
#'   cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
#'   cp <- linCmt()
#' })
#' ev <- pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
#'                      time = seq(0, 168, by = 12))
#' head(pmx_ode_solve(m, ev, params = c(tcl = 0, tv = 3, tka = -0.7)))
#' @export
pmx_ode_solve <- function(model, events, params = double()) {
  if (!inherits(model, "pmx_ode_model")) {
    pmx_abort("pmx_ode_error_argument", "`model` must be a model.", call = sys.call())
  }
  if (!is.numeric(params) || (length(params) && is.null(names(params)))) {
    pmx_abort("pmx_ode_error_argument", "`params` must be a named numeric vector.",
              call = sys.call())
  }
  rec <- as_records(events, model)
  opts <- solve_options(model)
  ids <- unique(rec$id)
  parts <- lapply(ids, function(id) {
    keep <- rec$id == id
    one <- lapply(rec[setdiff(names(rec), "covariates")], function(x) x[keep])
    one$covariates <- rec$covariates[keep, , drop = FALSE]
    rows <- which(keep)
    method <- subject_method(model, opts, one)
    res <- solve_one(model, one, params, method, opts)
    if (!is.null(res$failed_row) && res$failed_row >= 0) {
      pmx_warn("pmx_ode_warning_solver", paste0(
        "The ", method, " solver failed for subject ", id, " at time ",
        format(res$failed_at), "; its outputs are NA from there."
      ), id = id, time = res$failed_at)
    }
    if (!is.null(res$ss_unconverged) && res$ss_unconverged > 0) {
      pmx_warn("pmx_ode_warning_steady_state", paste0(
        "A steady state for subject ", id, " (row ", rows[res$ss_row + 1L],
        ") did not converge in ", opts$ss_cap, " dosing intervals."
      ), id = id, row = rows[res$ss_row + 1L])
    }
    if (!is.null(res$error)) {
      if (res$error == "events") abort_events(rows[res$row + 1L], res$message)
      pmx_abort("pmx_ode_error_solver", paste0("Row ", rows[res$row + 1L], ": ",
                                               res$message, "."),
                row = rows[res$row + 1L])
    }
    d <- as.data.frame(res$out)
    names(d) <- c("time", model$states, model$captures)
    cbind(id = rep(id, nrow(d)), d)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# The capture indices (0-based) that carry linCmt()'s rate constants, as
# zpx_solve_closed() reads them.
lincmt_spec <- function(model) {
  lc <- model$lincmt
  p <- lc$params
  at <- function(name) {
    if (is.na(name)) return(-1L)
    match(name, model$captures) - 1L
  }
  keys <- if (lc$micro) c("ka", "k", "v", "k12", "k21", "k13", "k31") else
    c("ka", "cl", "v", "q", "vp", "q2", "vp2")
  as.integer(c(lc$ncmt, lc$depot, lc$micro, vapply(keys, function(k) at(p[[k]]), 0L)))
}

solve_one <- function(model, rec, params, method, opts) {
  par <- record_params(model, rec, params)
  ny <- length(model$states)
  n <- length(rec$time)
  pre <- matrix(c(0, 1, 0, 0), nrow = n, ncol = 4L, byrow = TRUE)
  init <- double(ny)
  for (i in seq_len(n)) {
    if (i == 1L || rec$evid[i] %in% c(1L, 4L)) {
      out <- matrix(model_eval(model, "pre", par = par[i, ], y = double(ny))$out, nrow = ny)
      if (i == 1L) init <- out[, 5]
      c <- rec$cmt[i]
      if (rec$evid[i] %in% c(1L, 4L) && c >= 1L && c <= ny) pre[i, ] <- out[c, 1:4]
    }
  }
  if (method == "expm") {
    return(.Call(zpx_solve_linear, model_ptr(model), rec$time, rec$evid, rec$amt, rec$cmt,
                 rec$rate, rec$dur, rec$ii, rec$addl, rec$ss, pre[, 1], pre[, 2], pre[, 3],
                 pre[, 4], as.double(par), as.double(init)))
  }
  if (method == "closed") {
    return(.Call(zpx_solve_closed, model_ptr(model), rec$time, rec$evid, rec$amt, rec$cmt,
                 rec$rate, rec$dur, rec$ii, rec$addl, rec$ss, pre[, 1], pre[, 2], pre[, 3],
                 pre[, 4], as.double(par), lincmt_spec(model), as.double(init)))
  }
  .Call(zpx_solve_general, model_ptr(model), rec$time, rec$evid, rec$amt, rec$cmt,
        rec$rate, rec$dur, rec$ii, rec$addl, rec$ss, pre[, 1], pre[, 2], pre[, 3],
        pre[, 4], as.double(par), as.double(init),
        c(solver_codes[[method]], opts$rtol, opts$atol, opts$ss_tol, opts$ss_cap))
}
