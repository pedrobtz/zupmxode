# Solving (design §4.5): the deterministic core, parameters in, predictions
# out, one subject at a time: closed forms for closed-form models (D21),
# LSODA or Dormand-Prince for the rest (D5).

# The general solvers, by name (the tests loop over them).
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
  n <- length(rec$time)
  # Subjects contiguous, each keeping its rows' order; `orig` maps back to
  # the table's rows for messages.
  ids <- unique(rec$id)
  ord <- order(match(rec$id, ids), seq_len(n))
  orig <- ord
  for (k in setdiff(names(rec), "covariates")) rec[[k]] <- rec[[k]][ord]
  rec$covariates <- rec$covariates[ord, , drop = FALSE]
  par <- record_params(model, rec, params)
  sizes <- tabulate(match(rec$id, ids), length(ids))
  start <- c(0L, cumsum(sizes))
  first <- start[-length(start)] + 1L
  off <- rowsum(as.integer(rec$evid == 2L & rec$cmt < 0L), match(rec$id, ids))[, 1] > 0
  # Two choices at most, whatever the number of subjects: R allocates the
  # same in a solve of ten subjects as of ten thousand (test-alloc.R).
  plain <- subject_method(model, opts, list(evid = 0L, cmt = 1L))
  methods <- rep(plain, length(ids))
  if (any(off)) methods[off] <- subject_method(model, opts, list(evid = 2L, cmt = -1L))
  res <- .Call(zpx_solve_all, model_ptr(model), model$ir$code, model$ir$consts,
               rec$time, rec$evid, rec$amt, rec$cmt, rec$rate, rec$dur, rec$ii,
               rec$addl, rec$ss, as.integer(start), as.double(par),
               unname(method_codes[methods]),
               c(opts$rtol, opts$atol, opts$ss_tol, opts$ss_cap),
               if (model$class == "closed-form") lincmt_spec(model) else rep(-1L, 10L),
               pmx_ode_threads())
  # Errors first, in subject order, at the table's row.
  for (s in which(res$error != 0L)) {
    row <- orig[first[s] + max(res$error_row[s], 0L)]
    if (res$error[s] == 1L) abort_events(row, res$message[s])
    pmx_abort("pmx_ode_error_solver", paste0("Row ", row, ": ", res$message[s], "."), row = row)
  }
  for (s in which(res$failed_row >= 0L)) {
    pmx_warn("pmx_ode_warning_solver", paste0(
      "The ", methods[s], " solver failed for subject ", ids[s], " at time ",
      format(res$failed_at[s]), "; its outputs are NA from there."
    ), id = ids[s], time = res$failed_at[s])
  }
  for (s in which(res$ss_unconverged > 0L)) {
    row <- orig[first[s] + res$ss_row[s]]
    pmx_warn("pmx_ode_warning_steady_state", paste0(
      "A steady state for subject ", ids[s], " (row ", row, ") did not converge in ",
      opts$ss_cap, " dosing intervals."
    ), id = ids[s], row = row)
  }
  d <- as.data.frame(res$out)
  names(d) <- c("time", model$states, model$captures)
  # NaN from the C side's missing values is NA here, as R writes them.
  d[] <- lapply(d, function(x) {
    x[is.nan(x)] <- NA_real_
    x
  })
  nobs <- rowsum(as.integer(rec$evid == 0L), match(rec$id, ids))[, 1]
  out <- cbind(id = rep(ids, nobs), d)
  rownames(out) <- NULL
  out
}

# The solvers as src/core/solve.h numbers them.
method_codes <- c(closed = 0L, lsoda = 1L, dopri5 = 2L, expm = 3L)

#' The number of threads a solve uses
#'
#' Subjects are independent and are solved in parallel with OpenMP (design
#' §4.4, D8); the result is the same, bit for bit, at any thread count.
#' Without OpenMP (as on some macOS builds of R) everything runs on one
#' thread whatever is set here.
#'
#' @param n The number of threads, a positive whole number; leave it out to
#'   read the current setting.
#' @return The current setting, or, given `n`, the one it replaces,
#'   invisibly.
#' @examples
#' old <- pmx_ode_threads(2)
#' pmx_ode_threads()
#' pmx_ode_threads(old)
#' @export
pmx_ode_threads <- function(n) {
  if (missing(n)) return(.zpx$threads)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n != round(n) ||
        n > .Machine$integer.max) {
    pmx_abort("pmx_ode_error_argument", "`n` must be a whole number of threads, at least 1.",
              call = sys.call())
  }
  old <- .zpx$threads
  .zpx$threads <- as.integer(n)
  invisible(old)
}

.zpx <- new.env(parent = emptyenv())
.zpx$threads <- 1L

# The capture indices (0-based) that carry linCmt()'s rate constants, as
# src/core/solve.c's closed form reads them.
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
