# Solving (design §4.5): the deterministic core, parameters in, predictions
# out, one subject at a time. Closed-form models in this version of the
# package (roadmap Stage 5); the general and linear solvers follow.

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
#' value; `pmx_ode_error_solver` when a model's class has no solver yet.
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
  if (model$class != "closed-form") {
    pmx_abort("pmx_ode_error_solver", paste0(
      "This version of zupmxode solves closed-form models only; this model is ",
      model$class, "."
    ), call = sys.call())
  }
  rec <- as_records(events, model)
  ids <- unique(rec$id)
  parts <- lapply(ids, function(id) {
    keep <- rec$id == id
    one <- lapply(rec[setdiff(names(rec), "covariates")], function(x) x[keep])
    one$covariates <- rec$covariates[keep, , drop = FALSE]
    rows <- which(keep)
    res <- solve_closed_one(model, one, params)
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

solve_closed_one <- function(model, rec, params) {
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
  .Call(zpx_solve_closed, model_ptr(model), rec$time, rec$evid, rec$amt, rec$cmt,
        rec$rate, rec$dur, rec$ii, rec$addl, rec$ss, pre[, 1], pre[, 2], pre[, 3],
        pre[, 4], as.double(par), lincmt_spec(model), as.double(init))
}
