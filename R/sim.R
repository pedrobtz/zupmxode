# Population simulation (design §4.5, D7, D19): etas from omega and
# residual error drawn in R with R's RNG, every subject solved by one call
# to the deterministic core, so a simulation is reproducible from
# set.seed() at any thread count.

#' Simulate a population
#'
#' Draws each subject's etas from a multivariate normal with covariance
#' `omega`, solves every subject at `theta` and its etas with
#' [pmx_ode_solve()], and adds residual error to the outputs `sigma` names.
#' Every draw is made in R with R's random number generator, so
#' [set.seed()] makes a simulation reproducible, at any thread count.
#'
#' The draws are made in a fixed order: first the etas, as
#' `MASS::mvrnorm(n, 0, omega)` would draw them, one row per subject; then,
#' for each entry of `sigma` in turn, one standard normal per observation
#' record. A subject the solver fails for still takes its draws, so the
#' others' are unchanged.
#'
#' Four residual error models are built in: [add()], [prop()], their sum
#' `add(a) + prop(b)` (the combined model) and [lnorm()]. Any other is a
#' function applied to the predictions this returns.
#'
#' @inheritParams pmx_ode_solve
#' @param theta The model's fixed parameters, a named numeric vector.
#' @param omega The covariance of the etas: a symmetric positive
#'   semi-definite matrix, or a named vector of variances for a diagonal
#'   one. Its names, or with none the model's parameters whose names start
#'   with `eta`, in the order the model uses them, say which parameters are
#'   etas. `NULL` draws none.
#' @param sigma Residual error: a named list of error models, from [add()],
#'   [prop()] or [lnorm()], each name a compartment or an output of the
#'   model. `NULL` adds none.
#' @param n The number of subjects. A table of one subject is copied `n`
#'   times, the copies numbered 1 to `n`; a table of several subjects is
#'   simulated once per subject, and `n`, if given, must be their number.
#' @return A data frame with [pmx_ode_solve()]'s columns, then one column per
#'   eta holding the subject's draw, then, for each entry of `sigma`, the
#'   output with residual error, named after it with `_obs` appended.
#' @section Errors:
#' `pmx_ode_error_argument` for an argument that is not what is described
#' here; the errors of [pmx_ode_solve()].
#' @seealso [add()] for the error models.
#' @examples
#' m <- pmx_ode_model({
#'   ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
#'   cp <- linCmt()
#' })
#' ev <- pmx_ode_events(amt = 100, cmt = "depot", time = c(1, 2, 4, 8, 12, 24))
#' set.seed(1)
#' sim <- pmx_ode_sim(m, ev, theta = c(tka = 0.5, tcl = 1, tv = 3),
#'                    omega = diag(c(0.1, 0.1)), sigma = list(cp = prop(0.1)),
#'                    n = 10)
#' head(sim)
#' @export
pmx_ode_sim <- function(
  model,
  events,
  theta = double(),
  omega = NULL,
  sigma = NULL,
  n = NULL
) {
  call <- sys.call()
  if (!inherits(model, "pmx_ode_model")) {
    pmx_abort("pmx_ode_error_argument", "`model` must be a model.", call = call)
  }
  if (!is.numeric(theta) || (length(theta) && is.null(names(theta)))) {
    pmx_abort(
      "pmx_ode_error_argument",
      "`theta` must be a named numeric vector.",
      call = call
    )
  }
  omega <- check_omega(omega, model, call)
  etas <- colnames(omega)
  sigma <- check_sigma(sigma, model, call)
  if (any(names(theta) %in% etas)) {
    pmx_abort(
      "pmx_ode_error_argument",
      paste0(
        "`theta` gives ",
        paste0("`", intersect(names(theta), etas), "`", collapse = ", "),
        ", which `omega` draws."
      ),
      call = call
    )
  }
  if (!is.data.frame(events)) {
    pmx_abort(
      "pmx_ode_error_argument",
      "`events` must be a data frame.",
      call = call
    )
  }
  d <- as.data.frame(events, stringsAsFactors = FALSE)
  clash <- intersect(names(d), etas)
  if (length(clash)) {
    pmx_abort(
      "pmx_ode_error_argument",
      paste0(
        "`events` has a column ",
        paste0("`", clash, "`", collapse = ", "),
        ", which `omega` draws."
      ),
      call = call
    )
  }
  idcol <- which(toupper(names(d)) == "ID")[1]
  id <- if (is.na(idcol)) rep(1L, nrow(d)) else d[[idcol]]
  ids <- unique(id)
  if (!is.null(n)) {
    if (
      !is.numeric(n) ||
        length(n) != 1L ||
        is.na(n) ||
        n < 1 ||
        n != round(n) ||
        n > .Machine$integer.max
    ) {
      pmx_abort(
        "pmx_ode_error_argument",
        "`n` must be a whole number of subjects, at least 1.",
        call = call
      )
    }
    n <- as.integer(n)
    if (length(ids) > 1L && n != length(ids)) {
      pmx_abort(
        "pmx_ode_error_argument",
        paste0(
          "`events` holds ",
          length(ids),
          " subjects, so `n` must be ",
          length(ids),
          " or left out."
        ),
        call = call
      )
    }
  }
  if (length(ids) == 1L && !is.null(n) && n > 1L) {
    rows <- nrow(d)
    d <- d[rep(seq_len(rows), n), , drop = FALSE]
    id <- rep(seq_len(n), each = rows)
    if (is.na(idcol)) {
      d$ID <- id
    } else {
      d[[idcol]] <- id
    }
    rownames(d) <- NULL
    ids <- seq_len(n)
  }
  ns <- length(ids)
  # The draws, in the documented order: etas, then each sigma entry.
  eta <- draw_etas(ns, omega)
  subject <- match(id, ids)
  for (j in seq_along(etas)) {
    d[[etas[j]]] <- eta[subject, j]
  }
  out <- pmx_ode_solve(model, d, theta)
  row_subject <- match(out$id, ids)
  for (j in seq_along(etas)) {
    out[[etas[j]]] <- eta[row_subject, j]
  }
  for (name in names(sigma)) {
    e <- stats::rnorm(nrow(out))
    out[[paste0(name, "_obs")]] <- apply_error(sigma[[name]], out[[name]], e)
  }
  out
}

# omega as a matrix with the etas' names on both sides (zero columns for
# NULL).
check_omega <- function(omega, model, call) {
  bad <- function(why) {
    pmx_abort(
      "pmx_ode_error_argument",
      paste0("`omega` ", why, "."),
      call = call
    )
  }
  if (is.null(omega)) {
    return(matrix(0, 0, 0, dimnames = list(character(), character())))
  }
  if (!is.numeric(omega)) {
    bad("must be a numeric matrix or a named vector of variances")
  }
  if (!is.matrix(omega)) {
    if (is.null(names(omega))) {
      bad("given as a vector must be named")
    }
    v <- omega
    omega <- diag(as.double(v), length(v), length(v), names = FALSE)
    dimnames(omega) <- list(names(v), names(v))
  }
  k <- nrow(omega)
  if (ncol(omega) != k) {
    bad("must be square")
  }
  if (any(!is.finite(omega))) {
    bad("must hold finite numbers")
  }
  if (!isTRUE(all.equal(omega, t(omega), check.attributes = FALSE))) {
    bad("must be symmetric")
  }
  if (
    !is.null(rownames(omega)) &&
      !is.null(colnames(omega)) &&
      !identical(rownames(omega), colnames(omega))
  ) {
    bad("has row and column names that differ")
  }
  names <- rownames(omega)
  if (is.null(names)) {
    names <- colnames(omega)
  }
  if (is.null(names)) {
    cand <- grep("^eta", model$params, value = TRUE)
    if (length(cand) != k) {
      bad(paste0(
        "has ",
        k,
        " rows and no names, and the model has ",
        length(cand),
        " parameters named eta...; name its rows"
      ))
    }
    names <- cand
  }
  if (anyDuplicated(names)) {
    bad("names an eta twice")
  }
  unknown <- setdiff(names, model$params)
  if (length(unknown)) {
    bad(paste0(
      "names ",
      paste0("`", unknown, "`", collapse = ", "),
      ", not a parameter of the model"
    ))
  }
  # An event column's name is read as that column, never as a parameter.
  taken <- names[toupper(names) %in% c(event_columns, "ID")]
  if (length(taken)) {
    bad(paste0(
      "names ",
      paste0("`", taken, "`", collapse = ", "),
      ", which an event table reads as its own column"
    ))
  }
  if (k) {
    ev <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-6 * max(abs(ev), 1)) {
      bad("must be positive semi-definite")
    }
  }
  dimnames(omega) <- list(names, names)
  omega
}

# n draws from N(0, omega), one row per subject, as MASS::mvrnorm() makes
# them from the same seed (and nothing drawn when omega is empty).
draw_etas <- function(n, omega) {
  p <- ncol(omega)
  if (p == 0L) {
    return(matrix(0, n, 0))
  }
  es <- eigen(omega, symmetric = TRUE)
  x <- matrix(stats::rnorm(p * n), n)
  x <- 0 + es$vectors %*% diag(sqrt(pmax(es$values, 0)), p) %*% t(x)
  t(x)
}

check_sigma <- function(sigma, model, call) {
  if (is.null(sigma)) {
    return(list())
  }
  if (
    !is.list(sigma) ||
      inherits(sigma, "pmx_ode_resid") ||
      is.null(names(sigma)) ||
      any(names(sigma) == "")
  ) {
    pmx_abort(
      "pmx_ode_error_argument",
      "`sigma` must be a named list of error models, such as list(cp = prop(0.1)).",
      call = call
    )
  }
  if (anyDuplicated(names(sigma))) {
    pmx_abort(
      "pmx_ode_error_argument",
      "`sigma` names an output twice.",
      call = call
    )
  }
  outputs <- c(model$states, model$captures)
  unknown <- setdiff(names(sigma), outputs)
  if (length(unknown)) {
    pmx_abort(
      "pmx_ode_error_argument",
      paste0(
        "`sigma` names ",
        paste0("`", unknown, "`", collapse = ", "),
        ", not a compartment or output of the model (",
        paste(outputs, collapse = ", "),
        ")."
      ),
      call = call
    )
  }
  for (name in names(sigma)) {
    if (!inherits(sigma[[name]], "pmx_ode_resid")) {
      pmx_abort(
        "pmx_ode_error_argument",
        paste0(
          "`sigma$",
          name,
          "` must be an error model: add(), prop(), ",
          "add() + prop() or lnorm()."
        ),
        call = call
      )
    }
  }
  sigma
}

# The prediction f with the residual draw e under the error model r.
apply_error <- function(r, f, e) {
  switch(
    r$type,
    add = f + r$add * e,
    prop = f + f * r$prop * e,
    combined = f + sqrt(r$add^2 + r$prop^2 * f^2) * e,
    lnorm = f * exp(r$lnorm * e)
  )
}

#' Residual error models
#'
#' The four error models [pmx_ode_sim()] applies to a prediction `f`, each
#' with `e` a standard normal draw:
#'
#' * `add(a)`: additive, `f + a * e`;
#' * `prop(b)`: proportional, `f + f * b * e`;
#' * `add(a) + prop(b)`: combined, `f + sqrt(a^2 + b^2 * f^2) * e`, the
#'   variances added;
#' * `lnorm(s)`: exponential, `f * exp(s * e)`.
#'
#' The spellings are rxode2's.
#'
#' @param sd The standard deviation: a single non-negative finite number.
#' @return An error model, for the `sigma` argument of [pmx_ode_sim()].
#' @examples
#' prop(0.1)
#' add(0.5) + prop(0.1)
#' @export
add <- function(sd) resid_model("add", sd, sys.call())

#' @rdname add
#' @export
prop <- function(sd) resid_model("prop", sd, sys.call())

#' @rdname add
#' @export
lnorm <- function(sd) resid_model("lnorm", sd, sys.call())

resid_model <- function(type, sd, call) {
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd < 0) {
    pmx_abort(
      "pmx_ode_error_argument",
      "The standard deviation must be one non-negative finite number.",
      call = call
    )
  }
  r <- list(type = type)
  r[[type]] <- as.double(sd)
  structure(r, class = "pmx_ode_resid")
}

#' @export
Ops.pmx_ode_resid <- function(e1, e2) {
  if (missing(e2)) {
    pmx_abort(
      "pmx_ode_error_argument",
      "Error models combine only as add(a) + prop(b).",
      call = sys.call()
    )
  }
  types <- c(
    if (inherits(e1, "pmx_ode_resid")) e1$type,
    if (inherits(e2, "pmx_ode_resid")) e2$type
  )
  if (
    .Generic != "+" || !setequal(types, c("add", "prop")) || length(types) != 2L
  ) {
    pmx_abort(
      "pmx_ode_error_argument",
      "Error models combine only as add(a) + prop(b).",
      call = sys.call()
    )
  }
  a <- if (e1$type == "add") e1 else e2
  b <- if (e1$type == "prop") e1 else e2
  structure(
    list(type = "combined", add = a$add, prop = b$prop),
    class = "pmx_ode_resid"
  )
}

#' @export
print.pmx_ode_resid <- function(x, ...) {
  cat(
    switch(
      x$type,
      add = paste0("additive error, sd ", format(x$add)),
      prop = paste0("proportional error, sd ", format(x$prop)),
      combined = paste0(
        "combined error, additive sd ",
        format(x$add),
        ", proportional sd ",
        format(x$prop)
      ),
      lnorm = paste0("exponential error, sd ", format(x$lnorm))
    ),
    "\n",
    sep = ""
  )
  invisible(x)
}
