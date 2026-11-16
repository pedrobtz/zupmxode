#' The solver of a model
#'
#' Every model has a default solver, the most specific its class allows
#' (design §4.3, D6): the closed form for a `linCmt()` model, the matrix
#' exponential for a linear one of up to 20 compartments, LSODA otherwise.
#' `pmx_ode_method(model)` says which a model uses; `pmx_ode_method(model,
#' method)` returns a copy that uses `method`. [pmx_ode_info()] says why.
#'
#' @param model A model.
#' @param method `"auto"` (the default choice), `"closed"` (closed-form
#'   models only), `"expm"` (closed-form or linear models of up to 20
#'   compartments), `"lsoda"` or `"dopri5"` (Dormand-Prince 5(4)).
#' @return The method in use, or, given `method`, the model with it.
#' @section Errors:
#' `pmx_ode_error_argument` for a method the model's class cannot take.
#' @examples
#' m <- pmx_ode_model({
#'   d/dt(central) <- -k * central
#' })
#' pmx_ode_method(m)
#' pmx_ode_method(pmx_ode_method(m, "lsoda"))
#' @export
pmx_ode_method <- function(model, method) {
  if (!inherits(model, "pmx_ode_model")) {
    pmx_abort("pmx_ode_error_argument", "`model` must be a model.", call = sys.call())
  }
  if (missing(method)) {
    m <- solve_options(model)$method
    return(if (identical(m, "auto")) auto_method(model)$method else m)
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
        !method %in% c("auto", solver_methods)) {
    pmx_abort("pmx_ode_error_argument", paste0(
      "`method` must be one of ", paste0("\"", c("auto", solver_methods), "\"", collapse = ", "),
      "."), call = sys.call())
  }
  ok <- switch(method,
    closed = model$class == "closed-form",
    expm = model$class %in% c("closed-form", "linear") &&
      length(model$states) <= linear_max_states,
    TRUE
  )
  if (!ok) {
    pmx_abort("pmx_ode_error_argument", paste0(
      "A ", model$class, " model with ", length(model$states),
      " compartments cannot be solved by \"", method, "\"."), call = sys.call())
  }
  model$solve_opts$method <- method
  model
}

#' What a model is and how it is solved
#'
#' @param model A model.
#' @return A list of class `pmx_ode_info`: the compartments, parameters and
#'   outputs; the class the compiler found; the solver and why it was
#'   chosen; the execution tier; the seconds lowering took; and the IR's
#'   instructions per block.
#' @examples
#' m <- pmx_ode_model({
#'   cl <- exp(tcl); v <- exp(tv)
#'   cp <- linCmt()
#' })
#' pmx_ode_info(m)
#' @export
pmx_ode_info <- function(model) {
  if (!inherits(model, "pmx_ode_model")) {
    pmx_abort("pmx_ode_error_argument", "`model` must be a model.", call = sys.call())
  }
  chosen <- solve_options(model)$method
  why <- if (identical(chosen, "auto")) auto_method(model) else
    list(method = chosen, reason = "set with pmx_ode_method()")
  model_ptr(model)
  code <- model$ir$code
  structure(list(
    compartments = model$states, parameters = model$params, outputs = model$captures,
    class = model$class, solver = why$method, reason = why$reason,
    tier = pmx_ode_tier(model), lowering_seconds = model$state$lowering_seconds,
    instructions = stats::setNames(code[7 + 3L * (seq_along(ir_blocks) - 1L) + 2L], ir_blocks)
  ), class = "pmx_ode_info")
}

#' @export
print.pmx_ode_info <- function(x, ...) {
  cat("<pmx_ode_info>\n")
  cat("  class:        ", x$class, "\n", sep = "")
  cat("  solver:       ", x$solver, " (", x$reason, ")\n", sep = "")
  cat("  tier:         ", x$tier, sprintf(", lowered in %.1f ms", 1000 * x$lowering_seconds), "\n", sep = "")
  cat("  compartments: ", paste(x$compartments, collapse = ", "), "\n", sep = "")
  if (length(x$parameters)) cat("  parameters:   ", paste(x$parameters, collapse = ", "), "\n", sep = "")
  if (length(x$outputs)) cat("  outputs:      ", paste(x$outputs, collapse = ", "), "\n", sep = "")
  cat("  instructions: ", paste(names(x$instructions), x$instructions, sep = " ", collapse = ", "), "\n", sep = "")
  invisible(x)
}
