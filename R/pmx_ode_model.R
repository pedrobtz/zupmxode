#' Compile a model
#'
#' Reads a model written in the subset of rxode2's model language that
#' zupmxode runs (design §6.1), and compiles it to native code through
#' MIR, with no C compiler. The model is R syntax, read by R's parser:
#' assignments with `<-` or `=`, `d/dt(x)` for a derivative, `x(0)` for an
#' initial amount, `alag(x)`, `f(x)`, `rate(x)` and `dur(x)` for dosing
#' into `x`, `if`/`else`, comparison and logical operators, `linCmt()`, and
#' the functions `exp`, `log`, `sqrt`, `abs`, `floor`, `ceiling`, `sin`,
#' `cos`, `tan`, `gamma`, `lgamma`, `logit`, `expit`, `probit`,
#' `probitInv`, `min` and `max`. Every name the model reads but never
#' assigns is a parameter, given at solve time.
#'
#' @param model The model: a braced expression, a function whose body is
#'   the model, or its text.
#' @param tier The tier to run on: `"auto"`, `"native"` or `"interp"`; see
#'   [pmx_ode_tier()].
#' @return A model of class `pmx_ode_model`.
#' @section Errors:
#' `pmx_ode_error_compile` for a construct outside the subset, naming it,
#' its position, and the nearest supported alternative.
#' @examples
#' m <- pmx_ode_model({
#'   ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
#'   d/dt(depot) <- -ka * depot
#'   d/dt(central) <- ka * depot - cl / v * central
#'   cp <- central / v
#' })
#' m
#' @export
pmx_ode_model <- function(model, tier = "auto") {
  # A braced model is taken as written; anything else (text, a function, a
  # variable holding either) is evaluated first. Reading `model` before
  # substitute() would evaluate the braces.
  expr <- substitute(model)
  if (!(is.call(expr) && identical(expr[[1]], as.name("{")))) expr <- model
  compiled <- compile_model(expr)
  new_compiled_model(compiled, tier)
}

new_compiled_model <- function(compiled, tier = "auto") {
  m <- new_model(compiled$ir, tier)
  m$states <- compiled$states
  m$params <- compiled$params
  m$captures <- compiled$captures
  m$class <- compiled$class
  m$lincmt <- compiled$lincmt
  m$hoisted <- compiled$hoisted
  m$opts <- compiled$opts
  m
}

#' @export
print.pmx_ode_model <- function(x, ...) {
  cat("<pmx_ode_model> ", x$class, ", ", length(x$states), " compartment",
      if (length(x$states) != 1L) "s", " on the ", pmx_ode_tier(x), " tier\n",
      sep = "")
  if (length(x$states)) cat("  compartments: ", paste(x$states, collapse = ", "), "\n", sep = "")
  if (length(x$params)) cat("  parameters:   ", paste(x$params, collapse = ", "), "\n", sep = "")
  if (length(x$captures)) cat("  outputs:      ", paste(x$captures, collapse = ", "), "\n", sep = "")
  invisible(x)
}
