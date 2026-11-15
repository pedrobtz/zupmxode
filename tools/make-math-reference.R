#!/usr/bin/env Rscript
# Writes tests/testthat/fixtures/math-reference.csv: the math table's
# functions (design §4.1) on a grid of arguments, computed with Rmpfr at 200
# bits and rounded to the nearest double, the oracle test-math-table.R
# holds the table to. R's own functions are not that oracle: R's gamma()
# and lgamma() are off by up to 228 and 1e16 units in the last place on
# this grid (gamma(100.5), lgamma(1 - 2^-53)).
#
#   Rscript tools/make-math-reference.R      # needs Rmpfr
#
# Only finite arguments with a finite reference are written; poles,
# infinities and NaN follow R's rules, which the test checks against R.
suppressPackageStartupMessages(library(Rmpfr))
stopifnot(file.exists("DESCRIPTION"))
bits <- 200
grid <- c(0, 1, -1, 0.5, -0.5, 1e-300, -1e-300, 2.5, -2.5, 3, -3, 10, -10,
          170.5, 171.7, 1e10, -1e10, 1e300, .Machine$double.xmin, 1 - 2^-53,
          2^-1074, pi, 100.5, 7e-17, -37.5, 37.5, 1e-5, 0.999, 20.25, -20.25,
          5.5, -0.001, 0.1, 0.9, 0.25, 0.75, 1e-10, 1 - 1e-10,
          # either side of zpx_lgamma()'s switch to Stirling's series
          9.999999e9, 1.000001e10)
grid <- unique(grid)
M <- function(x) mpfr(x, bits)

# probit by Newton's method on pnorm() at 200 bits, from qnorm()'s start.
probit_ref <- function(p) {
  z <- M(stats::qnorm(p))
  P <- M(p)
  for (i in 1:8) z <- z - (pnorm(z) - P) / dnorm(z)
  z
}

ref1 <- list(
  exp = function(X) exp(X), log = function(X) log(X),
  sqrt = function(X) sqrt(X), abs = function(X) abs(X),
  floor = function(X) floor(X), ceiling = function(X) ceiling(X),
  sin = function(X) sin(X), cos = function(X) cos(X), tan = function(X) tan(X),
  gamma = function(X) gamma(X), lgamma = function(X) lgamma(X),
  logit = function(X) log(X / (1 - X)),
  expit = function(X) 1 / (1 + exp(-X)),
  probit = NULL,
  probitInv = function(X) pnorm(X)
)

# A double as the 16 hexadecimal digits of its IEEE bits, most significant
# first; tests/testthat/helper-math.R reads them back.
double_bits <- function(x) {
  b <- matrix(writeBin(as.double(x), raw(), endian = "big"), nrow = 8)
  apply(b, 2, function(col) paste(format(col), collapse = ""))
}

rows <- list()
add <- function(f, x, y, ref) {
  keep <- is.finite(ref) & is.finite(x) & is.finite(y)
  if (any(keep)) {
    rows[[length(rows) + 1L]] <<- data.frame(
      # IEEE bit patterns: R's decimal reader is not exact where long
      # double is double (macOS arm64), and it reads hexadecimal
      # subnormals as 0, so neither text form round-trips everywhere.
      f = f, x = double_bits(x[keep]), y = double_bits(y[keep]),
      ref = double_bits(ref[keep])
    )
  }
}
for (f in names(ref1)) {
  x <- grid
  if (f %in% c("logit", "probit")) x <- x[x > 0 & x < 1]
  if (f %in% c("log", "sqrt")) x <- x[x >= 0]
  ref <- if (f == "probit") {
    vapply(x, function(p) asNumeric(probit_ref(p)), 0)
  } else {
    suppressWarnings(asNumeric(ref1[[f]](M(x))))
  }
  add(f, x, rep(0, length(x)), ref)
}
g <- expand.grid(x = grid[grid > 0], y = c(0.5, 2, 3, -1, 1 / 3, 10.5, -2.5))
add("pow", g$x, g$y, suppressWarnings(asNumeric(M(g$x)^M(g$y))))

out <- file.path("tests", "testthat", "fixtures", "math-reference.csv")
utils::write.csv(do.call(rbind, rows), out, row.names = FALSE)
cat("wrote", out, "\n")
