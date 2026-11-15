# Shared by test-math-table.R: its oracles and bounds (the suite runs
# shuffled, so nothing a test uses may live at a test file's top level).

math_r <- list(
  exp = exp, log = log, sqrt = sqrt, abs = abs, floor = floor,
  ceiling = ceiling, sin = sin, cos = cos, tan = tan, gamma = gamma,
  lgamma = lgamma, logit = stats::qlogis, expit = stats::plogis,
  probit = stats::qnorm, probitInv = stats::pnorm
)

math_special <- c(NaN, -Inf, Inf, 0, -0, 1, -1, -3, -10, 2)

# Units in the last place of `ref` between got and ref.
ulps <- function(got, ref) {
  spacing <- 2^pmax(floor(log2(abs(ref))) - 52, -1074)
  abs(got - ref) / spacing
}

# The largest error allowed against the 200-bit value, per function: 0.5
# for what the C library rounds correctly, a few units for the rest.
# gamma's and lgamma's are wider: mingw-w64's tgamma() reaches 10 units,
# and lgamma's result goes through zero near x = -2.457, where a small
# absolute error is a large relative one (10 units at -2.5 on macOS).
math_ulps <- c(
  exp = 1, log = 1, sqrt = 0.5, abs = 0, floor = 0, ceiling = 0, sin = 1,
  cos = 1, tan = 2, gamma = 16, lgamma = 16, logit = 4, expit = 2, probit = 4,
  probitInv = 4, pow = 2
)

# The 16 hexadecimal digits of IEEE bits, as tools/make-math-reference.R
# writes them, back to doubles.
bits_double <- function(h) {
  bytes <- as.raw(strtoi(substring(rep(h, each = 8), seq(1, 15, 2), seq(2, 16, 2)), 16L))
  readBin(bytes, "double", n = length(h), size = 8, endian = "big")
}

reference <- function() {
  r <- utils::read.csv(test_path("fixtures", "math-reference.csv"),
                       stringsAsFactors = FALSE, colClasses = "character")
  r$x <- bits_double(r$x)
  r$y <- bits_double(r$y)
  r$ref <- bits_double(r$ref)
  r
}

