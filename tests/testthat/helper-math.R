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
# lgamma's is wider because its result goes through zero near x = -2.457,
# where a small absolute error is a large relative one (10 units at -2.5
# with macOS's library).
math_ulps <- c(
  exp = 1, log = 1, sqrt = 0.5, abs = 0, floor = 0, ceiling = 0, sin = 1,
  cos = 1, tan = 2, gamma = 8, lgamma = 16, logit = 4, expit = 2, probit = 4,
  probitInv = 4, pow = 2
)

reference <- function() {
  r <- utils::read.csv(test_path("fixtures", "math-reference.csv"),
                       stringsAsFactors = FALSE)
  r$x <- as.numeric(r$x)
  r$y <- as.numeric(r$y)
  r$ref <- as.numeric(r$ref)
  r
}

