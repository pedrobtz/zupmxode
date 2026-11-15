# The math table (design §4.1). Two oracles: at poles, infinities, NaN and
# signed zero the table follows R's own functions exactly, since those are
# the values a model written for rxode2 sees; for finite arguments it is
# held to the true value, computed with Rmpfr at 200 bits
# (tools/make-math-reference.R), because R itself is off there (gamma(100.5)
# by 228 units in the last place, lgamma(1 - 2^-53) by 1e16).

test_that("each function is within its error of the 200-bit value", {
  r <- reference()
  for (f in unique(r$f)) {
    s <- r[r$f == f, ]
    got <- .Call(zpx_math_r, ir_math[[f]], s$x,
                 if (ir_math[[f]] >= 32L) s$y else double())
    err <- ulps(got, s$ref)
    expect_lte(max(err), math_ulps[[f]], label = paste(f, "ulps"))
  }
  # Every function of the table has reference values.
  expect_setequal(unique(r$f), setdiff(names(ir_math), c("min", "max")))
})

test_that("at poles, infinities, NaN and zeros each function is R's", {
  for (f in names(math_r)) {
    want <- suppressWarnings(math_r[[f]](math_special))
    got <- .Call(zpx_math_r, ir_math[[f]], math_special, double())
    # Where R's answer is special: non-finite, zero, or at a pole of gamma
    # (the negative integers). A finite result is the reference test's.
    # identical() takes -0 and 0 as equal: the sign of a zero result is
    # C's, as rxode2's compiled code gives it (R's sin(-0) is +0, C's -0).
    special <- !is.finite(want) | want == 0
    expect_identical(got[special], want[special], label = f)
  }
})

test_that("pow follows R's ^ wherever an operand or the result is special", {
  g <- expand.grid(a = c(math_special, 0.5, -0.5, 1e300),
                   b = c(math_special, 0.5, -0.5, 1 / 3, -2, 3, 1e300))
  got <- .Call(zpx_math_r, ir_math[["pow"]], g$a, g$b)
  # R warns of lost accuracy for (-Inf)^1e300; the answer is still R's.
  want <- suppressWarnings(g$a^g$b)
  special <- !is.finite(g$a) | !is.finite(g$b) | !is.finite(want) |
    want == 0 | g$a == 0 | g$a == 1 | g$b == 0 | g$b == 2
  expect_identical(got[special], want[special])
  zero <- special & !is.na(want) & want == 0
  expect_identical(sign(1 / got[zero]), sign(1 / want[zero]))
})

test_that("min and max are R's, NaN and signed zero included", {
  g <- expand.grid(a = c(math_special, 0.5), b = c(math_special, 0.5))
  expect_identical(.Call(zpx_math_r, ir_math[["min"]], g$a, g$b),
                   mapply(min, g$a, g$b))
  expect_identical(.Call(zpx_math_r, ir_math[["max"]], g$a, g$b),
                   mapply(max, g$a, g$b))
})

test_that("every tier calls the table as the table computes it", {
  x <- c(math_special, reference()$x[1:40])
  for (tier in available_tiers()) {
    for (f in names(ir_math)) {
      two <- ir_math[[f]] >= 32L
      ir <- ir_asm(list(rhs = paste0(
        "loadp r0 0\nloadp r1 1\n",
        if (two) paste("call2 r2", f, "r0 r1") else paste("call1 r2", f, "r0"),
        "\nstoreo 0 r2"
      )), npar = 2, nout = c(rhs = 1))
      m <- new_model(ir, tier)
      b <- if (two) 2.5 else 0
      got <- vapply(x, function(v) model_eval(m, "rhs", par = c(v, b))$out, 0)
      want <- .Call(zpx_math_r, ir_math[[f]], x,
                    if (two) rep(b, length(x)) else double())
      expect_identical(writeBin(got, raw()), writeBin(want, raw()),
                       label = paste(tier, f))
    }
  }
})
