# The tiers (design §4.7, §6.3): every hand-built model's blocks give the
# same doubles on the native tier, the interpreter tier and the reference
# evaluator, and those doubles are R's own arithmetic on the same inputs.

test_that("a one-compartment model gives R's answer on every tier", {
  hb <- hb_one_cpt()
  for (tier in available_tiers()) {
    m <- new_model(hb$ir, tier)
    expect_identical(pmx_ode_tier(m), tier)
    for (case in hb$cases) {
      got <- hb$run(hb_eval_model(m), case$par, case$y)
      ref <- hb$run(hb_eval_ref(hb$ir), case$par, case$y)
      expect_identical(got, ref, label = tier)
      expect_equal(got, hb$expect(case$par, case$y), tolerance = 0,
                   label = tier)
    }
  }
})

test_that("branches and comparisons follow C on every tier", {
  hb <- hb_branches()
  for (tier in available_tiers()) {
    m <- new_model(hb$ir, tier)
    for (par in hb$cases) {
      got <- model_eval(m, "rhs", t = 7, par = par)$out
      expect_identical(got, ir_eval(hb$ir, "rhs", t = 7, par = par)$out)
      expect_identical(got, hb$expect(par, 7),
                       label = paste(tier, paste(par, collapse = ", ")))
    }
  }
})

test_that("constants, two-argument calls, signed zero and NaN agree", {
  hb <- hb_ops()
  for (tier in available_tiers()) {
    m <- new_model(hb$ir, tier)
    for (par in hb$cases) {
      got <- model_eval(m, "rhs", par = par)
      expect_identical(got, ir_eval(hb$ir, "rhs", par = par))
      expect_identical(got$out, hb$expect(par),
                       label = paste(tier, paste(par, collapse = ", ")))
      # neg(-0) is +0: identical() cannot see the sign, 1 / x can.
      expect_identical(1 / got$out[4], Inf)
    }
  }
})

test_that("the tiers agree with each other bit for bit", {
  tiers <- available_tiers()
  skip_if(length(tiers) < 2L, "only one tier runs here")
  hb <- hb_ops()
  a <- new_model(hb$ir, "native")
  b <- new_model(hb$ir, "interp")
  set.seed(1)
  for (i in 1:200) {
    par <- c(stats::rnorm(1, sd = 10), stats::rnorm(1, sd = 3))
    x <- model_eval(a, "rhs", par = par)$out
    y <- model_eval(b, "rhs", par = par)$out
    expect_identical(writeBin(x, raw()), writeBin(y, raw()))
  }
})

test_that("empty blocks run and leave their outputs at zero", {
  ir <- ir_asm(list(), ny = 1, npar = 1, nout = c(post = 2))
  for (tier in available_tiers()) {
    m <- new_model(ir, tier)
    expect_identical(model_eval(m, "post", y = 1, par = 1)$out, c(0, 0))
  }
})

test_that("pmx_ode_tier() refuses what is not a model or a tier", {
  m <- new_model(hb_ops()$ir)
  expect_error(pmx_ode_tier(list()), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_tier(m, "gpu"), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_tier(m, NA_character_), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_tier(m, c("native", "interp")),
               class = "pmx_ode_error_argument")
  expect_error(pmx_ode_tier(m, 1), class = "pmx_ode_error_argument")
  expect_true(pmx_ode_tier(m) %in% c("native", "interp"))
})

test_that("a session where no tier runs is a classed error (D15)", {
  local_mocked_bindings(zupmxode_tiers = function() {
    c(native = FALSE, interp = FALSE)
  })
  expect_error(new_model(hb_ops()$ir), class = "pmx_ode_error_no_tier")
  expect_error(resolve_tier("interp"), class = "pmx_ode_error_no_tier")
})

test_that("auto takes the interpreter where native code does not run", {
  skip_if_not("interp" %in% available_tiers())
  local_mocked_bindings(zupmxode_tiers = function() {
    c(native = FALSE, interp = TRUE)
  })
  expect_identical(resolve_tier("auto"), "interp")
  expect_error(resolve_tier("native"), class = "pmx_ode_error_no_tier")
})
