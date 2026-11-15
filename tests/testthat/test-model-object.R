# The compiled model owns a MIR context through an external pointer
# (design §4.7): freed by its finalizer once, never while in use, and not
# stranded when an error unwinds past it.

test_that("dropping a model frees its context, once", {
  gc()
  live <- .Call(zpx_live_models)
  m <- new_model(hb_ops()$ir)
  expect_identical(.Call(zpx_live_models), live + 1L)
  rm(m)
  gc()
  expect_identical(.Call(zpx_live_models), live)
  gc()
  expect_identical(.Call(zpx_live_models), live)
})

test_that("copies share one context, freed when the last goes", {
  gc()
  live <- .Call(zpx_live_models)
  m <- new_model(hb_ops()$ir)
  m2 <- m
  rm(m)
  gc()
  expect_identical(.Call(zpx_live_models), live + 1L)
  expect_identical(model_eval(m2, "rhs", par = c(2, 3))$out[1], 8)
  rm(m2)
  gc()
  expect_identical(.Call(zpx_live_models), live)
})

test_that("pmx_ode_tier() makes a new model and leaves the old one working", {
  hb <- hb_ops()
  m <- new_model(hb$ir)
  for (tier in available_tiers()) {
    m2 <- pmx_ode_tier(m, tier)
    expect_identical(pmx_ode_tier(m2), tier)
    expect_identical(model_eval(m, "rhs", par = c(2, 3))$out,
                     model_eval(m2, "rhs", par = c(2, 3))$out)
  }
})

test_that("an error after lowering does not strand the context", {
  gc()
  live <- .Call(zpx_live_models)
  ir <- hb_ops()$ir
  expect_error(.Call(zpx_test_unwind, ir$code, ir$consts,
                     tier_codes[[available_tiers()[1]]]),
               "test hook")
  gc()
  expect_identical(.Call(zpx_live_models), live)
})

test_that("a forked child runs the parent's model", {
  skip_on_os("windows")
  m <- new_model(hb_one_cpt()$ir)
  par <- c(log(0.8), log(2), log(20))
  parent <- model_eval(m, "rhs", par = par, y = c(100, 5),
                       work = exp(par))$out
  child <- parallel::mccollect(parallel::mcparallel(
    model_eval(m, "rhs", par = par, y = c(100, 5), work = exp(par))$out
  ))[[1]]
  expect_identical(child, parent)
})

test_that("a block refuses arrays of the wrong length", {
  m <- new_model(hb_one_cpt()$ir)
  expect_error(model_eval(m, "rhs", par = 1, y = c(1, 2)), "`par`")
  expect_error(model_eval(m, "rhs", par = c(1, 2, 3), y = 1), "`y`")
  expect_error(model_eval(m, "rhs", par = c(1, 2, 3), y = c(1, 2), work = 1),
               "`work`")
})
