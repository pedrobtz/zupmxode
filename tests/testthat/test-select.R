# Solver selection (design §4.3, D6): each class gets its most specific
# solver, pmx_ode_method() overrides it within what the class allows, and
# pmx_ode_info() says what was chosen and why.

test_that("each class gets its documented solver", {
  closed <- pmx_ode_model("cl <- tcl\nv <- tv\ncp <- linCmt()")
  linear <- pmx_ode_model("d/dt(x) <- -k * x")
  general <- pmx_ode_model("d/dt(x) <- -vmax * x / (km + x)")
  expect_identical(pmx_ode_method(closed), "closed")
  expect_identical(pmx_ode_method(linear), "expm")
  expect_identical(pmx_ode_method(general), "lsoda")
  expect_match(pmx_ode_info(closed)$reason, "linCmt")
  expect_match(pmx_ode_info(general)$reason, "not linear")
})

test_that("pmx_ode_method() overrides within what the class allows", {
  linear <- pmx_ode_model("d/dt(x) <- -k * x")
  general <- pmx_ode_model("d/dt(x) <- -x * x")
  for (method in c("expm", "lsoda", "dopri5", "auto")) {
    m <- pmx_ode_method(linear, method)
    expect_identical(pmx_ode_method(m), if (method == "auto") "expm" else method)
    got <- pmx_ode_solve(m, pmx_ode_events(amt = 1, time = 2), c(k = 0.5))
    expect_equal(got$x, exp(-1), tolerance = 1e-5, label = method)
  }
  expect_identical(pmx_ode_info(pmx_ode_method(linear, "dopri5"))$reason,
                   "set with pmx_ode_method()")
  expect_error(pmx_ode_method(linear, "closed"), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_method(general, "expm"), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_method(linear, "rk4"), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_method(linear, NA_character_), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_method(list()), class = "pmx_ode_error_argument")
  expect_error(pmx_ode_info(list()), class = "pmx_ode_error_argument")
})

test_that("a compartment switched off takes auto to LSODA, and refuses the others", {
  closed <- pmx_ode_model("cl <- tcl\nv <- tv\ncp <- linCmt()")
  off <- data.frame(TIME = c(0, 1, 2), EVID = c(1, 2, 0), AMT = c(10, 0, 0), CMT = c(1, -1, 1))
  got <- pmx_ode_solve(closed, off, c(tcl = 1, tv = 1))
  expect_identical(got$central, 0)
  expect_error(pmx_ode_solve(pmx_ode_method(closed, "expm"), off, c(tcl = 1, tv = 1)),
               class = "pmx_ode_error_solver")
})

test_that("pmx_ode_info() reports the model, its solver and its tier", {
  m <- pmx_ode_model({
    ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / v * central
    cp <- central / v
  })
  info <- pmx_ode_info(m)
  expect_s3_class(info, "pmx_ode_info")
  expect_identical(info$compartments, c("depot", "central"))
  expect_identical(info$parameters, c("tka", "tcl", "tv"))
  expect_identical(info$solver, "expm")
  expect_true(info$tier %in% c("native", "interp"))
  expect_true(is.numeric(info$lowering_seconds) && info$lowering_seconds >= 0)
  expect_identical(names(info$instructions), ir_blocks)
  expect_output(print(info), "solver:       expm")
})
