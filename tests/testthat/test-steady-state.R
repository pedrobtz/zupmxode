# Steady state for models with no closed form (design §4.2): repeated
# dosing until the state before a dose stops moving, with a cap.

test_that("a general model's steady state is where repeated dosing goes", {
  m <- with_solver(pmx_ode_model("cp <- central / v\nd/dt(central) <- -vmax * cp / (km + cp)"),
                   rtol = 1e-10, atol = 1e-12)
  p <- c(v = 20, vmax = 10, km = 2)
  ss <- pmx_ode_solve(m, pmx_ode_events(amt = 100, ii = 24, ss = 1L, cmt = "central",
                                         time = c(0, 6, 23.9)), p)
  long <- pmx_ode_solve(m, pmx_ode_events(amt = 100, ii = 24, addl = 299, cmt = "central",
                                           time = 299 * 24 + c(0, 6, 23.9)), p)
  expect_equal(ss$central, long$central, tolerance = 1e-7)
})

test_that("the cap: converging on the last interval allowed is no warning", {
  # LSODA: a linear model would take the matrix exponential, whose steady
  # state is solved directly, with no intervals to count.
  m <- pmx_ode_method(pmx_ode_model("d/dt(central) <- -k * central"), "lsoda")
  ev <- pmx_ode_events(amt = 100, ii = 12, ss = 1L, cmt = "central", time = 1)
  warns <- function(cap) {
    w <- FALSE
    withCallingHandlers(
      pmx_ode_solve(with_solver(m, ss_cap = cap), ev, c(k = 0.05)),
      pmx_ode_warning_steady_state = function(cnd) {
        w <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    w
  }
  # The fewest intervals that converge, found from below.
  cap <- 1
  while (warns(cap) && cap < 200) cap <- cap + 1
  expect_lt(cap, 200)
  expect_false(warns(cap))
  expect_true(warns(cap - 1))
  # The steady state there is the closed form's, within ss_tol's reach.
  closed <- pmx_ode_solve(pmx_ode_model("k <- tk\nv <- 1\ncp <- linCmt()"), ev, c(tk = 0.05))
  got <- pmx_ode_solve(with_solver(m, ss_cap = cap, rtol = 1e-10, atol = 1e-12), ev, c(k = 0.05))
  expect_equal(got$central, closed$central, tolerance = 1e-6)
})
