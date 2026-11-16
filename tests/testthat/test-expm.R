# The matrix exponential (design §4.3): linear systems solved exactly
# between events, against LSODA at tight tolerance, the closed forms, and
# both oracles.

test_that("every linear corpus case: the matrix exponential against LSODA and the oracles", {
  cases <- oracle_cases_of_class("linear")
  expect_gte(length(cases), 25L)
  for (case in cases) {
    c <- oracle_case(case)
    m <- pmx_ode_model(c$rx)
    expm <- pmx_ode_solve(pmx_ode_method(m, "expm"), c$events, c$params)
    ode <- pmx_ode_solve(with_solver(pmx_ode_method(m, "lsoda"), rtol = 1e-11, atol = 1e-13,
                                     ss_tol = 1e-12),
                         c$events, c$params)
    expect_table_equal(expm, ode, rtol = 1e-8, atol = 1e-10, label = paste(case, "lsoda"))
    want <- oracle_read("rxode2", paste0(case, ".csv"))
    expect_table_equal(expm[names(want)], want, label = paste(case, "rxode2"))
  }
})

test_that("on closed-form models the matrix exponential is the closed form", {
  for (case in oracle_cases_of_class("closed-form")) {
    c <- oracle_case(case)
    m <- pmx_ode_model(c$rx)
    closed <- pmx_ode_solve(m, c$events, c$params)
    expm <- pmx_ode_solve(pmx_ode_method(m, "expm"), c$events, c$params)
    expect_table_equal(expm, closed, rtol = 1e-11, atol = 1e-13, label = case)
  }
})

test_that("a constant input and steady state: SS = 1 and SS = 2 with it running", {
  m <- pmx_ode_model("d/dt(central) <- -k * central + kin")
  p <- c(k = 0.1, kin = 2)
  # SS = 1: the steady state of the dose and the input together, against
  # 400 doses.
  ss <- pmx_ode_solve(m, pmx_ode_events(amt = 100, ii = 12, ss = 1L, cmt = "central",
                                         time = c(0, 6)), p)
  long <- pmx_ode_solve(m, pmx_ode_events(amt = 100, ii = 12, addl = 399, cmt = "central",
                                           time = 399 * 12 + c(0, 6)), p)
  expect_equal(ss$central, long$central, tolerance = 1e-10)
  # SS = 2 adds the dose's own steady state to what the input has built.
  base <- 2 / 0.1
  ev <- data.frame(TIME = c(0, 1000, 1000), EVID = c(2, 1, 0), AMT = c(0, 100, 0),
                   II = c(0, 12, 0), SS = c(0, 2, 0), CMT = 1)
  init <- pmx_ode_model("central(0) <- c0\nd/dt(central) <- -k * central + kin")
  got <- pmx_ode_solve(init, ev, c(p, c0 = base))
  dose_ss <- 100 * exp(-0.1 * 12) / (1 - exp(-0.1 * 12))
  expect_equal(got$central, base + dose_ss + 100, tolerance = 1e-12)
  # A constant infusion's steady state: A y + c + r = 0.
  ci <- pmx_ode_solve(m, pmx_ode_events(amt = 0, rate = 5, ss = 1L, cmt = "central",
                                         time = 3), p)
  expect_equal(ci$central, (2 + 5) / 0.1, tolerance = 1e-13)
})

test_that("the state limit: 20 compartments take the matrix exponential, 21 do not", {
  chain <- function(n) paste(c("d/dt(x1) <- -k * x1",
                               sprintf("d/dt(x%d) <- k * x%d - k * x%d", 2:n, 1:(n - 1), 2:n)),
                             collapse = "\n")
  m20 <- pmx_ode_model(chain(20))
  m21 <- pmx_ode_model(chain(21))
  expect_identical(c(m20$class, m21$class), c("linear", "linear"))
  expect_identical(pmx_ode_method(m20), "expm")
  expect_identical(pmx_ode_method(m21), "lsoda")
  expect_error(pmx_ode_method(m21, "expm"), class = "pmx_ode_error_argument")
  # 20 compartments by the matrix exponential: the chain's last amount is
  # the Erlang density, D (k t)^19 e^(-k t) / 19!.
  got <- pmx_ode_solve(m20, pmx_ode_events(amt = 1, cmt = "x1", time = c(5, 20)), c(k = 1))
  expect_equal(got$x20, (c(5, 20))^19 * exp(-c(5, 20)) / factorial(19), tolerance = 1e-10)
})
