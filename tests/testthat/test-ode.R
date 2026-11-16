# The general solvers (design §4.3, D5): LSODA and Dormand-Prince against
# rxode2 and mrgsolve on every corpus case without a closed form, and
# against the closed forms on the cases that have one.

test_that("LSODA matches both oracles on every general and linear case", {
  cases <- c(oracle_cases_of_class("linear"), oracle_cases_of_class("general"))
  expect_gte(length(cases), 60L)
  for (case in cases) {
    c <- oracle_case(case)
    m <- with_solver(pmx_ode_model(c$rx), method = "lsoda", rtol = 1e-10, atol = 1e-12)
    if (sub("__.*", "", case) %in% zupmxode_disagreements) {
      expect_warning(pmx_ode_solve(m, c$events, c$params), class = "pmx_ode_warning_solver")
      next
    }
    got <- suppressWarnings(pmx_ode_solve(m, c$events, c$params))
    for (engine in c("rxode2", "mrgsolve")) {
      path <- oracle_path(engine, paste0(case, ".csv"))
      if (!file.exists(path)) next
      want <- oracle_read(engine, paste0(case, ".csv"))
      expect_table_equal(got[names(want)], want, label = paste(case, engine, "lsoda"))
    }
  }
})

test_that("Dormand-Prince matches both oracles where it can start", {
  cases <- c(oracle_cases_of_class("linear"), oracle_cases_of_class("general"))
  for (case in cases) {
    c <- oracle_case(case)
    m <- with_solver(pmx_ode_model(c$rx), method = "dopri5", rtol = 1e-10, atol = 1e-12)
    if (sub("__.*", "", case) %in% c(dopri_stiff_cases, zupmxode_disagreements)) {
      expect_warning(pmx_ode_solve(m, c$events, c$params), class = "pmx_ode_warning_solver")
      next
    }
    got <- pmx_ode_solve(m, c$events, c$params)
    want <- oracle_read("rxode2", paste0(case, ".csv"))
    expect_table_equal(got[names(want)], want, label = paste(case, "dopri5"))
  }
})

test_that("at the default tolerances both solvers stay near the oracles", {
  for (case in c("pk2_oral__bolus_addl", "pk1_mm__infusion_rate", "irm1__bolus_addl")) {
    c <- oracle_case(case)
    want <- oracle_read("rxode2", paste0(case, ".csv"))
    for (method in names(solver_codes)) {
      got <- pmx_ode_solve(with_solver(pmx_ode_model(c$rx), method = method),
                           c$events, c$params)
      expect_table_equal(got[names(want)], want, rtol = 1e-4, atol = 1e-6,
                         label = paste(case, method))
    }
  }
})

test_that("every closed-form case: the closed form against LSODA on its equations", {
  for (case in oracle_cases_of_class("closed-form")) {
    c <- oracle_case(case)
    m <- pmx_ode_model(c$rx)
    closed <- pmx_ode_solve(m, c$events, c$params)
    # ss_tol too: repeated dosing stops within it of the steady state.
    ode <- pmx_ode_solve(with_solver(m, closed = FALSE, rtol = 1e-10, atol = 1e-12,
                                     ss_tol = 1e-12),
                         c$events, c$params)
    expect_table_equal(ode, closed, rtol = 1e-7, atol = 1e-9, label = case)
  }
})

test_that("dense output gives what stepping to each time gives", {
  m <- with_solver(pmx_ode_model("d/dt(depot) <- -ka * depot\nd/dt(central) <- ka * depot - k * central"),
                   method = "dopri5", rtol = 1e-10, atol = 1e-12)
  times <- seq(0.25, 24, by = 0.25)
  p <- c(ka = 0.8, k = 0.1)
  dense <- pmx_ode_solve(m, pmx_ode_events(amt = 100, cmt = "depot", time = times), p)
  direct <- vapply(times, function(t) {
    pmx_ode_solve(m, pmx_ode_events(amt = 100, cmt = "depot", time = t), p)$central
  }, 0)
  exact <- 100 * 0.8 / 0.7 * (exp(-0.1 * times) - exp(-0.8 * times))
  expect_equal(dense$central, direct, tolerance = 1e-9)
  expect_equal(dense$central, exact, tolerance = 1e-9)
})

test_that("a compartment switched off is emptied and held until a dose", {
  m <- pmx_ode_model("d/dt(central) <- -k * central")
  ev <- data.frame(TIME = c(0, 2, 3, 4, 5, 6), EVID = c(1, 2, 0, 1, 0, 0),
                   AMT = c(100, 0, 0, 10, 0, 0), CMT = c(1, -1, 1, 1, 1, 1))
  got <- pmx_ode_solve(with_solver(m, rtol = 1e-10, atol = 1e-12), ev, c(k = 0.1))
  expect_equal(got$central, c(0, 10 * exp(-0.1), 10 * exp(-0.2)), tolerance = 1e-9)
})

test_that("each listed disagreement is R's own arithmetic, and has its note", {
  # The r0rmaxcrmax models: IC50 is exactly zero as R computes it.
  c <- oracle_case("lib_indirect_1cpt_inhi_kin_r0rmaxcrmax")
  p <- exp(c$params)
  imax <- (p[["lr0"]] - p[["lrmax"]]) / p[["lr0"]]
  ic50 <- p[["lcrmax"]] * (p[["lrmax"]] - (1 - imax) * p[["lr0"]]) / (p[["lr0"]] - p[["lrmax"]])
  expect_identical(ic50, 0)
  notes <- readLines(oracle_path("disagreements.md"))
  headed <- sub("^#+ +", "", grep("^#{2,} ", notes, value = TRUE))
  expect_true(all(zupmxode_disagreements %in% headed))
})

test_that("a compartment switched off takes nothing in while off", {
  # Drug flows from the depot into central; with central off, what the
  # depot loses goes nowhere (PREDPP), and central stays at 0.
  m <- pmx_ode_model("d/dt(depot) <- -ka * depot\nd/dt(central) <- ka * depot - k * central")
  ev <- data.frame(TIME = c(0, 1, 2, 3), EVID = c(1, 2, 0, 0), AMT = c(100, 0, 0, 0),
                   CMT = c(1, -2, 1, 1))
  got <- pmx_ode_solve(with_solver(m, rtol = 1e-10, atol = 1e-12), ev, c(ka = 0.5, k = 0.1))
  expect_identical(got$central, c(0, 0))
  expect_equal(got$depot, 100 * exp(-0.5 * c(2, 3)), tolerance = 1e-9)
})
