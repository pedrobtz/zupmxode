# Per-subject evaluation (design §4.6, D12): pmx_ode_subject_eval() and
# zupmxode_subject_eval() give pmx_ode_solve()'s predictions, bit for bit,
# and the error model's variance.

test_that("every corpus model's predictions are pmx_ode_solve()'s, and var its formula", {
  skip_under_gctorture()
  man <- oracle_read("manifest.csv")
  evaluated <- 0L
  for (case in unique(man$case[man$status == "ok"])) {
    c <- oracle_case(case)
    m <- tryCatch(pmx_ode_model(c$rx), pmx_ode_error_compile = function(e) NULL)
    if (is.null(m)) next
    want <- with_warnings(pmx_ode_solve(m, c$events, c$params))
    subs <- pmx_ode_subjects(m, c$events)
    th <- c$params[attr(subs, "theta")]
    et <- c$params[attr(subs, "eta")]
    for (out in c(m$states, m$captures)) {
      sigma <- stats::setNames(list(add(0.5) + prop(0.1)), out)
      got <- with_warnings(pmx_ode_subject_eval(m, subs[[1]], th, et, sigma))
      f <- want$value[[out]]
      expect_identical(got$value$pred, f, label = paste(case, out))
      expect_identical(got$value$var, 0.5^2 + 0.1^2 * f^2, label = paste(case, out))
      expect_identical(length(got$warnings), length(want$warnings), label = paste(case, out))
      evaluated <- evaluated + 1L
    }
  }
  expect_gt(evaluated, 150L)
})

test_that("each error model's variance is its formula", {
  m <- sim_model()
  subs <- pmx_ode_subjects(m, sim_events())
  eta <- c(eta.cl = 0.1, eta.v = -0.2)
  f <- pmx_ode_solve(m, sim_events(), c(sim_theta(), eta))$cp
  forms <- list(
    list(add(0.7), rep(0.7^2, 4)),
    list(prop(0.2), 0.2^2 * f^2),
    list(add(0.7) + prop(0.2), 0.7^2 + 0.2^2 * f^2),
    list(lnorm(0.3), rep(0.3^2, 4))
  )
  for (k in forms) {
    r <- pmx_ode_subject_eval(m, subs[[1]], sim_theta(), eta, list(cp = k[[1]]))
    expect_identical(r$pred, f)
    expect_identical(r$var, k[[2]])
    expect_identical(r$time, c(1, 4, 12, 24))
  }
})

test_that("each subject of a table is its own, with its covariates", {
  m <- thread_models()$closed
  ev <- population_events(5)
  p <- c(tcl = 2, tv = 20, tka = 0.8)
  want <- pmx_ode_solve(m, ev, p)
  subs <- pmx_ode_subjects(m, ev)
  expect_identical(attr(subs, "theta"), c("tcl", "tv", "tka"))
  expect_identical(vapply(subs, `[[`, 0, "id"), as.double(1:5))
  for (i in 1:5) {
    r <- pmx_ode_subject_eval(m, subs[[i]], p[c(3, 1, 2)], sigma = list(cp = add(1)))
    expect_identical(r$pred, want$cp[want$id == i])
  }
})

test_that("the same call from many threads on different subjects gives the serial results", {
  skip_under_gctorture()
  ev <- population_events(40)
  p <- c(tcl = 2, tv = 20, tka = 0.8, ka = 0.8, k = 0.1, vmax = 10, km = 2)
  n <- test_threads()
  for (name in names(thread_models())) {
    for (tier in available_tiers()) {
      m <- pmx_ode_tier(thread_models()[[name]], tier)
      subs <- pmx_ode_subjects(m, ev)
      th <- p[attr(subs, "theta")]
      sigma <- list(central = prop(0.1))
      one <- test_eval_threads(m, subs, th, double(), sigma, 1)
      many <- test_eval_threads(m, subs, th, double(), sigma, n)
      expect_identical(many, one, label = paste(name, tier))
      expect_true(all(one$status == 0L), label = paste(name, tier))
      serial <- lapply(subs, function(s) {
        pmx_ode_subject_eval(m, s, th, sigma = sigma)$pred
      })
      expect_identical(one$pred, serial, label = paste(name, tier))
    }
  }
})

test_that("tol and the method reach the solver", {
  m <- thread_models()$general
  ev <- population_events(1)
  p <- c(ka = 0.8, vmax = 10, km = 2)
  subs <- pmx_ode_subjects(m, ev)
  loose <- pmx_ode_subject_eval(m, subs[[1]], p, sigma = list(central = add(1)), tol = c(1e-3, 1e-3))
  want <- pmx_ode_solve(with_solver(m, rtol = 1e-3, atol = 1e-3), ev, p)
  expect_identical(loose$pred, want$central)
  dp <- pmx_ode_method(m, "dopri5")
  subs <- pmx_ode_subjects(dp, ev)
  r <- pmx_ode_subject_eval(dp, subs[[1]], p, sigma = list(central = add(1)))
  expect_identical(r$pred, pmx_ode_solve(dp, ev, p)$central)
})

test_that("event and solver failures are raised and warned as pmx_ode_solve() does", {
  m <- pmx_ode_method(pmx_ode_model("d/dt(central) <- -k * central^2\ncp <- central"), "dopri5")
  ev <- pmx_ode_events(amt = 1, cmt = 1, time = c(0.01, 1, 2))
  subs <- pmx_ode_subjects(m, ev)
  expect_warning(r <- pmx_ode_subject_eval(m, subs[[1]], c(k = -50), sigma = list(cp = add(1))),
                 class = "pmx_ode_warning_solver")
  expect_identical(is.na(r$pred), c(FALSE, TRUE, TRUE))
  expect_identical(is.na(r$var), c(FALSE, TRUE, TRUE))
  # A lag the model makes negative breaks an event rule at these parameters.
  lag <- pmx_ode_model("alag(central) <- lg\nd/dt(central) <- -central\ncp <- central")
  ev <- pmx_ode_events(amt = 1, cmt = 1, time = c(1, 2))
  subs <- pmx_ode_subjects(lag, ev)
  expect_error(pmx_ode_subject_eval(lag, subs[[1]], c(lg = -1), sigma = list(cp = add(1))),
               class = "pmx_ode_error_events")
  # A closed form with a rate constant that is not positive.
  cf <- pmx_ode_model("cl <- tcl\nv <- tv\ncp <- linCmt()")
  ev <- pmx_ode_events(amt = 1, cmt = "central", time = c(1, 2))
  subs <- pmx_ode_subjects(cf, ev)
  expect_error(pmx_ode_subject_eval(cf, subs[[1]], c(tcl = -1, tv = 1), sigma = list(cp = add(1))),
               class = "pmx_ode_error_solver")
})
