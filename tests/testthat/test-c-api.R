# The C API's refusals (design §4.6, §8): a bad call returns its status
# code and never raises an R error; and the R wrapper's own refusals.

test_that("bad calls return their status codes, never an R error", {
  m <- thread_models()$general
  big <- pmx_ode_subjects(m, population_events(1))
  small <- pmx_ode_subjects(m, population_events(1, times = 1))
  want <- c(ok = 0L, version = 16L, version0 = 16L, buffer = 18L, arena = 19L, null = 17L,
            null_theta = 17L, output = 17L, output_neg = 17L, error_model = 17L, sd = 17L,
            sd_nan = 17L, method = 17L, method_range = 17L)
  for (what in names(want)) {
    expect_no_error(st <- test_c_api(m, big[[1]], small[[1]], what))
    expect_identical(st, want[[what]], label = what)
  }
  # A subject of another model with the same dimensions, in this model's arena.
  other <- pmx_ode_model("d/dt(depot) <- -ka * depot\nd/dt(central) <- ka * depot - vmax * central / (km + central) * wt / 71")
  foreign <- pmx_ode_subjects(other, population_events(1, times = 1))
  expect_identical(test_c_api(m, big[[1]], foreign[[1]], "other_model"), 17L)
  # Each code has its own description.
  msgs <- vapply(c(0:2, 16:21), function(s) .Call(zpx_status_message_r, s), "")
  expect_false(anyDuplicated(msgs) > 0)
  expect_identical(.Call(zpx_status_message_r, 99L), "unknown status")
})

test_that("a subject prepared for another model is refused", {
  a <- thread_models()$general
  b <- thread_models()$linear
  subs <- pmx_ode_subjects(a, population_events(1))
  expect_error(pmx_ode_subject_eval(b, subs[[1]], c(ka = 1, k = 1), sigma = list(central = add(1))),
               class = "pmx_ode_error_argument")
  # Same compartments, parameters and outputs, other equations.
  c2 <- pmx_ode_model("d/dt(depot) <- -ka * depot\nd/dt(central) <- ka * depot - vmax * central / (km + central) * wt / 71")
  expect_error(pmx_ode_subject_eval(c2, subs[[1]], c(ka = 1, vmax = 1, km = 1),
                                    sigma = list(central = add(1))),
               class = "pmx_ode_error_argument")
  # The same model on the other tier is the same model.
  for (tier in available_tiers()) {
    r <- pmx_ode_subject_eval(pmx_ode_tier(a, tier), subs[[1]], c(ka = 0.8, vmax = 10, km = 2),
                              sigma = list(central = add(1)))
    expect_identical(r$pred, pmx_ode_solve(a, population_events(1), c(ka = 0.8, vmax = 10, km = 2))$central)
  }
})

test_that("a subject saved and loaded again must be prepared again", {
  m <- thread_models()$general
  subs <- pmx_ode_subjects(m, population_events(1))
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(subs, path)
  back <- readRDS(path)
  expect_error(pmx_ode_subject_eval(m, back[[1]], c(ka = 0.8, vmax = 10, km = 2),
                                    sigma = list(central = add(1))),
               class = "pmx_ode_error_argument")
  # The original still works, twice, with the same bits.
  r1 <- pmx_ode_subject_eval(m, subs[[1]], c(ka = 0.8, vmax = 10, km = 2), sigma = list(central = add(1)))
  r2 <- pmx_ode_subject_eval(m, subs[[1]], c(ka = 0.8, vmax = 10, km = 2), sigma = list(central = add(1)))
  expect_identical(r1, r2)
})

test_that("a subject outlives the model object it was prepared with", {
  subs <- local({
    m <- thread_models()$general
    pmx_ode_subjects(m, population_events(1))
  })
  gc()
  m <- thread_models()$general
  r <- pmx_ode_subject_eval(m, subs[[1]], c(ka = 0.8, vmax = 10, km = 2), sigma = list(central = add(1)))
  expect_identical(r$pred, pmx_ode_solve(m, population_events(1), c(ka = 0.8, vmax = 10, km = 2))$central)
})

test_that("pmx_ode_subjects() and pmx_ode_subject_eval() refuse bad arguments", {
  m <- sim_model()
  ev <- sim_events()
  bad <- function(expr) expect_error(expr, class = "pmx_ode_error_argument")
  bad(pmx_ode_subjects("m", ev))
  bad(pmx_ode_subjects(m, ev, eta = "nope"))
  bad(pmx_ode_subjects(m, ev, eta = c("eta.cl", "eta.cl")))
  bad(pmx_ode_subjects(m, ev, eta = NA_character_))
  ev2 <- ev
  ev2$eta.cl <- 0
  bad(pmx_ode_subjects(m, ev2, eta = "eta.cl"))
  subs <- pmx_ode_subjects(m, ev)
  s <- subs[[1]]
  eta <- c(eta.cl = 0, eta.v = 0)
  ok <- function(...) {
    args <- utils::modifyList(list(model = m, subject = s, theta = sim_theta(), eta = eta,
                                   sigma = list(cp = add(1))), list(...))
    do.call(pmx_ode_subject_eval, args)
  }
  expect_no_error(ok())
  bad(ok(model = "m"))
  bad(ok(subject = "s"))
  bad(ok(theta = sim_theta()[1:2]))
  bad(ok(theta = unname(sim_theta())))
  bad(ok(theta = c(tcl = 0, tv = 3, nope = 1)))
  bad(ok(eta = c(eta.cl = 0)))
  bad(ok(eta = "a"))
  bad(ok(sigma = list(cp = add(1), central = add(1))))
  bad(ok(sigma = list(nope = add(1))))
  bad(ok(tol = c(1e-6)))
  bad(ok(tol = c(-1, 1e-6)))
  bad(ok(tol = c(NA, 1e-6)))
  # A covariate-split: wt in the table is not theta.
  mc <- thread_models()$closed
  sc <- pmx_ode_subjects(mc, population_events(1))
  expect_identical(attr(sc, "theta"), c("tcl", "tv", "tka"))
  expect_identical(attr(sc, "eta"), character())
})
