# Population simulation (design §4.5, §5, D7, D19): etas and residual
# error drawn in R with R's RNG, the solve deterministic, so a simulation
# is reproducible from set.seed() at any thread count.

test_that("the same seed gives the same bits at 1 and many threads", {
  m <- sim_model()
  run <- function(threads) {
    withr::local_seed(42)
    with_threads(threads, function() {
      pmx_ode_sim(
        m,
        sim_events(),
        theta = sim_theta(),
        omega = diag(c(0.1, 0.2)),
        sigma = list(cp = prop(0.1)),
        n = 200
      )
    })
  }
  one <- run(1)
  many <- run(test_threads())
  expect_identical(same_bits(many, one), TRUE)
  expect_identical(same_bits(run(1), one), TRUE)
  withr::local_seed(43)
  other <- pmx_ode_sim(
    m,
    sim_events(),
    theta = sim_theta(),
    omega = diag(c(0.1, 0.2)),
    sigma = list(cp = prop(0.1)),
    n = 200
  )
  expect_false(identical(other$cp_obs, one$cp_obs))
})

test_that("each error model is its formula applied to the solve with the same draws", {
  m <- pmx_ode_model({
    cl <- exp(tcl)
    v <- exp(tv)
    ka <- exp(tka)
    cp <- linCmt()
  })
  ev <- sim_events(c(0.5, 1, 2, 4, 8, 12, 24))
  f <- pmx_ode_solve(m, ev, sim_theta())$cp
  models <- list(
    add = list(add(0.7), function(e) f + 0.7 * e),
    prop = list(prop(0.2), function(e) f + f * 0.2 * e),
    combined = list(add(0.7) + prop(0.2), function(e) {
      f + sqrt(0.7^2 + 0.2^2 * f^2) * e
    }),
    lnorm = list(lnorm(0.3), function(e) f * exp(0.3 * e))
  )
  for (name in names(models)) {
    withr::local_seed(7)
    s <- pmx_ode_sim(
      m,
      ev,
      theta = sim_theta(),
      sigma = list(cp = models[[name]][[1]])
    )
    withr::local_seed(7)
    e <- stats::rnorm(length(f))
    expect_identical(s$cp, f, label = name)
    expect_identical(s$cp_obs, models[[name]][[2]](e), label = name)
  }
})

test_that("the combined model is the same whichever way round it is written", {
  expect_identical(add(0.7) + prop(0.2), prop(0.2) + add(0.7))
  expect_error(add(1) + add(2), class = "pmx_ode_error_argument")
  expect_error(lnorm(1) + prop(2), class = "pmx_ode_error_argument")
  expect_error(add(1) + 1, class = "pmx_ode_error_argument")
})

test_that("residual error is drawn for each output in the order sigma lists them", {
  m <- pmx_ode_model({
    cl <- exp(tcl)
    v <- exp(tv)
    ka <- exp(tka)
    cp <- linCmt()
    auc <- 2 * cp
  })
  ev <- sim_events()
  f <- pmx_ode_solve(m, ev, sim_theta())
  withr::local_seed(3)
  s <- pmx_ode_sim(
    m,
    ev,
    theta = sim_theta(),
    sigma = list(auc = add(1), cp = add(2))
  )
  withr::local_seed(3)
  e1 <- stats::rnorm(nrow(f))
  e2 <- stats::rnorm(nrow(f))
  expect_identical(s$auc_obs, f$auc + 1 * e1)
  expect_identical(s$cp_obs, f$cp + 2 * e2)
  # A state takes residual error too.
  withr::local_seed(3)
  s <- pmx_ode_sim(m, ev, theta = sim_theta(), sigma = list(central = add(1)))
  expect_identical(s$central_obs, f$central + e1)
})

test_that("10^4 subjects' etas have omega's covariance, within sampling error", {
  skip_on_slow_checker()
  m <- sim_model()
  omega <- matrix(c(0.09, 0.03, 0.03, 0.16), 2)
  withr::local_seed(11)
  n <- 1e4
  s <- pmx_ode_sim(m, sim_events(24), theta = sim_theta(), omega = omega, n = n)
  eta <- as.matrix(s[c("eta.cl", "eta.v")])
  expect_identical(nrow(eta), as.integer(n))
  # The standard error of a sample covariance is about
  # sqrt((s_ii s_jj + s_ij^2) / n); four of them.
  se <- sqrt((outer(diag(omega), diag(omega)) + omega^2) / n)
  expect_true(all(abs(stats::cov(eta) - omega) < 4 * se))
  expect_true(all(abs(colMeans(eta)) < 4 * sqrt(diag(omega) / n)))
})

test_that("etas are MASS::mvrnorm()'s draws from the same seed", {
  skip_if_not_installed("MASS")
  m <- sim_model()
  omega <- matrix(c(0.09, 0.03, 0.03, 0.16), 2)
  withr::local_seed(5)
  s <- pmx_ode_sim(
    m,
    sim_events(24),
    theta = sim_theta(),
    omega = omega,
    n = 50
  )
  withr::local_seed(5)
  ref <- MASS::mvrnorm(50, c(0, 0), omega)
  expect_identical(unname(as.matrix(s[c("eta.cl", "eta.v")])), unname(ref))
})

test_that("each subject is solved at theta plus its own etas", {
  m <- sim_model()
  ev <- sim_events()
  withr::local_seed(9)
  s <- pmx_ode_sim(m, ev, theta = sim_theta(), omega = diag(c(0.1, 0.2)), n = 5)
  expect_identical(s$id, rep(1:5, each = 4))
  for (i in 1:5) {
    one <- s[s$id == i, ]
    p <- c(sim_theta(), eta.cl = one$eta.cl[1], eta.v = one$eta.v[1])
    expect_identical(
      one$cp,
      pmx_ode_solve(m, ev, p)$cp,
      label = paste("subject", i)
    )
    expect_identical(length(unique(one$eta.cl)), 1L)
  }
})

test_that("omega's names say which parameters are etas, in any order", {
  m <- sim_model()
  ev <- sim_events()
  om <- diag(c(0.1, 0.2))
  withr::local_seed(1)
  a <- pmx_ode_sim(m, ev, theta = sim_theta(), omega = om, n = 3)
  dimnames(om) <- list(c("eta.v", "eta.cl"), c("eta.v", "eta.cl"))
  withr::local_seed(1)
  b <- pmx_ode_sim(m, ev, theta = sim_theta(), omega = om, n = 3)
  # The first column drawn goes to the first name.
  expect_identical(b$eta.v, a$eta.cl)
  expect_identical(b$eta.cl, a$eta.v)
  # A named vector is a diagonal omega.
  withr::local_seed(1)
  d <- pmx_ode_sim(
    m,
    ev,
    theta = sim_theta(),
    omega = c(eta.cl = 0.1, eta.v = 0.2),
    n = 3
  )
  expect_identical(d, a)
  # An eta held at zero variance draws zero.
  s <- pmx_ode_sim(m, ev, theta = sim_theta(), omega = diag(c(0, 0.2)), n = 3)
  expect_true(all(s$eta.cl == 0))
})

test_that("a table of several subjects simulates each once, keeping its records", {
  m <- sim_model()
  ev <- rbind(
    pmx_ode_events(amt = 100, cmt = "depot", time = c(1, 2), id = 7),
    pmx_ode_events(amt = 50, cmt = "depot", time = c(1, 2, 3), id = 3)
  )
  withr::local_seed(2)
  s <- pmx_ode_sim(m, ev, theta = sim_theta(), omega = diag(c(0.1, 0.2)))
  expect_identical(s$id, c(7, 7, 3, 3, 3))
  p <- c(sim_theta(), eta.cl = s$eta.cl[3], eta.v = s$eta.v[3])
  expect_identical(s$cp[3:5], pmx_ode_solve(m, ev[ev$ID == 3, ], p)$cp)
  expect_identical(
    nrow(pmx_ode_sim(
      m,
      ev,
      theta = sim_theta(),
      omega = diag(c(0.1, 0.2)),
      n = 2
    )),
    5L
  )
  expect_error(
    pmx_ode_sim(m, ev, theta = sim_theta(), omega = diag(c(0.1, 0.2)), n = 3),
    class = "pmx_ode_error_argument"
  )
})

test_that("covariates in the table reach every simulated subject", {
  m <- pmx_ode_model({
    cl <- exp(tcl + eta.cl) * (wt / 70)^0.75
    v <- exp(tv)
    ka <- exp(tka)
    cp <- linCmt()
  })
  ev <- sim_events()
  ev$wt <- c(50, 50, 90, 90, 90)
  withr::local_seed(4)
  s <- pmx_ode_sim(m, ev, theta = sim_theta(), omega = c(eta.cl = 0.1), n = 2)
  p <- c(sim_theta(), eta.cl = s$eta.cl[5])
  expect_identical(s$cp[5:8], pmx_ode_solve(m, ev, p)$cp)
})

test_that("no omega and no sigma is the solve at theta, drawing nothing", {
  m <- pmx_ode_model({
    cl <- exp(tcl)
    v <- exp(tv)
    ka <- exp(tka)
    cp <- linCmt()
  })
  ev <- sim_events()
  withr::local_seed(1)
  before <- .Random.seed
  s <- pmx_ode_sim(m, ev, theta = sim_theta(), n = 2)
  expect_identical(.Random.seed, before)
  expect_identical(s$cp, rep(pmx_ode_solve(m, ev, sim_theta())$cp, 2))
})

test_that("a subject the solver fails is NA, and the draws stay where they were", {
  # dy/dt = 50 y^2 from y = 1 has a pole at t = 0.02: subject 1's solver
  # fails there (dopri5; LSODA's silence is #31).
  m <- pmx_ode_method(
    pmx_ode_model("d/dt(central) <- -k * central^2\ncp <- central"),
    "dopri5"
  )
  ev <- rbind(
    pmx_ode_events(amt = 1, cmt = 1, time = c(0.01, 1, 2), id = 1),
    pmx_ode_events(amt = 1, cmt = 1, time = c(0.01, 1, 2), id = 2)
  )
  ev$k <- rep(c(-50, 0.1), each = 4)
  withr::local_seed(8)
  expect_warning(
    s <- pmx_ode_sim(m, ev, sigma = list(cp = add(1))),
    class = "pmx_ode_warning_solver"
  )
  withr::local_seed(8)
  e <- stats::rnorm(6)
  expect_identical(is.na(s$cp_obs), c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(s$cp_obs, s$cp + 1 * e)
  expect_identical(s$cp[4:6], pmx_ode_solve(m, ev[ev$ID == 2, ])$cp)
})

test_that("bad arguments are refused with their class", {
  m <- sim_model()
  ev <- sim_events()
  om <- diag(c(0.1, 0.2))
  sim <- function(...) {
    args <- utils::modifyList(
      list(model = m, events = ev, theta = sim_theta(), omega = om),
      list(...)
    )
    do.call(pmx_ode_sim, args)
  }
  bad <- function(...) expect_error(sim(...), class = "pmx_ode_error_argument")
  bad(model = "a")
  bad(theta = "a")
  bad(theta = c(1, 2))
  bad(omega = matrix(1:6 / 10, 2))
  bad(omega = matrix(c(0.1, 0.05, 0.02, 0.2), 2))
  bad(omega = matrix(c(0.1, 0.5, 0.5, 0.2), 2))
  bad(omega = diag(c(0.1, NA)))
  bad(omega = diag(c(0.1, Inf)))
  bad(omega = "a")
  bad(omega = diag(3) / 10)
  bad(omega = c(eta.cl = 0.1, eta.nope = 0.2))
  bad(omega = c(eta.cl = 0.1, eta.cl = 0.2))
  bad(omega = c(eta.cl = -0.1, eta.v = 0.2))
  bad(theta = c(sim_theta(), eta.cl = 0))
  bad(n = 0)
  bad(n = NA)
  bad(n = 1.5)
  bad(n = Inf)
  bad(n = .Machine$integer.max + 1)
  bad(n = c(1, 2))
  bad(n = "2")
  bad(sigma = prop(0.1))
  bad(sigma = list(prop(0.1)))
  bad(sigma = list(nope = prop(0.1)))
  bad(sigma = list(cp = 0.1))
  bad(sigma = list(cp = prop(0.1), cp = add(1)))
  e2 <- ev
  e2$eta.cl <- 0
  bad(events = e2)
  # Row and column names that disagree.
  om2 <- diag(c(0.1, 0.2))
  dimnames(om2) <- list(c("eta.cl", "eta.v"), c("eta.v", "eta.cl"))
  bad(omega = om2)
  # Unary operators on an error model.
  expect_error(-add(1), class = "pmx_ode_error_argument")
  expect_error(!prop(1), class = "pmx_ode_error_argument")
  expect_error(add(-1), class = "pmx_ode_error_argument")
  expect_error(prop(NA), class = "pmx_ode_error_argument")
  expect_error(lnorm("a"), class = "pmx_ode_error_argument")
  expect_error(add(c(1, 2)), class = "pmx_ode_error_argument")
  expect_error(prop(Inf), class = "pmx_ode_error_argument")
  # A model with a parameter neither theta nor omega gives.
  expect_error(
    sim(theta = c(tcl = 0, tv = 3)),
    class = "pmx_ode_error_argument"
  )
})

test_that("an unnamed omega takes the model's eta parameters in order", {
  m <- pmx_ode_model({
    v <- exp(tv + eta.v)
    cl <- exp(tcl + eta.cl)
    ka <- exp(tka)
    cp <- linCmt()
  })
  withr::local_seed(1)
  s <- pmx_ode_sim(
    m,
    sim_events(),
    theta = sim_theta(),
    omega = diag(c(0, 0.2)),
    n = 2
  )
  expect_true(all(s$eta.v == 0))
  expect_false(all(s$eta.cl == 0))
})

test_that("an empty omega draws nothing, as NULL does", {
  m <- pmx_ode_model({
    cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
    cp <- linCmt()
  })
  withr::local_seed(1)
  a <- pmx_ode_sim(m, sim_events(), theta = sim_theta(), omega = matrix(0, 0, 0), n = 2)
  withr::local_seed(1)
  b <- pmx_ode_sim(m, sim_events(), theta = sim_theta(), n = 2)
  expect_identical(a, b)
})

test_that("an eta named like an event column is refused", {
  m <- pmx_ode_model({
    cl <- exp(tcl + ii); v <- exp(tv); ka <- exp(tka)
    cp <- linCmt()
  })
  expect_error(pmx_ode_sim(m, sim_events(), theta = sim_theta(), omega = c(ii = 0.1)),
               class = "pmx_ode_error_argument")
  # Refused for its name, not later for a value the column swallowed.
  expect_snapshot(
    pmx_ode_sim(m, sim_events(), theta = sim_theta(), omega = c(ii = 0.1)),
    error = TRUE
  )
})
