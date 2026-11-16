# The closed-form solver (design §4.3, D21) against rxode2's linCmt() and
# mrgsolve's PKMODEL (design §7), on every closed-form case of the corpus.

test_that("every closed-form corpus case matches both oracles", {
  skip_under_gctorture()
  cases <- oracle_cases_of_class("closed-form")
  expect_gte(length(cases), 24L)
  for (case in cases) {
    got <- solve_case(case)
    for (engine in c("rxode2", "mrgsolve")) {
      path <- oracle_path(engine, paste0(case, ".csv"))
      if (!file.exists(path)) next
      want <- oracle_read(engine, paste0(case, ".csv"))
      expect_table_equal(got[names(want)], want, label = paste(case, engine))
    }
    # rxode2's linCmt() is a closed form too: the two agree far inside
    # the design's tolerance.
    want <- oracle_read("rxode2", paste0(case, ".csv"))
    expect_table_equal(got[names(want)], want, rtol = 1e-12, atol = 1e-14,
                       label = paste(case, "rxode2, tight"))
  }
})

test_that("1,000 doses stay within 1e-12 of 200-bit superposition", {
  ref <- utils::read.csv(test_path("fixtures", "closed-form-reference.csv"),
                         colClasses = c(ref = "character"))
  ref$ref <- hex_double(ref$ref)
  models <- list(
    pk1_oral = list(text = "ka <- tka\ncl <- tcl\nv <- tv\ncp <- linCmt()",
                    params = c(tka = 0.8, tcl = 2, tv = 20), cmt = "depot"),
    pk2_iv = list(text = "cl <- tcl\nv <- tv\nq <- tq\nvp <- tvp\ncp <- linCmt()",
                  params = c(tcl = 2, tv = 20, tq = 4, tvp = 60), cmt = "central")
  )
  for (name in names(models)) {
    s <- models[[name]]
    r <- ref[ref$model == name, ]
    ev <- pmx_ode_events(amt = 100, ii = 1, addl = 999, cmt = s$cmt, time = r$time)
    got <- pmx_ode_solve(pmx_ode_model(s$text), ev, s$params)
    rel <- abs(got$central - r$ref) / abs(r$ref)
    expect_lt(max(rel), 1e-12, label = name)
  }
})

test_that("SS = 1 is the state repeated dosing converges to", {
  models <- list(
    "ka <- tka\ncl <- tcl\nv <- tv\ncp <- linCmt()",
    "cl <- tcl\nv <- tv\nq <- tq\nvp <- tvp\ncp <- linCmt()",
    "cl <- tcl\nv <- tv\nq <- tq\nvp <- tvp\nq2 <- tq2\nvp2 <- tvp2\ncp <- linCmt()"
  )
  params <- c(tka = 0.8, tcl = 2, tv = 20, tq = 4, tvp = 60, tq2 = 1, tvp2 = 100)
  for (text in models) {
    m <- pmx_ode_model(text)
    cmt <- m$states[1]
    for (rate in c(0, 25)) {
      ss <- pmx_ode_solve(m, pmx_ode_events(amt = 100, ii = 12, ss = 1L, rate = rate,
                                             cmt = cmt, time = c(0, 3, 11.99)),
                          params[m$params])
      # 600 doses reach steady state for the slowest of these systems.
      long <- pmx_ode_solve(m, pmx_ode_events(amt = 100, ii = 12, addl = 599, rate = rate,
                                               cmt = cmt, time = 599 * 12 + c(0, 3, 11.99)),
                            params[m$params])
      expect_equal(ss$cp, long$cp, tolerance = 1e-10, label = paste(text, rate))
    }
  }
  # A constant infusion's steady state: rate / clearance in plasma.
  m <- pmx_ode_model(models[[2]])
  ss <- pmx_ode_solve(m, pmx_ode_events(amt = 0, rate = 10, ss = 1L, cmt = "central",
                                         time = c(0, 50)), params[m$params])
  expect_equal(ss$cp, c(10 / 2, 10 / 2), tolerance = 1e-13)
})

test_that("ka equal to the elimination rate gives the confluent solution", {
  # With ka = k the depot's eigenvalue meets the central one; the solver
  # takes the matrix exponential, and central(t) = D ka t e^(-k t) exactly.
  m <- pmx_ode_model("ka <- tka\nk <- tk\nv <- tv\ncp <- linCmt()")
  t <- c(0.5, 1, 4, 10, 30)
  got <- pmx_ode_solve(m, pmx_ode_events(amt = 100, cmt = "depot", time = t),
                       c(tka = 0.3, tk = 0.3, tv = 10))
  expect_equal(got$central, 100 * 0.3 * t * exp(-0.3 * t), tolerance = 1e-13)
  # Just apart, the eigenvector path gives the same within the gap.
  near <- pmx_ode_solve(m, pmx_ode_events(amt = 100, cmt = "depot", time = t),
                        c(tka = 0.3 * (1 + 1e-5), tk = 0.3, tv = 10))
  expect_equal(near$central, got$central, tolerance = 1e-4)
})

test_that("a covariate changes the rate constants from its record on", {
  m <- pmx_ode_model("cl <- tcl * wt / 70\nv <- tv\ncp <- linCmt()")
  ev <- data.frame(TIME = c(0, 12, 12.01, 24), EVID = c(1, 2, 0, 0),
                   AMT = c(100, 0, 0, 0), CMT = "central", wt = c(70, 140, 140, 140))
  got <- pmx_ode_solve(m, ev, c(tcl = 2, tv = 20))
  a12 <- 100 * exp(-0.1 * 12)
  expect_equal(got$central, c(a12 * exp(-0.2 * 0.01), a12 * exp(-0.2 * 12)),
               tolerance = 1e-13)
})

test_that("several subjects are solved apart", {
  m <- pmx_ode_model("cl <- tcl\nv <- tv\ncp <- linCmt()")
  ev <- rbind(pmx_ode_events(amt = 100, cmt = "central", time = 1, id = 1),
              pmx_ode_events(amt = 50, cmt = "central", time = 1, id = 2))
  got <- pmx_ode_solve(m, ev, c(tcl = 2, tv = 20))
  expect_identical(got$id, c(1, 2))
  expect_equal(got$central, c(100, 50) * exp(-0.1), tolerance = 1e-14)
})

test_that("the solver refuses what it cannot solve, naming why", {
  lin <- pmx_ode_model("cl <- tcl\nv <- tv\ncp <- linCmt()")
  ode <- pmx_ode_model("d/dt(x) <- -k * x")
  # A model with no closed form goes to the general solvers (Stage 6).
  expect_equal(pmx_ode_solve(ode, pmx_ode_events(amt = 1, time = 1), c(k = 1))$x,
               exp(-1), tolerance = 1e-5)
  expect_error(pmx_ode_solve(lin, pmx_ode_events(amt = 1, time = 1), c(tcl = -1, tv = 1)),
               class = "pmx_ode_error_solver")
  # A compartment switched off is the general solvers' (D21): auto takes
  # the subject to LSODA, and the closed form, forced, refuses it.
  off <- data.frame(TIME = 0:1, EVID = c(1, 2), AMT = c(1, 0), CMT = c(1, -1))
  expect_error(pmx_ode_solve(pmx_ode_method(lin, "closed"), off, c(tcl = 1, tv = 1)),
               class = "pmx_ode_error_solver")
  expect_error(pmx_ode_solve(lin, pmx_ode_events(amt = 1, time = 1), c(tcl = 1)),
               class = "pmx_ode_error_argument")
  expect_error(pmx_ode_solve(lin, pmx_ode_events(amt = 1, time = 1), 1),
               class = "pmx_ode_error_argument")
  expect_error(pmx_ode_solve(list(), pmx_ode_events(amt = 1, time = 1)),
               class = "pmx_ode_error_argument")
  bad <- pmx_ode_events(amt = 1, time = 1, addl = 2)
  err <- expect_error(pmx_ode_solve(lin, bad, c(tcl = 1, tv = 1)),
                      class = "pmx_ode_error_events")
  expect_identical(err$row, 1L)
})

test_that("a reset stops an infusion, and its end then changes nothing", {
  m <- pmx_ode_model("cl <- tcl\nv <- tv\ncp <- linCmt()")
  ev <- data.frame(TIME = c(0, 5, 6, 12, 20), EVID = c(1, 3, 0, 0, 0),
                   AMT = c(100, 0, 0, 0, 0), RATE = c(10, 0, 0, 0, 0), CMT = "central")
  got <- pmx_ode_solve(m, ev, c(tcl = 2, tv = 20))
  expect_identical(got$central, c(0, 0, 0))
  # EVID 4: reset and dose; the earlier infusion's end must not take its
  # rate from the new dose's.
  ev <- data.frame(TIME = c(0, 5, 12, 20), EVID = c(1, 4, 0, 0),
                   AMT = c(100, 50, 0, 0), RATE = c(10, 0, 0, 0), CMT = "central")
  got <- pmx_ode_solve(m, ev, c(tcl = 2, tv = 20))
  expect_equal(got$central, 50 * exp(-0.1 * c(7, 15)), tolerance = 1e-13)
})
