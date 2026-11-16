# The scheduler against timelines worked out by hand, one per rule of
# PREDPP's event semantics that it implements (design §6.2).

test_that("a bolus lands at its time, scaled by bioavailability", {
  t <- tl(ev(TIME = c(0, 2), EVID = c(1, 0), AMT = c(100, 0), CMT = c("depot", "depot")),
          fdepot = 0.8)
  expect_identical(t$kind, c("bolus", "obs"))
  expect_identical(t$amount[1], 80)
  expect_identical(t$row, c(1L, 2L))
})

test_that("a lag shifts the dose, and only the dose", {
  t <- tl(ev(TIME = c(0, 1), EVID = c(1, 0), AMT = c(100, 0), CMT = "depot"), lag = 1.5)
  expect_identical(tl_str(t), c("1.0 obs 1", "1.5 bolus 1"))
})

test_that("ADDL adds doses at k * II from the dose's time", {
  for (addl in c(0L, 1L, 1000L)) {
    t <- tl(ev(TIME = 0, EVID = 1, AMT = 10, CMT = "depot", II = 0.1, ADDL = addl))
    expect_identical(nrow(t), addl + 1L)
    # Computed as time + k * II, not summed: the last of 1000 doses at
    # 0.1 is 100 exactly as R computes 1000 * 0.1, with no drift.
    expect_identical(t$time, (0:addl) * 0.1)
  }
})

test_that("RATE gives an infusion whose duration bioavailability changes", {
  t <- tl(ev(TIME = 0, EVID = 1, AMT = 100, CMT = "central", RATE = 25), fdepot = 0.5)
  expect_identical(t$kind, c("inf_start", "inf_end"))
  expect_identical(t$time, c(0, 4))
  expect_identical(t$rate, c(25, 25))
})

test_that("bioavailability scales an infusion into a depot by its duration", {
  t <- tl(ev(TIME = 0, EVID = 1, AMT = 100, CMT = "depot", RATE = 25), fdepot = 0.5)
  expect_identical(t$time, c(0, 2))
})

test_that("RATE -1 takes the model's rate, RATE -2 its duration", {
  t <- tl(ev(TIME = 0, EVID = 1, AMT = 100, CMT = "central", RATE = -1), r0 = 50)
  expect_identical(t$time, c(0, 2))
  expect_identical(t$rate[1], 50)
  t <- tl(ev(TIME = 0, EVID = 1, AMT = 100, CMT = "central", RATE = -2), d0 = 4)
  expect_identical(t$time, c(0, 4))
  expect_identical(t$rate[1], 25)
})

test_that("RATE -2 keeps the model's duration; bioavailability changes the rate", {
  m <- pmx_ode_model("f(depot) <- fd\ndur(depot) <- d0\nrate(depot) <- r0\nd/dt(depot) <- -depot")
  d <- data.frame(TIME = 0, EVID = 1, AMT = 100, CMT = "depot", RATE = -2)
  t <- sched_timeline(m, d, c(fd = 0.5, d0 = 4, r0 = 1))
  expect_identical(t$time, c(0, 4))
  expect_identical(t$rate[1], 12.5)
  # And RATE -1 keeps the model's rate: F changes the duration.
  d$RATE <- -1
  t <- sched_timeline(m, d, c(fd = 0.5, d0 = 4, r0 = 10))
  expect_identical(t$time, c(0, 5))
})

test_that("DUR keeps its duration; bioavailability changes the rate", {
  t <- tl(ev(TIME = 0, EVID = 1, AMT = 100, CMT = "depot", DUR = 4), fdepot = 0.5)
  expect_identical(t$time, c(0, 4))
  expect_identical(t$rate[1], 12.5)
})

test_that("a steady-state dose is a steady-state entry, then the dose", {
  for (ss in 1:2) {
    t <- tl(ev(TIME = 0, EVID = 1, AMT = 100, CMT = "depot", II = 12, SS = ss))
    expect_identical(t$kind, c("ss", "bolus"))
    expect_identical(t$ss[1], ss)
    expect_identical(c(t$amount[1], t$ii[1]), c(100, 12))
  }
  # A steady-state infusion with II 0 runs on: no end.
  t <- tl(ev(TIME = 0, EVID = 1, AMT = 0, CMT = "central", RATE = 10, SS = 1))
  expect_identical(t$kind, c("ss", "inf_start"))
  # With II, the infusion of each interval ends.
  t <- tl(ev(TIME = 0, EVID = 1, AMT = 100, CMT = "central", RATE = 20, II = 24, SS = 1))
  expect_identical(tl_str(t), c("0 ss 2", "0 inf_start 2", "5 inf_end 2"))
})

test_that("EVID 2 marks a time, switching a compartment off or on", {
  t <- tl(ev(TIME = c(0, 1, 2, 3), EVID = c(1, 2, 2, 2), AMT = c(10, 0, 0, 0),
             CMT = c(1, -2, 2, 0)))
  expect_identical(t$kind, c("bolus", "off", "on", "mark"))
  expect_identical(t$cmt, c(1L, 2L, 2L, 0L))
})

test_that("EVID 3 resets; EVID 4 resets and doses at once", {
  t <- tl(ev(TIME = c(0, 24, 30), EVID = c(1, 3, 4), AMT = c(100, 0, 50), CMT = "depot"))
  expect_identical(tl_str(t), c(" 0 bolus 1", "24 reset 0", "30 reset 0", "30 bolus 1"))
})

test_that("at one time: infusion ends, then generated doses, then records in order", {
  # An observation written before a dose at its time is read before it.
  t <- tl(ev(TIME = c(5, 5), EVID = c(0, 1), AMT = c(0, 10), CMT = "depot"))
  expect_identical(t$kind, c("obs", "bolus"))
  # A dose its lag moves onto a record's time comes before the record.
  t <- tl(ev(TIME = c(0, 2), EVID = c(1, 0), AMT = c(10, 0), CMT = "depot"), lag = 2)
  expect_identical(t$kind, c("bolus", "obs"))
  # An infusion ending when the next starts: the end first.
  t <- tl(ev(TIME = c(0, 2), EVID = c(1, 1), AMT = c(100, 100), CMT = "central",
             RATE = c(50, 50)))
  expect_identical(tl_str(t), c("0 inf_start 2", "2 inf_end 2", "2 inf_start 2", "4 inf_end 2"))
})

test_that("a covariate sets the lag at the record that carries it", {
  m <- pmx_ode_model("alag(depot) <- wt / 100\nd/dt(depot) <- -ka * depot")
  d <- ev(TIME = c(0, 10), EVID = 1, AMT = 10, CMT = "depot", wt = c(50, 150))
  t <- sched_timeline(m, d, c(ka = 1))
  expect_identical(t$time, c(0.5, 11.5))
})

test_that("timeline entries: 1e6 are allowed, 1e6 + 1 refused", {
  d <- function(addl) ev(TIME = 0, EVID = 1, AMT = 1, CMT = "depot", II = 1e-3, ADDL = addl)
  expect_identical(nrow(tl(d(1e6 - 2))), as.integer(1e6 - 1))
  expect_identical(nrow(tl(d(1e6 - 1))), as.integer(1e6))
  expect_error(tl(d(1e6)), class = "pmx_ode_error_events")
  # The limit counts every entry, infusion ends included.
  inf <- ev(TIME = 0, EVID = 1, AMT = 1, CMT = "central", RATE = 1, II = 2, ADDL = 5e5)
  expect_error(tl(inf), class = "pmx_ode_error_events")
})

test_that("pmx_ode_events() builds the table a person would write", {
  e <- pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
                      time = c(0, 12, 24))
  expect_s3_class(e, "pmx_ode_events")
  expect_identical(e$EVID, c(1L, 0L, 0L, 0L))
  expect_identical(e$TIME, c(0, 0, 12, 24))
  expect_identical(e$ADDL, c(6L, 0L, 0L, 0L))
  t <- tl(e)
  expect_identical(sum(t$kind == "bolus"), 7L)
  e2 <- pmx_ode_events(amt = c(100, 50), dose_time = c(0, 12), cmt = "central",
                       rate = c(0, 10), time = 6)
  expect_identical(e2$TIME, c(0, 6, 12))
  expect_identical(e2$RATE, c(0, 0, 10))
})
