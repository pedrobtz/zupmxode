# Event-table errors (design §8): each names the row and the rule, as
# pmx_ode_error_events.

test_that("each rule is refused at its row", {
  cases <- list(
    list(TIME = c(0, 2, 1), EVID = 0, AMT = 0, CMT = 1, row = 3L, rule = "TIME"),
    list(TIME = c(0, NaN), EVID = 0, AMT = 0, CMT = 1, row = 2L, rule = "TIME"),
    list(TIME = 0:1, EVID = c(0, 7), AMT = 0, CMT = 1, row = 2L, rule = "EVID"),
    list(TIME = 0:1, EVID = c(0, 0), AMT = c(0, 5), CMT = 1, row = 2L, rule = "AMT"),
    list(TIME = 0, EVID = 1, AMT = -1, CMT = 1, row = 1L, rule = "AMT"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 3, row = 1L, rule = "CMT"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 1, RATE = -3, row = 1L, rule = "RATE"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 2, RATE = -1, row = 1L, rule = "rate()"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 2, RATE = -2, row = 1L, rule = "dur()"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 1, RATE = 1, DUR = 2, row = 1L, rule = "DUR"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 1, ADDL = 3, row = 1L, rule = "ADDL"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 1, II = -1, row = 1L, rule = "II"),
    list(TIME = 0, EVID = 1, AMT = 1, CMT = 1, SS = 1, row = 1L, rule = "SS"),
    list(TIME = 0, EVID = 0, AMT = 0, CMT = 1, SS = 1, II = 12, row = 1L, rule = "SS")
  )
  for (cs in cases) {
    row <- cs$row
    rule <- cs$rule
    cs$row <- NULL
    cs$rule <- NULL
    err <- do.call(bad, cs)
    expect_identical(err$row, row, label = rule)
    expect_match(err$rule, rule, fixed = TRUE, label = rule)
  }
})

test_that("the model's lag and bioavailability are checked at the dose", {
  d <- data.frame(TIME = 0, EVID = 1, AMT = 1, CMT = "depot")
  expect_identical(bad(TIME = 0, EVID = 1, AMT = 1, CMT = 1,
                       params = sched_params(lag = -1))$row, 1L)
  expect_identical(bad(TIME = 0, EVID = 1, AMT = 1, CMT = 1,
                       params = sched_params(fdepot = NaN))$row, 1L)
})

test_that("an unknown compartment name names the model's compartments", {
  err <- bad(TIME = 0, EVID = 1, AMT = 1, CMT = "gut")
  expect_match(err$rule, "depot, central", fixed = TRUE)
})

test_that("a missing column value and a missing TIME column are refused", {
  err <- bad(TIME = c(0, 1), EVID = 0, AMT = c(0, NA), CMT = 1)
  expect_identical(err$row, 2L)
  expect_error(sched_timeline(sched_model(), data.frame(EVID = 0), sched_params()),
               class = "pmx_ode_error_events")
  expect_error(sched_timeline(sched_model(), data.frame(TIME = 0, AMT = "a"),
                              sched_params()), class = "pmx_ode_error_events")
})

test_that("parameters the model needs and the table does not carry are an error", {
  d <- data.frame(TIME = 0, EVID = 0, AMT = 0, CMT = 1)
  expect_error(sched_timeline(sched_model(), d, c(ka = 1)),
               class = "pmx_ode_error_argument")
  expect_error(sched_timeline(sched_model(), list(), sched_params()),
               class = "pmx_ode_error_argument")
  two <- data.frame(ID = 1:2, TIME = 0, EVID = 0, AMT = 0, CMT = 1)
  expect_error(sched_timeline(sched_model(), two, sched_params()),
               class = "pmx_ode_error_argument")
})

test_that("the C entry refuses columns of the wrong type or length", {
  d <- c(0, 1)
  i <- c(0L, 0L)
  args <- list(d, i, d, i, d, d, d, i, i, d, d, d, d, 2L)
  expect_error(do.call(.Call, c(list(zpx_schedule_r), replace(args, 1, list(1:2)))),
               "bad double column")
  expect_error(do.call(.Call, c(list(zpx_schedule_r), replace(args, 2, list(c(0, 0))))),
               "bad integer column")
  expect_error(do.call(.Call, c(list(zpx_schedule_r), replace(args, 3, list(0)))),
               "bad double column")
})
