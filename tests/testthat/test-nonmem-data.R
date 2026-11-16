# pmx_ode_solve() on a NONMEM-style data frame (design §5, §6.2): the same
# records as pmx_ode_events() makes give the same output, and NONMEM's own
# conventions (DV, MDV, lowercase names, numbered compartments) are read
# as NONMEM reads them.

test_that("a corpus data set gives what the same doses from pmx_ode_events() give", {
  for (case in c(
    "pk1_iv__bolus_addl",
    "pk1_iv__infusion_rate",
    "pk2_oral__bolus_addl"
  )) {
    c <- oracle_case(case)
    m <- pmx_ode_model(c$rx)
    d <- c$events
    dose <- d[d$EVID == 1, ]
    obs <- d[d$EVID == 0, ]
    ev <- pmx_ode_events(
      amt = dose$AMT,
      cmt = dose$CMT,
      time = obs$TIME,
      dose_time = dose$TIME,
      ii = dose$II,
      addl = dose$ADDL,
      rate = dose$RATE,
      ss = dose$SS
    )
    expect_identical(
      pmx_ode_solve(m, d, c$params),
      pmx_ode_solve(m, ev, c$params),
      label = case
    )
  }
})

test_that("NONMEM's DV, MDV, lowercase names and numbered compartments read as NONMEM reads them", {
  c <- oracle_case("pk1_iv__bolus_addl")
  m <- pmx_ode_model(c$rx)
  want <- pmx_ode_solve(m, c$events, c$params)
  d <- c$events
  nm <- data.frame(
    id = d$ID,
    time = d$TIME,
    amt = d$AMT,
    cmt = match(d$CMT, m$states),
    rate = d$RATE,
    ii = d$II,
    addl = d$ADDL,
    ss = d$SS,
    dv = ifelse(d$EVID == 0, 1.5, 0),
    mdv = as.integer(d$EVID != 0),
    site = "A",
    stringsAsFactors = FALSE
  )
  expect_identical(pmx_ode_solve(m, nm, c$params), want)
  # An MDV record with no dose marks a time and is not an observation.
  extra <- nm[1, ]
  extra$time <- 0.1
  extra$amt <- 0
  before <- nm$time <= extra$time
  nm2 <- rbind(nm[before, ], extra, nm[!before, ])
  expect_identical(pmx_ode_solve(m, nm2, c$params), want)
  # Given EVID, MDV is not read.
  nm3 <- nm
  nm3$evid <- d$EVID
  nm3$mdv <- 1L
  expect_identical(pmx_ode_solve(m, nm3, c$params), want)
})

test_that("a NONMEM table of several subjects keeps each subject's records", {
  c <- oracle_case("pk1_iv__bolus_addl")
  m <- pmx_ode_model(c$rx)
  one <- pmx_ode_solve(m, c$events, c$params)
  d <- replicate_subjects(c$events, 3)
  d$ID <- c(30, 10, 20)[d$ID]
  got <- pmx_ode_solve(m, d, c$params)
  expect_identical(got$id, rep(c(30, 10, 20), each = nrow(one)))
  expect_identical(got$cp, rep(one$cp, 3))
})

test_that("with no EVID, a steady-state infusion with no amount is a dose, as NONMEM reads it", {
  m <- pmx_ode_model("d/dt(central) <- -cl / v * central\ncp <- central / v")
  p <- c(cl = 1, v = 10)
  with_evid <- data.frame(ID = 1, TIME = 0:2, AMT = 0, RATE = c(10, 0, 0),
                          SS = c(1, 0, 0), II = 0, CMT = 1, EVID = c(1, 0, 0))
  nm <- with_evid
  nm$EVID <- NULL
  nm$MDV <- c(1, 0, 0)
  expect_identical(pmx_ode_solve(m, nm, p), pmx_ode_solve(m, with_evid, p))
})
