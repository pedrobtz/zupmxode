# Models and checks for the scheduler's tests (roadmap Stage 4).

# A two-compartment oral model whose lag, bioavailability, modelled rate
# and duration are all parameters, so each rule can be set from a test.
sched_model <- function() {
  pmx_ode_model("
    alag(depot) <- lag
    f(depot) <- fdepot
    rate(central) <- r0
    dur(central) <- d0
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k * central
  ")
}

sched_params <- function(...) {
  p <- c(lag = 0, fdepot = 1, r0 = 0, d0 = 0, ka = 1, k = 0.1)
  over <- c(...)
  p[names(over)] <- over
  p
}

# The timeline's rows as "time kind cmt" strings, for comparing with what a
# person worked out by hand.
tl <- function(events, ..., model = sched_model()) {
  t <- sched_timeline(model, events, sched_params(...))
  t
}

tl_str <- function(t) paste(format(t$time), t$kind, t$cmt)

# From test-events-rules.R.
ev <- function(...) {
  d <- data.frame(...)
  for (col in c("EVID", "AMT", "RATE", "DUR", "II", "ADDL", "SS")) {
    if (!col %in% names(d)) d[[col]] <- 0
  }
  d
}

# From test-events-errors.R.
bad <- function(..., params = sched_params()) {
  d <- data.frame(...)
  err <- expect_error(sched_timeline(sched_model(), d, params),
                      class = "pmx_ode_error_events")
  invisible(err)
}
