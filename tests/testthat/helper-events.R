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

# A population of n subjects, the same doses, each with its own weight,
# for the parallel runtime's tests.
population_events <- function(n, times = c(0.5, 2, 6, 12, 24, 36)) {
  ev <- do.call(rbind, lapply(seq_len(n), function(i) {
    e <- pmx_ode_events(amt = 100, ii = 12, addl = 2, cmt = 1L, time = times, id = i)
    e$wt <- 50 + i
    e
  }))
  ev
}

# Runs f with the thread count set to n, then restores it.
with_threads <- function(n, f) {
  old <- pmx_ode_threads(n)
  on.exit(pmx_ode_threads(old))
  f()
}

# Allocations R records while expr runs (test-alloc.R).
profmem_count <- function(expr) {
  path <- withr::local_tempfile()
  utils::Rprofmem(path, threshold = 0)
  force(expr)
  utils::Rprofmem(NULL)
  lines <- readLines(path, warn = FALSE)
  length(lines)
}

# A model of each class, every one using a covariate (test-threads.R).
thread_models <- function() {
  list(
    closed = pmx_ode_model("cl <- tcl * (wt / 70)^0.75\nv <- tv\nka <- tka\ncp <- linCmt()"),
    linear = pmx_ode_model("d/dt(depot) <- -ka * depot\nd/dt(central) <- ka * depot - k * wt / 70 * central"),
    general = pmx_ode_model("d/dt(depot) <- -ka * depot\nd/dt(central) <- ka * depot - vmax * central / (km + central) * wt / 70")
  )
}

# Threads for the tests: 8, or 2 where CRAN's two-core limit applies.
test_threads <- function() {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) 8L else 2L
}

# TRUE when two tables hold the same doubles, bit for bit (NA and NaN
# included), else where they first differ.
same_bits <- function(a, b) {
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) return("shape")
  for (k in names(a)) {
    if (!identical(writeBin(as.double(a[[k]]), raw()), writeBin(as.double(b[[k]]), raw()))) {
      return(k)
    }
  }
  TRUE
}

# A table's subjects copied k times, each copy with its own IDs.
replicate_subjects <- function(events, k) {
  id <- if ("ID" %in% names(events)) "ID" else if ("id" %in% names(events)) "id" else NULL
  if (is.null(id)) {
    events$ID <- 1L
    id <- "ID"
  }
  ids <- unique(events[[id]])
  do.call(rbind, lapply(seq_len(k), function(i) {
    e <- events
    e[[id]] <- match(e[[id]], ids) + (i - 1L) * length(ids)
    e
  }))
}

# expr's value and the messages of the warnings it raised, muffled.
with_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(expr, warning = function(x) {
    w[length(w) + 1L] <<- conditionMessage(x)
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

# The thread tests that solve thousands of subjects or the whole corpus
# many times over: past the budget under gctorture (10 minutes) and
# valgrind (an hour, with its threads run one at a time). Under both,
# the 40-subject test and the subject-order test still cross the region.
skip_on_slow_checker <- function() {
  testthat::skip_if(Sys.getenv("R_ACTIONS_WORKFLOW") %in% c("gctorture", "valgrind"),
                    "thousands of subjects; the smaller thread tests run here")
}
