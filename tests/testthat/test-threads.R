# The parallel runtime (design §4.4, D8): subjects split across threads
# give, bit for bit, what one thread gives, on every solver and tier.

test_that("many threads give the bits 1 thread gives, every solver, every tier", {
  skip_under_gctorture()
  ev <- population_events(40)
  p <- c(tcl = 2, tv = 20, tka = 0.8, ka = 0.8, k = 0.1, vmax = 10, km = 2)
  n <- test_threads()
  for (name in names(thread_models())) {
    base <- thread_models()[[name]]
    methods <- if (name == "general") c("lsoda", "dopri5") else pmx_ode_method(base)
    for (method in methods) {
      for (tier in available_tiers()) {
        m <- pmx_ode_method(pmx_ode_tier(base, tier), method)
        one <- with_threads(1, function() pmx_ode_solve(m, ev, p[m$params]))
        many <- with_threads(n, function() pmx_ode_solve(m, ev, p[m$params]))
        expect_identical(same_bits(many, one), TRUE, label = paste(name, method, tier))
      }
    }
  }
})

test_that("every corpus case, copied over subjects, is the same at 1 and many threads", {
  skip_on_slow_checker()
  man <- oracle_read("manifest.csv")
  n <- test_threads()
  solved <- 0L
  for (case in unique(man$case[man$status == "ok"])) {
    if (sub("__.*$", "", case) %in% zupmxode_disagreements) next
    c <- oracle_case(case)
    # The corpus holds models outside the subset too (test-oracle-*.R);
    # only those that compile are solved.
    m <- tryCatch(pmx_ode_model(c$rx), pmx_ode_error_compile = function(e) NULL)
    if (is.null(m)) next
    ev <- replicate_subjects(c$events, 2L * n)
    for (tier in available_tiers()) {
      mt <- pmx_ode_tier(m, tier)
      # A solver failure (lib_tgi_sat_genVonBertalanffy fails LSODA at the
      # default tolerances) must be the same failure at any thread count.
      one <- with_warnings(with_threads(1, function() pmx_ode_solve(mt, ev, c$params)))
      many <- with_warnings(with_threads(n, function() pmx_ode_solve(mt, ev, c$params)))
      expect_identical(same_bits(many$value, one$value), TRUE, label = paste(case, tier))
      expect_identical(many$warnings, one$warnings, label = paste(case, tier))
      solved <- solved + 1L
    }
  }
  # Most of the corpus compiles: a sweep that solved nothing proves nothing.
  expect_gt(solved, 100L)
})

test_that("the subjects come back in order, each its own", {
  m <- thread_models()$closed
  ev <- population_events(30, times = 12)
  got <- with_threads(3, function() pmx_ode_solve(m, ev, c(tcl = 2, tv = 20, tka = 0.8)))
  expect_identical(got$id, 1:30)
  alone <- vapply(1:30, function(i) {
    pmx_ode_solve(m, ev[ev$ID == i, ], c(tcl = 2, tv = 20, tka = 0.8))$cp
  }, 0)
  expect_identical(got$cp, alone)
})

test_that("an interrupt between chunks returns control, and the next solve works", {
  skip_on_slow_checker()
  m <- pmx_ode_method(thread_models()$general, "lsoda")
  ev <- population_events(2000, times = seq(0.5, 36, by = 0.5))
  p <- c(ka = 0.8, vmax = 10, km = 2)
  # setTimeLimit() is checked where an interrupt is: in
  # R_CheckUserInterrupt(), between chunks of subjects. With the limit at
  # a quarter of the whole solve, the solve stops near that quarter; were
  # the C side not checking, R would raise the same error only once the
  # .Call had run to the end.
  took <- system.time(with_threads(2, function() pmx_ode_solve(m, ev, p)))[["elapsed"]]
  skip_if(took < 0.2, "the solve is too fast here to interrupt reliably")
  started <- proc.time()[["elapsed"]]
  stopped <- tryCatch({
    setTimeLimit(elapsed = took / 4, transient = TRUE)
    with_threads(2, function() pmx_ode_solve(m, ev, p))
    FALSE
  }, error = function(e) grepl("time limit", conditionMessage(e)))
  setTimeLimit(elapsed = Inf)
  waited <- proc.time()[["elapsed"]] - started
  expect_true(stopped)
  expect_lt(waited, 0.7 * took)
  gc()
  small <- with_threads(2, function() pmx_ode_solve(m, population_events(3), p))
  expect_identical(nrow(small), 18L)
})

test_that("pmx_ode_threads() takes a whole number of threads", {
  old <- pmx_ode_threads()
  expect_identical(old, 1L)
  expect_identical(pmx_ode_threads(3), old)
  expect_identical(pmx_ode_threads(), 3L)
  pmx_ode_threads(old)
  for (bad in list(0, -1, 1.5, NA, "2", c(1, 2), Inf)) {
    expect_error(pmx_ode_threads(bad), class = "pmx_ode_error_argument")
  }
  expect_identical(pmx_ode_threads(), old)
})

test_that("Linux builds run in parallel: the tests above are not one thread", {
  # R on Linux has OpenMP in every build CI uses; elsewhere it may not
  # (the CRAN macOS build), and there everything runs on one thread.
  skip_if_not(identical(Sys.info()[["sysname"]], "Linux"), "OpenMP is not certain here")
  expect_gt(.Call(zpx_omp_max), 0L)
})
