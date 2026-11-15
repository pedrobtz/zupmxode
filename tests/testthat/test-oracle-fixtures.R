# The oracle fixtures are the answers later stages are checked against
# (design §7). These tests check the fixtures themselves: that they are the
# whole of what the manifest lists, that every engine disagreement has its
# PREDPP note, and that rxode2 still gives the answers they record. The
# suite runs shuffled, top-level code too, so each test reads what it uses.

test_that("every ok case has its fixture, and every fixture is listed", {
  manifest <- oracle_read("manifest.csv")
  ok <- manifest[manifest$status == "ok", ]
  expected <- sort(file.path(ok$engine, paste0(ok$case, ".csv")))
  found <- sort(c(
    file.path("rxode2", list.files(oracle_path("rxode2"))),
    file.path("mrgsolve", list.files(oracle_path("mrgsolve")))
  ))
  expect_identical(found, expected)
})

test_that("every case solved has its model, parameters and events", {
  manifest <- oracle_read("manifest.csv")
  ok <- unique(manifest$case[manifest$status == "ok"])
  for (case in ok) {
    model <- sub("__.*$", "", case)
    expect_true(file.exists(oracle_path("models", paste0(model, ".R"))),
                label = model)
    expect_true(file.exists(oracle_path("models",
                                        paste0(model, ".params.csv"))),
                label = model)
    expect_true(file.exists(oracle_path("events", paste0(case, ".csv"))),
                label = case)
  }
})

test_that("every engine disagreement has a PREDPP note", {
  cases <- as.character(unique(oracle_read("disagreements.csv")$case))
  notes <- readLines(oracle_path("disagreements.md"))
  headed <- sub("^#+ +", "", grep("^#{2,} ", notes, value = TRUE))
  expect_identical(setdiff(cases, headed), character())
})

# Three cases, one per solver class the corpus holds: a linear ODE with
# additional doses, a closed form at steady state under infusion, and a
# general model with two outputs. Re-solved here with rxode2 by the same
# helper the generator uses; a difference means the fixtures are stale or
# rxode2 changed, and either way they must be regenerated and the change
# read before they are trusted again.
test_that("rxode2 reproduces the fixtures of three corpus cases", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  freshness_cases <- c(
    "pk1_oral__bolus_addl", "pk2_iv_lin__ss1_infusion", "irm1__bolus_single"
  )
  for (case in freshness_cases) {
    c <- oracle_case(case)
    fixture <- oracle_read("rxode2", paste0(case, ".csv"))
    outputs <- setdiff(names(fixture), c("id", "time"))
    res <- suppressMessages(
      oracle_solve_rxode2(c$rx, c$events, c$params, outputs)
    )
    expect_table_equal(res, fixture, rtol = 1e-8, atol = 1e-10)
  }
})
