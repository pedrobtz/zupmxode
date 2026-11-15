# The table helper is what every three-engine test rests on: a helper that
# passed a wrong table would turn each of those failures into a pass. Each
# test builds its own table: the suite runs shuffled, top-level code too.

test_that("a cell further than rtol from expected fails", {
  expected <- table_example()
  object <- expected
  object$cp[2] <- 5 * (1 + 1.1e-6)
  expect_length(table_differences(object, expected, c("id", "time"),
                                  1e-6, 1e-8), 1L)
})

test_that("a cell within rtol of expected passes", {
  expected <- table_example()
  object <- expected
  object$cp[2] <- 5 * (1 + 0.9e-6)
  expect_length(table_differences(object, expected, c("id", "time"),
                                  1e-6, 1e-8), 0L)
})

test_that("near zero, atol decides", {
  expected <- table_example()
  object <- expected
  object$cp[3] <- 0.9e-8
  expect_length(table_differences(object, expected, c("id", "time"),
                                  1e-6, 1e-8), 0L)
  object$cp[3] <- 1.1e-8
  expect_length(table_differences(object, expected, c("id", "time"),
                                  1e-6, 1e-8), 1L)
})

test_that("NA against a number fails, NA against NA passes", {
  expected <- table_example()
  object <- expected
  object$cp[1] <- NA
  expect_length(table_differences(object, expected, c("id", "time"),
                                  1e-6, 1e-8), 1L)
  both <- expected
  both$cp[1] <- NA
  expect_length(table_differences(both, both, c("id", "time"),
                                  1e-6, 1e-8), 0L)
})

test_that("a missing row, column or a duplicate key fails", {
  expected <- table_example()
  expect_match(table_differences(expected[-2, ], expected, c("id", "time"),
                                 1e-6, 1e-8), "only in expected")
  expect_match(table_differences(expected["id"], expected, c("id", "time"),
                                 1e-6, 1e-8), "missing", all = FALSE)
  expect_match(table_differences(expected[c(1, 1, 2, 3), ], expected,
                                 c("id", "time"), 1e-6, 1e-8),
               "duplicate key")
})

test_that("row order does not matter", {
  expected <- table_example()
  expect_length(table_differences(expected[3:1, ], expected,
                                  c("id", "time"), 1e-6, 1e-8), 0L)
})

test_that("the expectation fails on a perturbed table on either route", {
  expected <- table_example()
  object <- expected
  object$cp[2] <- 6
  withr::local_options(zupmxode.test.local_tables = TRUE)
  expect_failure(expect_table_equal(object, expected))
  expect_success(expect_table_equal(expected, expected))
  skip_if_not_installed("dazdiff")
  withr::local_options(zupmxode.test.local_tables = FALSE)
  expect_failure(expect_table_equal(object, expected))
  expect_success(expect_table_equal(expected, expected))
})
