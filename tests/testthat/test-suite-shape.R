# The suite runs shuffled, top-level code included (tests/testthat.R), so a
# test file may hold nothing at its top level but test_that() calls: data
# or a function defined there may run after the tests that use it. Shared
# code lives in helper-*.R, which is not shuffled.

test_that("test files hold only test_that() calls at their top level", {
  for (f in list.files(test_path(), pattern = "^test-.*[.]R$", full.names = TRUE)) {
    calls <- as.list(parse(f, keep.source = FALSE))
    other <- Filter(function(e) !(is.call(e) && identical(e[[1]], as.name("test_that"))), calls)
    expect_identical(length(other), 0L,
                     label = paste(basename(f), "top-level statements other than test_that()"))
  }
})
