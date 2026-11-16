# Prediction tables compared at the design's tolerances (§6.3), the way the
# three-engine tests (§7, §9) compare them. A value is the same when it is
# within `atol` absolutely or within `rtol` of the larger magnitude of the
# two: dazdiff's rule, so the two routes agree (design D16). dazdiff's
# expectation is used when it is installed; otherwise this local one.
#
# `by` names the identity columns; every other numeric column in `expected`
# is compared. Both tables must have the same keys; a key in one only, a
# missing compared column, or a duplicate key fails the expectation.
expect_table_equal <- function(object, expected, by = c("id", "time"),
                               rtol = 1e-6, atol = 1e-8, label = NULL) {
  if (is.null(label)) label <- deparse1(substitute(object))
  # NaN and NA both mean no number: rxode2 writes NA where zupmxode's
  # arithmetic gives NaN, and dazdiff tells the two apart.
  # A column read back from CSV with nothing but NA is logical; it is a
  # column of missing numbers, so it is compared as one.
  nan_to_na <- function(d) {
    for (n in names(d)) {
      if (is.logical(d[[n]]) && all(is.na(d[[n]]))) d[[n]] <- as.double(d[[n]])
      if (is.double(d[[n]])) d[[n]][is.nan(d[[n]])] <- NA
    }
    d
  }
  object <- nan_to_na(object)
  expected <- nan_to_na(expected)
  if (requireNamespace("dazdiff", quietly = TRUE) &&
        !isTRUE(getOption("zupmxode.test.local_tables"))) {
    dazdiff::expect_table_equal(object, expected, by = by,
                                rel_tol = rtol, tolerance = atol)
    return(invisible(object))
  }
  problems <- table_differences(object, expected, by, rtol, atol)
  testthat::expect(
    length(problems) == 0L,
    sprintf("%s differs from the expected table:\n%s", label,
            paste(utils::head(problems, 10L), collapse = "\n"))
  )
  invisible(object)
}

# The differences between two tables, as lines of text; empty when they
# agree. Kept apart from the expectation so its own test can see it.
table_differences <- function(object, expected, by, rtol, atol) {
  missing <- setdiff(names(expected), names(object))
  if (length(missing)) {
    return(sprintf("column `%s` is missing", missing))
  }
  key <- function(d) do.call(paste, c(unname(as.list(d[by])), sep = "\r"))
  ko <- key(object)
  ke <- key(expected)
  dups <- c(ko[duplicated(ko)], ke[duplicated(ke)])
  if (length(dups)) {
    return(sprintf("duplicate key %s", unique(gsub("\r", ", ", dups))))
  }
  only <- c(sprintf("row %s only in object", setdiff(ko, ke)),
            sprintf("row %s only in expected", setdiff(ke, ko)))
  if (length(only)) {
    return(gsub("\r", ", ", only))
  }
  object <- object[match(ke, ko), , drop = FALSE]
  measures <- setdiff(names(expected)[vapply(expected, is.numeric, NA)], by)
  out <- character()
  for (col in measures) {
    a <- object[[col]]
    e <- expected[[col]]
    gap <- abs(a - e)
    same <- (is.na(a) & is.na(e)) |
      (!is.na(gap) & (gap <= atol | gap <= rtol * pmax(abs(a), abs(e)))) |
      (!is.na(a) & !is.na(e) & a == e)
    bad <- which(!same)
    out <- c(out, sprintf("%s at %s: %.17g, expected %.17g", col,
                          gsub("\r", ", ", ke[bad]), a[bad], e[bad]))
  }
  out
}

# A small table for the helper's own tests.
table_example <- function() {
  data.frame(id = c(1L, 1L, 2L), time = c(0, 1, 0), cp = c(10, 5, 0))
}
