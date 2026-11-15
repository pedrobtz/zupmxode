# The oracle fixtures of tests/testthat/fixtures/oracle/ (design §7), written
# by tools/make-corpus.R, which sources this file so the generator and the
# tests solve with rxode2 the same way.

oracle_path <- function(...) {
  file.path(testthat::test_path("fixtures", "oracle"), ...)
}

oracle_read <- function(...) {
  utils::read.csv(oracle_path(...), stringsAsFactors = FALSE)
}

# A model's text, parameters and a case's events, as the fixtures hold them.
oracle_case <- function(case) {
  model <- sub("__.*$", "", case)
  p <- oracle_read("models", paste0(model, ".params.csv"))
  list(
    model = model,
    rx = paste(readLines(oracle_path("models", paste0(model, ".R"))),
               collapse = "\n"),
    params = stats::setNames(p$value, p$name),
    events = oracle_read("events", paste0(case, ".csv"))
  )
}

# rxode2 on a NONMEM-style table with compartments by name: id, time and the
# outputs at every observation record, in record order. Covariates are the
# table's extra columns, carried forward from the record that sets them.
oracle_solve_rxode2 <- function(rx_text, events, params, outputs,
                                rtol = 1e-10, atol = 1e-12) {
  mod <- rxode2::rxode2(rx_text)
  d <- events
  names(d) <- tolower(names(d))
  params <- params[setdiff(names(params), names(d))]
  res <- rxode2::rxSolve(mod, events = d, params = params,
                         rtol = rtol, atol = atol, ssRtol = rtol,
                         ssAtol = atol, returnType = "data.frame",
                         addDosing = FALSE, covsInterpolation = "locf")
  res <- as.data.frame(res)
  if (!"id" %in% names(res)) res$id <- 1L
  # rxode2 returns EVID 2 records as rows too; keep the observations.
  kept <- d$evid %in% c(0, 2)
  stopifnot(nrow(res) == sum(kept))
  res <- res[d$evid[kept] == 0, , drop = FALSE]
  out <- data.frame(id = as.integer(res$id), time = res$time,
                    res[outputs], check.names = FALSE)
  rownames(out) <- NULL
  out
}
