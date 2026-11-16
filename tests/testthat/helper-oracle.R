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

# Every ok case of the oracle corpus whose model compiles to `class`.
oracle_cases_of_class <- function(class) {
  man <- oracle_read("manifest.csv")
  cases <- unique(man$case[man$status == "ok"])
  keep <- vapply(cases, function(case) {
    m <- tryCatch(pmx_ode_model(oracle_case(case)$rx), pmx_ode_error_compile = function(e) NULL)
    !is.null(m) && identical(m$class, class)
  }, NA)
  cases[keep]
}

# zupmxode's solution of a case, with the fixture's columns.
solve_case <- function(case) {
  c <- oracle_case(case)
  m <- pmx_ode_model(c$rx)
  pmx_ode_solve(m, c$events, c$params)
}

# The 16 hexadecimal digits of IEEE bits back to doubles.
hex_double <- function(h) {
  bytes <- as.raw(strtoi(substring(rep(h, each = 8), seq(1, 15, 2), seq(2, 16, 2)), 16L))
  readBin(bytes, "double", n = length(h), size = 8, endian = "big")
}

# A model with its solver options set, as pmx_ode_method() will (Stage 7).
with_solver <- function(model, ...) {
  model$solve_opts <- utils::modifyList(as.list(model$solve_opts), list(...))
  model
}

# Corpus models where zupmxode's reading differs from an oracle's, each
# with its section in fixtures/oracle/disagreements.md. Their parameters
# make IC50 exactly 0, as R computes it too, so the effect's derivative is
# 0 / 0 at t = 0: zupmxode reports a solver failure where rxode2's
# reordered C arithmetic found a tiny IC50 and went on.
zupmxode_disagreements <- c("lib_indirect_1cpt_inhi_kin_r0rmaxcrmax",
                            "lib_indirect_1cpt_inhi_kout_r0rmaxcrmax")

# Corpus models Dormand-Prince cannot start, besides those: a start its
# step-size control cannot get past (rxode2 warned on both when the
# fixtures were made), or a model that is NaN throughout (both engines
# agree on NA there). LSODA solves the first and agrees on the second.
dopri_stiff_cases <- c("lib_tgi_no_sat_powerLaw", "lib_tgi_sat_genVonBertalanffy")
