#!/usr/bin/env Rscript
# Writes tests/c/replay/<case>.txt: corpus cases replayed by
# tests/c/test_replay.c through src/core/ alone (the scheduler, every
# solver, the reference evaluator for the model's blocks), each compared
# with pmx_ode_solve()'s output: to 1e-12 relative, 1e-8 for LSODA and
# Dormand-Prince, since the platforms' maths libraries differ in the last
# bit and an adaptive solver's steps follow it (tests/c/test_replay.c). The R suite ties that output
# to rxode2 and mrgsolve (tests/testthat/test-closed-form.R, test-ode.R,
# test-expm.R), so a mutant of the C core that changes a number is caught
# here, by the binary the mutation job runs (roadmap Stage 12).
#
#   Rscript tools/make-c-replay.R      # needs the package loadable
#
# Every double is written as the 16 hexadecimal digits of its IEEE bits.
pkgload::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-oracle.R"))
oracle_path <- function(...) file.path("tests", "testthat", "fixtures", "oracle", ...)
out_dir <- file.path("tests", "c", "replay")
unlink(out_dir, recursive = TRUE)
dir.create(out_dir, recursive = TRUE)

bits <- function(x) {
  vapply(as.double(x), function(v) {
    paste(as.character(writeBin(v, raw(), endian = "big")), collapse = "")
  }, "")
}

# Each case and the solver it is replayed with: every solver, steady
# state, infusions, lags, modelled rates, covariates and resets.
cases <- list(
  c("pk1_oral__bolus_addl", "auto"), c("pk2_oral__ss1_bolus", "auto"),
  c("pk3_iv_lin__ss1_infusion", "auto"), c("pk1_oral_lagf__lag_bolus", "auto"),
  c("pk1_iv__infusion_rate", "auto"), c("pk1_iv__ss2_add", "auto"),
  c("pk1_iv_ratedur__rate_modelled", "auto"), c("pk1_iv_wt__covariate_locf", "auto"),
  c("pk1_oral__reset_evid3", "auto"), c("pk1_oral__reset_dose_evid4", "auto"),
  c("pk1_mm__bolus_addl", "lsoda"), c("pk1_mm__ss1_bolus", "lsoda"),
  c("pk1_mm__infusion_rate", "dopri5"), c("pk1_mm__ss1_infusion", "dopri5"),
  c("pk2_iv_lin__ss1_infusion", "expm"), c("pk1_oral__bolus_addl", "lsoda"),
  c("pk1_oral_lin__bolus_addl", "auto"), c("pk2_iv_lin__infusion_rate", "auto"),
  c("pk1_oral_lin__reset_dose_evid4", "auto")
)

for (cs in cases) {
  case <- cs[1]
  c <- oracle_case(case)
  m <- pmx_ode_model(c$rx)
  if (cs[2] != "auto") m <- pmx_ode_method(m, cs[2])
  method <- pmx_ode_method(m)
  rec <- as_records(c$events, m)
  par <- record_params(m, rec, c$params)
  opts <- solve_options(m)
  lc <- if (m$class == "closed-form") lincmt_spec(m) else rep(-1L, 10L)
  out <- pmx_ode_solve(m, c$events, c$params)
  vals <- as.matrix(out[-1])
  ncap <- length(m$captures)
  lines <- c(
    paste("case", paste0(case, "/", method)),
    paste("code", length(m$ir$code), paste(m$ir$code, collapse = " ")),
    paste("consts", length(m$ir$consts), paste(bits(m$ir$consts), collapse = " ")),
    paste("ncap", ncap),
    paste("method", method_codes[[method]]),
    paste("opts", paste(bits(c(opts$rtol, opts$atol, opts$ss_tol)), collapse = " "),
          as.integer(opts$ss_cap)),
    paste("lc", paste(lc, collapse = " ")),
    paste("records", length(rec$time)),
    vapply(seq_along(rec$time), function(i) {
      paste(paste(bits(c(rec$time[i], rec$amt[i], rec$rate[i], rec$dur[i], rec$ii[i])),
                  collapse = " "), rec$evid[i], rec$cmt[i], rec$addl[i], rec$ss[i])
    }, ""),
    paste("par", length(par), paste(bits(par), collapse = " ")),
    paste("out", nrow(vals), ncol(vals), paste(bits(vals), collapse = " "))
  )
  name <- paste0(case, "-", method, ".txt")
  writeLines(lines, file.path(out_dir, name))
  cat(name, ": ", nrow(vals), " observations, ", m$class, "\n", sep = "")
}
