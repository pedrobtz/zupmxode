#!/usr/bin/env Rscript
# Writes the seed corpus of tests/fuzz/fuzz_ir_verify.c to tests/fuzz/corpus/:
# the hand-built models of tests/testthat/helper-models.R, in the harness's
# input format (a uint32 constant count, the doubles, then the int32 code).
#
#   Rscript tools/make-fuzz-corpus.R      # needs the package loadable
pkgload::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-models.R"))
out <- file.path("tests", "fuzz", "corpus")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# One call of every function of the math table, so that replaying the
# seeds alone reaches the whole table and the reach floor does not depend
# on what a 30-second run happens to generate.
math_calls <- paste(vapply(names(ir_math), function(f) {
  if (ir_math[[f]] >= 32L) paste("call2 r2", f, "r0 r1") else paste("call1 r2", f, "r0")
}, ""), collapse = "\nstoreo 0 r2\n")
seeds <- list(one_cpt = hb_one_cpt()$ir, branches = hb_branches()$ir,
              ops = hb_ops()$ir, empty = ir_asm(list()),
              math = ir_asm(list(rhs = paste0("loadp r0 0\nloadp r1 1\n", math_calls,
                                              "\nstoreo 0 r2")),
                            npar = 2, nout = c(rhs = 1)))
for (n in names(seeds)) {
  ir <- seeds[[n]]
  con <- file(file.path(out, paste0(n, ".bin")), "wb")
  writeBin(length(ir$consts), con, size = 4, endian = "little")
  writeBin(ir$consts, con, size = 8, endian = "little")
  writeBin(ir$code, con, size = 4, endian = "little")
  close(con)
}
cat("wrote", length(seeds), "seeds to", out, "\n")

# Seeds for tests/fuzz/fuzz_scheduler.c: one byte of compartments, then 16
# bytes a record, in the harness's decoding (small(b) is b / 4 below 128).
sched_out <- file.path("tests", "fuzz", "corpus-scheduler")
dir.create(sched_out, recursive = TRUE, showWarnings = FALSE)
record <- function(time = 0, evid = 0, amt = 0, cmt = 1, rate = 0, dur = 0,
                   ii = 0, addl = 0, ss = 0, lag = 0, f = 1, mrate = 0,
                   mdur = 0, chain = 0) {
  q <- function(x) as.integer(round(x * 4))
  c(q(time), evid, amt, cmt + 9L, if (rate == 0) 0L else if (rate < 0) 252L + (rate == -1) else 128L + q(rate),
    if (dur == 0) 0L else 192L + q(dur), if (ii == 0) 0L else 64L + q(ii),
    if (addl == 0) 0L else 128L + addl, ss, if (lag == 0) 0L else 192L + q(lag),
    if (f == 1) 0L else 128L + q(f), q(mrate), q(mdur), 0L, chain, 0L)
}
sched_seeds <- list(
  bolus_obs = list(2, record(0, 1, 100, 1), record(1, 0, cmt = 2, chain = 1)),
  addl_lag = list(2, record(0, 1, 50, 1, ii = 3, addl = 4, lag = 1.5)),
  infusions = list(2, record(0, 1, 100, 2, rate = 10), record(1, 1, 100, 2, rate = -1, mrate = 20, chain = 1),
                   record(1, 1, 100, 2, rate = -2, mdur = 3, chain = 1), record(1, 1, 60, 1, dur = 2, chain = 1)),
  steady = list(2, record(0, 1, 100, 1, ii = 12, ss = 1), record(6, 1, 50, 2, rate = 10, ss = 2, ii = 12),
                record(1, 1, 0, 2, rate = 5, ss = 1, chain = 1)),
  resets = list(3, record(0, 1, 100, 1), record(5, 3, cmt = 0), record(6, 4, 30, 2), record(7, 2, cmt = -3),
                record(8, 2, cmt = 3), record(9, 2, cmt = 0))
)
for (n in names(sched_seeds)) {
  s <- sched_seeds[[n]]
  bytes <- as.raw(c(s[[1]] - 1L, unlist(s[-1])) %% 256L)
  writeBin(bytes, file.path(sched_out, paste0(n, ".bin")))
}
cat("wrote", length(sched_seeds), "scheduler seeds to", sched_out, "\n")
