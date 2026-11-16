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
# Every one-argument function at arguments in its tails and poles (the
# harness's parameters are all 0.25, which reach one branch of each): the
# probit's two tails, lgamma's reflection below -15 and its series above
# 15, and the rest, so that a cold cache meets the reach floor too.
tail_args <- c(-20.5, -3, -1e-18, 1e-20, 1e-6, 0.01, 0.99, 1 - 1e-6, 20.5, 1e300)
one_arg <- names(ir_math)[ir_math < 32L]
math_tails <- paste(c(
  unlist(lapply(tail_args, function(x) {
    c(paste("const r0", format(x, digits = 17)),
      vapply(one_arg, function(f) paste0("call1 r2 ", f, " r0\nstoreo 0 r2"), ""))
  }))
), collapse = "\n")
# pow (and the other two-argument functions) at signed zeros, negative
# bases and infinities.
inf_bits <- function(x) if (is.finite(x)) format(x, digits = 17) else ir_bits(x)
two_arg <- names(ir_math)[ir_math >= 32L]
pairs <- expand.grid(x = c(-Inf, -2, 0, 2, Inf), y = c(-Inf, -1, 0.5, 2, 3, Inf))
pow_tails <- paste(unlist(lapply(seq_len(nrow(pairs)), function(i) {
  c(paste("const r0", inf_bits(pairs$x[i])), paste("const r1", inf_bits(pairs$y[i])),
    vapply(two_arg, function(f) paste0("call2 r2 ", f, " r0 r1\nstoreo 0 r2"), ""))
})), collapse = "\n")
# Every opcode once, reading a work slot.
ops_all <- ir_asm(list(rhs = paste(
  "loadt r0", "loady r1 0", "loadp r2 0", "storew 0 r2", "loadw r3 0", "mov r4 r3",
  "add r5 r0 r1", "sub r5 r5 r2", "mul r5 r5 r3", "div r5 r5 r4", "neg r5 r5",
  "lt r6 r0 r1", "le r6 r6 r1", "gt r6 r6 r1", "ge r6 r6 r1", "eq r6 r6 r1",
  "ne r6 r6 r1", "and r6 r6 r1", "or r6 r6 r1", "not r6 r6",
  "jmpz r6 skip", "const r5 2", "skip:", "jmp end", "const r5 3", "end:",
  "storeo 0 r5", sep = "\n")), ny = 1, npar = 1, nwork = 1, nout = c(rhs = 1))
rhs_at <- match("rhs", ir_blocks)
ops_insn0 <- 6L + 3L * length(ir_blocks) + 1L
refusals <- list(
  refuse_count = c(3L, -1L),
  refuse_block = c(6L + 3L * (rhs_at - 1L) + 1L, -1L),
  refuse_opcode = c(ops_insn0, 99L),
  refuse_jump = c(ops_insn0 + 22L * ir_width + 1L, 0L)
)
seeds <- list(one_cpt = hb_one_cpt()$ir, branches = hb_branches()$ir,
              ops = hb_ops()$ir, empty = ir_asm(list()),
              math = ir_asm(list(rhs = paste0("loadp r0 0\nloadp r1 1\n", math_calls,
                                              "\nstoreo 0 r2")),
                            npar = 2, nout = c(rhs = 1)),
              math_tails = ir_asm(list(rhs = math_tails), nout = c(rhs = 1)),
              pow_tails = ir_asm(list(rhs = pow_tails), nout = c(rhs = 1)),
              ops_all = ops_all)
# Refusals, each a copy of ops_all with one word broken: a count in the
# header, a block's register count, an opcode, and a jump backwards.
for (k in names(refusals)) {
  ir <- ops_all
  ir$code[refusals[[k]][1]] <- refusals[[k]][2]
  seeds[[k]] <- ir
}
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
