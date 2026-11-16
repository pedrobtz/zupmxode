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
