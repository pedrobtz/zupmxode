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
seeds <- list(one_cpt = hb_one_cpt()$ir, branches = hb_branches()$ir,
              ops = hb_ops()$ir, empty = ir_asm(list()))
for (n in names(seeds)) {
  ir <- seeds[[n]]
  con <- file(file.path(out, paste0(n, ".bin")), "wb")
  writeBin(length(ir$consts), con, size = 4, endian = "little")
  writeBin(ir$consts, con, size = 8, endian = "little")
  writeBin(ir$code, con, size = 4, endian = "little")
  close(con)
}
cat("wrote", length(seeds), "seeds to", out, "\n")
