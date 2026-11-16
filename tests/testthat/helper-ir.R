# Shared by test-ir-verify.R and test-ir-serialise.R (the suite runs
# shuffled, so nothing a test uses may live at a test file's top level).

verdict <- function(code, consts = double()) {
  .Call(zpx_ir_verify_r, as.integer(code), as.double(consts))
}

# An IR with one block `rhs` of the given instructions (a matrix of 5-int
# rows) and header counts.
raw_ir <- function(insns = matrix(integer(), 0, 5), ny = 0L, npar = 0L,
                   nwork = 0L, nconst = 0L, nreg = 1L, nout = 1L) {
  desc <- c(0L, 0L, 0L, 0L, 0L, 0L, nreg, nout, nrow(insns), 0L, 0L, 0L,
            0L, 0L, 0L)
  c(ir_magic, ir_version, ny, npar, nwork, nconst, desc, as.vector(t(insns)))
}

insn <- function(...) matrix(as.integer(c(...)), ncol = 5, byrow = TRUE)


save_restore <- function(x) {
  path <- withr::local_tempfile(fileext = ".rds", .local_envir = parent.frame())
  saveRDS(x, path)
  readRDS(path)
}

