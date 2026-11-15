# The verifier (src/core/ir.c) is what makes IR safe to lower: each error it
# reports, reached, and each limit straddled at limit - 1, limit and
# limit + 1 (assurance Limits rows).

test_that("a well-formed IR is accepted", {
  expect_identical(verdict(raw_ir())$status, "ok")
  expect_identical(verdict(raw_ir(insn(1, 0, 0, 0, 0), nconst = 1L), 1)$status,
                   "ok")
})

test_that("each kind of malformed IR is refused with its error", {
  msg <- function(code, consts = double()) verdict(code, consts)$message
  expect_match(msg(integer()), "shorter")
  expect_match(msg(raw_ir()[-21]), "shorter")
  x <- raw_ir(); x[1] <- 1L
  expect_match(msg(x), "magic")
  x <- raw_ir(); x[2] <- 0L
  expect_match(msg(x), "version")
  expect_match(msg(raw_ir(ny = -1L)), "count")
  expect_match(msg(c(raw_ir(), 0L)), "length")
  expect_match(msg(raw_ir(nconst = 2L), 1), "constant pool")
  expect_match(msg(raw_ir(insn(99, 0, 0, 0, 0))), "unknown opcode")
  # An operand past its range, one per kind of operand.
  expect_match(msg(raw_ir(insn(8, 1, 0, 0, 0))), "operand")   # register
  expect_match(msg(raw_ir(insn(1, 0, 0, 0, 0))), "operand")   # constant
  expect_match(msg(raw_ir(insn(3, 0, 0, 0, 0))), "operand")   # state
  expect_match(msg(raw_ir(insn(6, 1, 0, 0, 0))), "operand")   # output
  expect_match(msg(raw_ir(insn(20, 0, 15, 0, 0))), "operand") # function
  expect_match(msg(raw_ir(insn(21, 0, 0, 0, 0))), "operand") # 1-arg as 2
  expect_match(msg(raw_ir(insn(8, 0, 0, 0, 1))), "operand")   # unused != 0
  # Jumps: to itself, backward, past the end; forward to the end is fine.
  expect_match(msg(raw_ir(insn(40, 0, 0, 0, 0))), "jump")
  expect_match(msg(raw_ir(insn(8, 0, 0, 0, 0, 40, 0, 0, 0, 0))), "jump")
  expect_match(msg(raw_ir(insn(40, 2, 0, 0, 0))), "jump")
  expect_identical(verdict(raw_ir(insn(40, 1, 0, 0, 0)))$status, "ok")
  expect_match(msg(raw_ir(insn(41, 0, 0, 0, 0))), "jump")
})

test_that("the verdict names the block and instruction", {
  v <- verdict(raw_ir(insn(8, 0, 0, 0, 0, 99, 0, 0, 0, 0)))
  expect_identical(c(v$block, v$insn), c(2L, 1L))
})

test_that("instructions per block: 1e6 is accepted, 1e6 + 1 refused", {
  mov <- c(8L, 0L, 0L, 0L, 0L)
  for (n in c(1e6 - 1, 1e6, 1e6 + 1)) {
    code <- raw_ir()
    code[15] <- as.integer(n)
    code <- c(code, rep(mov, n))
    expect_identical(verdict(code)$status == "ok", n <= 1e6, label = n)
  }
})

test_that("registers per block: 65536 is accepted, 65537 refused", {
  for (n in c(65535L, 65536L, 65537L)) {
    code <- raw_ir(insn(8, n - 1L, 0, 0, 0), nreg = n)
    expect_identical(verdict(code)$status == "ok", n <= 65536L, label = n)
  }
})

test_that("slots: 65536 states, parameters, work and outputs, not 65537", {
  for (n in c(65535L, 65536L, 65537L)) {
    ok <- n <= 65536L
    expect_identical(verdict(raw_ir(ny = n))$status == "ok", ok)
    expect_identical(verdict(raw_ir(npar = n))$status == "ok", ok)
    expect_identical(verdict(raw_ir(nwork = n))$status == "ok", ok)
    expect_identical(verdict(raw_ir(nout = n))$status == "ok", ok)
  }
})

test_that("constants: 1e6 is accepted, 1e6 + 1 refused", {
  for (n in c(1e6 - 1, 1e6, 1e6 + 1)) {
    code <- raw_ir(nconst = as.integer(n))
    expect_identical(verdict(code, double(n))$status == "ok", n <= 1e6,
                     label = n)
  }
})

test_that("what is not an integer and a double vector is refused", {
  expect_identical(.Call(zpx_ir_verify_r, as.double(raw_ir()), double())$status,
                   "ir")
  expect_identical(.Call(zpx_ir_verify_r, raw_ir(), integer())$status, "ir")
})

test_that("the assembler's output is what the verifier accepts", {
  for (hb in list(hb_one_cpt(), hb_branches(), hb_ops())) {
    expect_identical(ir_verify(hb$ir)$status, "ok")
  }
  expect_error(ir_asm(list(rhs = "frob r0")), "unknown IR opcode")
  expect_error(ir_asm(list(rhs = "mov r0")), "takes 2 operands")
  expect_error(ir_asm(list(rhs = "jmp nowhere")), "bad operand")
})
