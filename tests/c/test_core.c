/* Builds src/core/ with no R headers on the include path and checks it.
 * Exit status 0 on success. */
#include <math.h>
#include <stdio.h>

#include "core.h"
#include "ir.h"

static int failed = 0;

static void check(int ok, const char *what) {
  if (!ok) {
    fprintf(stderr, "FAIL: %s\n", what);
    failed = 1;
  }
}

/* rhs computes exp(par[0]) * y[0] into register 3 (no store, so the test
 * reads the register file): ny 1, npar 1, one constant, four registers. */
static int32_t ir_code[] = {
  ZPX_IR_MAGIC, ZPX_IR_VERSION, 1, 1, 0, 1,
  0, 0, 0,  0, 0, 0,  4, 1, 4,  0, 0, 0,  0, 0, 0,
  ZPX_OP_LOADP, 0, 0, 0, 0,
  ZPX_OP_CALL1, 1, ZPX_FN_EXP, 0, 0,
  ZPX_OP_LOADY, 2, 0, 0, 0,
  ZPX_OP_MUL, 3, 1, 2, 0,
};
static const double ir_consts[] = {0.5};

int main(void) {
  zpx_ir ir;
  int block, insn;
  size_t n = sizeof(ir_code) / sizeof(ir_code[0]);

  check(zupmxode_core_abi() == ZUPMXODE_CORE_ABI, "core ABI");

  /* Accepted and evaluated, then altered one operand at a time. */
  check(zpx_ir_verify(ir_code, n, ir_consts, 1, &ir, &block, &insn) == ZPX_IR_OK,
        "a well-formed IR is accepted");
  double y = 3, par = 0, out = -1, regs[4];
  zpx_ir_eval(&ir, ZPX_RHS, 0, &y, &out, &par, NULL, regs);
  check(regs[3] == 3.0, "the reference evaluator computes exp(0) * 3");

  ir_code[n - 5 + 1] = 4; /* mul's destination: register 4 of 4 */
  check(zpx_ir_verify(ir_code, n, ir_consts, 1, &ir, &block, &insn) == ZPX_IR_E_OPERAND,
        "a register past nreg is refused");
  check(block == ZPX_RHS && insn == 3, "the refusal names block and instruction");
  ir_code[n - 5 + 1] = 3;

  check(zpx_ir_verify(ir_code, n - 1, ir_consts, 1, &ir, &block, &insn) == ZPX_IR_E_LENGTH,
        "a truncated IR is refused");
  check(isnan(zpx_gamma(0.0)) && isnan(zpx_gamma(-2.0)), "gamma is NaN at its poles");
  check(zpx_pow(0.0, -1.0) == INFINITY && zpx_pow(-0.0, -1.0) == INFINITY,
        "a zero base has no sign in pow");

  if (!failed) printf("core: ok (ABI %d, IR %d)\n", ZUPMXODE_CORE_ABI, ZPX_IR_VERSION);
  return failed;
}
