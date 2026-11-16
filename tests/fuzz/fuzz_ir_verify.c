/* libFuzzer harness for the IR verifier (src/core/ir.c), built without R by
 * r-actions' fuzz.yml (.github/workflows/fuzz.yaml).
 *
 * Input: a little-endian uint32 count n (taken modulo 4096), then n doubles
 * for the constant pool, then the rest as int32 header and code.
 *
 * Oracle, beyond "no crash": the verifier's promise is that an accepted IR
 * reads and writes only inside the arrays its header declares. So every
 * accepted IR is evaluated, every block, with each array allocated at
 * exactly its declared length (ASan reports any access outside), twice,
 * and the two evaluations must agree bit for bit. A refused IR must name a
 * known error. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "ir.h"

static double *alloc_doubles(int n, double fill) {
  /* At least one byte, so a zero-length array is still a distinct
   * allocation ASan guards. */
  double *p = malloc(n > 0 ? (size_t) n * sizeof(double) : 1);
  if (p == NULL) abort();
  for (int i = 0; i < n; i++) p[i] = fill + i;
  return p;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size < 4) return 0;
  uint32_t n = (uint32_t) data[0] | (uint32_t) data[1] << 8 | (uint32_t) data[2] << 16 |
               (uint32_t) data[3] << 24;
  n %= 4096;
  data += 4;
  size -= 4;
  if ((size_t) n * sizeof(double) > size) n = (uint32_t) (size / sizeof(double));
  double *consts = alloc_doubles((int) n, 0);
  memcpy(consts, data, (size_t) n * sizeof(double));
  data += (size_t) n * sizeof(double);
  size -= (size_t) n * sizeof(double);
  size_t ncode = size / sizeof(int32_t);
  int32_t *code = malloc(ncode > 0 ? ncode * sizeof(int32_t) : 1);
  if (code == NULL) abort();
  memcpy(code, data, ncode * sizeof(int32_t));

  zpx_ir ir;
  int block, insn;
  int e = zpx_ir_verify(code, ncode, consts, n, &ir, &block, &insn);
  if (e != ZPX_IR_OK) {
    if (e < ZPX_IR_E_SHORT || e > ZPX_IR_E_JUMP) abort();
  } else {
    double *y = alloc_doubles(ir.ny, 1.5);
    double *par = alloc_doubles(ir.npar, 0.25);
    for (int b = 0; b < ZPX_IR_NBLOCKS; b++) {
      double *out1 = alloc_doubles(ir.block[b].nout, 0);
      double *out2 = alloc_doubles(ir.block[b].nout, 0);
      double *w1 = alloc_doubles(ir.nwork, 2);
      double *w2 = alloc_doubles(ir.nwork, 2);
      double *regs = alloc_doubles(ir.block[b].nreg, 0);
      zpx_ir_eval(&ir, b, 0.5, y, out1, par, w1, regs);
      zpx_ir_eval(&ir, b, 0.5, y, out2, par, w2, regs);
      if (memcmp(out1, out2, (size_t) ir.block[b].nout * sizeof(double)) != 0 ||
          memcmp(w1, w2, (size_t) ir.nwork * sizeof(double)) != 0)
        abort();
      free(out1);
      free(out2);
      free(w1);
      free(w2);
      free(regs);
    }
    free(y);
    free(par);
  }
  free(code);
  free(consts);
  return 0;
}
