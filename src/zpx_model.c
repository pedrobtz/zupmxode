/* The compiled model as R sees it (design §4.5, §4.7): an external pointer
 * owning the lowered code, freed by its finalizer. Every entry returns a
 * status list rather than raising, where the failure is the IR's or MIR's,
 * so R/model.R raises the classed condition (design §8). */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <stdlib.h>
#include <string.h>

#include "zpx_model.h"

/* Live models, for the tests of the finalizer. Main thread only. */
static int live_models = 0;

static SEXP model_tag(void) {
  static SEXP tag = NULL;
  if (tag == NULL) tag = Rf_install("zupmxode_model");
  return tag;
}

static void model_free(zpx_model *m) {
  if (m == NULL) return;
  zpx_lower_free(&m->lowered);
  free(m->code);
  free(m->consts);
  free(m);
  live_models--;
}

static void model_finalize(SEXP x) {
  zpx_model *m = R_ExternalPtrAddr(x);
  R_ClearExternalPtr(x);
  model_free(m);
}

static SEXP status_list(const char *kind, const char *message, int block, int insn) {
  const char *names[] = {"status", "message", "block", "insn", ""};
  SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 0, Rf_mkString(kind));
  SET_VECTOR_ELT(res, 1, Rf_mkString(message != NULL ? message : ""));
  SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(block));
  SET_VECTOR_ELT(res, 3, Rf_ScalarInteger(insn));
  UNPROTECT(1);
  return res;
}

static int verify_sexp(SEXP code, SEXP consts, zpx_ir *ir, int *block, int *insn) {
  if (TYPEOF(code) != INTSXP || TYPEOF(consts) != REALSXP) return ZPX_IR_E_SHORT;
  return zpx_ir_verify(INTEGER(code), (size_t) XLENGTH(code), REAL(consts),
                       (size_t) XLENGTH(consts), ir, block, insn);
}

SEXP zpx_ir_verify_r(SEXP code, SEXP consts) {
  zpx_ir ir;
  int block, insn;
  int e = verify_sexp(code, consts, &ir, &block, &insn);
  return status_list(e == ZPX_IR_OK ? "ok" : "ir", zpx_ir_error_message(e), block, insn);
}

/* tier: 1 native, 2 interp. The external pointer and its finalizer exist
 * before the context does, so an R error from here on cannot strand it.
 * `fail_after` is a test hook: nonzero raises an R error once the model is
 * lowered and owned, which the unwind test catches. */
static SEXP model_new(SEXP code, SEXP consts, SEXP s_tier, int fail_after) {
  zpx_ir ir;
  int block, insn, tier = Rf_asInteger(s_tier);
  int e = verify_sexp(code, consts, &ir, &block, &insn);
  if (e != ZPX_IR_OK) return status_list("ir", zpx_ir_error_message(e), block, insn);
  const char *why = zpx_lower_api();
  if (why != NULL) return status_list("api", why, -1, -1);
  if (!zpx_lower_tier_available(tier)) return status_list("no_tier", "", -1, -1);
  SEXP x = PROTECT(R_MakeExternalPtr(NULL, model_tag(), R_NilValue));
  R_RegisterCFinalizerEx(x, model_finalize, TRUE);
  zpx_model *m = calloc(1, sizeof(*m));
  if (m == NULL) Rf_error("cannot allocate a model");
  live_models++;
  R_SetExternalPtrAddr(x, m);
  m->dims = ir;
  m->ncode = (size_t) XLENGTH(code);
  m->nconsts = (size_t) XLENGTH(consts);
  m->code = malloc(m->ncode * sizeof(int) + 1);
  m->consts = malloc(m->nconsts * sizeof(double) + 1);
  if (m->code == NULL || m->consts == NULL) Rf_error("cannot allocate a model");
  memcpy(m->code, INTEGER(code), m->ncode * sizeof(int));
  memcpy(m->consts, REAL(consts), m->nconsts * sizeof(double));
  unsigned long long h = 1469598103934665603ULL;
  const unsigned char *bytes = (const unsigned char *) m->code;
  for (size_t i = 0; i < m->ncode * sizeof(int); i++) h = (h ^ bytes[i]) * 1099511628211ULL;
  bytes = (const unsigned char *) m->consts;
  for (size_t i = 0; i < m->nconsts * sizeof(double); i++) h = (h ^ bytes[i]) * 1099511628211ULL;
  m->fingerprint = h;
  void *scratch = R_alloc(zpx_lower_scratch_size(&ir) + 1, 1);
  const char *message;
  int st = zpx_lower(&ir, tier, &m->lowered, scratch, &message, &block);
  if (st != ZPX_LOWER_OK) {
    SEXP res = PROTECT(status_list(st == ZPX_LOWER_E_NO_TIER ? "no_tier" : "lowering",
                                   message, block, -1));
    R_ClearExternalPtr(x);
    model_free(m);
    UNPROTECT(2);
    return res;
  }
  if (fail_after) Rf_error("zupmxode test hook: an error after lowering");
  UNPROTECT(1);
  return x;
}

SEXP zpx_model_new(SEXP code, SEXP consts, SEXP tier) {
  return model_new(code, consts, tier, 0);
}

SEXP zpx_test_unwind(SEXP code, SEXP consts, SEXP tier) {
  return model_new(code, consts, tier, 1);
}

zpx_model *zpx_model_get(SEXP x) {
  if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != model_tag()) return NULL;
  return R_ExternalPtrAddr(x);
}

static zpx_model *model_get(SEXP x) { return zpx_model_get(x); }

SEXP zpx_model_valid(SEXP x) { return Rf_ScalarLogical(model_get(x) != NULL); }

SEXP zpx_model_tier(SEXP x) {
  zpx_model *m = model_get(x);
  return Rf_ScalarInteger(m == NULL ? NA_INTEGER : m->lowered.tier);
}

SEXP zpx_live_models(void) { return Rf_ScalarInteger(live_models); }

SEXP zpx_tiers(void) {
  const char *names[] = {"native", "interp", ""};
  SEXP res = PROTECT(Rf_mkNamed(LGLSXP, names));
  int ok = zpx_lower_api() == NULL;
  LOGICAL(res)[0] = ok && zpx_lower_tier_available(ZPX_TIER_NATIVE);
  LOGICAL(res)[1] = ok && zpx_lower_tier_available(ZPX_TIER_INTERP);
  UNPROTECT(1);
  return res;
}

/* Checks the lengths of the arrays a block reads and writes, as the
 * lowered code and the reference evaluator trust them. */
static void check_args(const zpx_ir *ir, int b, SEXP y, SEXP par, SEXP work) {
  if (b < 0 || b >= ZPX_IR_NBLOCKS) Rf_error("no block %d", b);
  if (TYPEOF(y) != REALSXP || XLENGTH(y) != ir->ny) Rf_error("`y` must be %d doubles", ir->ny);
  if (TYPEOF(par) != REALSXP || XLENGTH(par) != ir->npar)
    Rf_error("`par` must be %d doubles", ir->npar);
  if (TYPEOF(work) != REALSXP || XLENGTH(work) != ir->nwork)
    Rf_error("`work` must be %d doubles", ir->nwork);
}

static SEXP eval_result(SEXP out, SEXP work) {
  const char *names[] = {"out", "work", ""};
  SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 0, out);
  SET_VECTOR_ELT(res, 1, work);
  UNPROTECT(1);
  return res;
}

/* Runs block b of a lowered model once: list(out, work), work updated on
 * a copy. */
SEXP zpx_model_eval(SEXP x, SEXP s_block, SEXP t, SEXP y, SEXP par, SEXP work) {
  zpx_model *m = model_get(x);
  if (m == NULL) Rf_error("the model has no code: it was saved and not restored");
  int b = Rf_asInteger(s_block);
  check_args(&m->dims, b, y, par, work);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, m->dims.block[b].nout));
  SEXP w = PROTECT(Rf_duplicate(work));
  memset(REAL(out), 0, sizeof(double) * (size_t) XLENGTH(out));
  m->lowered.fn[b](Rf_asReal(t), REAL(y), REAL(out), REAL(par), REAL(w));
  SEXP res = eval_result(out, w);
  UNPROTECT(2);
  return res;
}

/* The same through the reference evaluator, on an IR that is verified
 * here first. */
SEXP zpx_ir_eval_r(SEXP code, SEXP consts, SEXP s_block, SEXP t, SEXP y, SEXP par, SEXP work) {
  zpx_ir ir;
  int block, insn;
  int e = verify_sexp(code, consts, &ir, &block, &insn);
  if (e != ZPX_IR_OK) Rf_error("%s", zpx_ir_error_message(e));
  int b = Rf_asInteger(s_block);
  check_args(&ir, b, y, par, work);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, ir.block[b].nout));
  SEXP w = PROTECT(Rf_duplicate(work));
  double *regs = (double *) R_alloc((size_t) ir.block[b].nreg + 1, sizeof(double));
  memset(REAL(out), 0, sizeof(double) * (size_t) XLENGTH(out));
  zpx_ir_eval(&ir, b, Rf_asReal(t), REAL(y), REAL(out), REAL(par), REAL(w), regs);
  SEXP res = eval_result(out, w);
  UNPROTECT(2);
  return res;
}

/* The math table's function f on each element of x (and y): for the test
 * that compares it with R's own functions. */
SEXP zpx_math_r(SEXP s_f, SEXP x, SEXP y) {
  int f = Rf_asInteger(s_f);
  R_xlen_t n = XLENGTH(x);
  SEXP res = PROTECT(Rf_allocVector(REALSXP, n));
  zpx_fn1 f1 = zpx_math1(f);
  zpx_fn2 f2 = zpx_math2(f);
  if (f1 == NULL && f2 == NULL) Rf_error("no math function %d", f);
  if (f2 != NULL && XLENGTH(y) != n) Rf_error("`y` must have the length of `x`");
  for (R_xlen_t i = 0; i < n; i++)
    REAL(res)[i] = f1 != NULL ? f1(REAL(x)[i]) : f2(REAL(x)[i], REAL(y)[i]);
  UNPROTECT(1);
  return res;
}
