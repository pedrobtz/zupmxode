/* The R glue: the only translation unit that includes R's headers. The
 * engine itself lives in src/core/, which builds without them. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "core/core.h"

static SEXP zupmxode_core_abi_r(void) {
  return Rf_ScalarInteger(zupmxode_core_abi());
}

SEXP zpx_ir_verify_r(SEXP, SEXP);
SEXP zpx_model_new(SEXP, SEXP, SEXP);
SEXP zpx_test_unwind(SEXP, SEXP, SEXP);
SEXP zpx_model_valid(SEXP);
SEXP zpx_model_tier(SEXP);
SEXP zpx_live_models(void);
SEXP zpx_tiers(void);
SEXP zpx_model_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP zpx_ir_eval_r(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP zpx_math_r(SEXP, SEXP, SEXP);
SEXP zpx_solve_all(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP zpx_omp_max(void);
SEXP zpx_schedule_r(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
  {"zupmxode_core_abi_r", (DL_FUNC) &zupmxode_core_abi_r, 0},
  {"zpx_ir_verify_r", (DL_FUNC) &zpx_ir_verify_r, 2},
  {"zpx_model_new", (DL_FUNC) &zpx_model_new, 3},
  {"zpx_test_unwind", (DL_FUNC) &zpx_test_unwind, 3},
  {"zpx_model_valid", (DL_FUNC) &zpx_model_valid, 1},
  {"zpx_model_tier", (DL_FUNC) &zpx_model_tier, 1},
  {"zpx_live_models", (DL_FUNC) &zpx_live_models, 0},
  {"zpx_tiers", (DL_FUNC) &zpx_tiers, 0},
  {"zpx_model_eval", (DL_FUNC) &zpx_model_eval, 6},
  {"zpx_ir_eval_r", (DL_FUNC) &zpx_ir_eval_r, 7},
  {"zpx_math_r", (DL_FUNC) &zpx_math_r, 3},
  {"zpx_solve_all", (DL_FUNC) &zpx_solve_all, 18},
  {"zpx_omp_max", (DL_FUNC) &zpx_omp_max, 0},
  {"zpx_schedule_r", (DL_FUNC) &zpx_schedule_r, 14},
  {NULL, NULL, 0}
};

attribute_visible void R_init_zupmxode(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
