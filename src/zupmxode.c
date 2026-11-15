/* The R glue: the only translation unit that includes R's headers. The
 * engine itself lives in src/core/, which builds without them. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "core/core.h"

static SEXP zupmxode_core_abi_r(void) {
  return Rf_ScalarInteger(zupmxode_core_abi());
}

static const R_CallMethodDef call_methods[] = {
  {"zupmxode_core_abi_r", (DL_FUNC) &zupmxode_core_abi_r, 0},
  {NULL, NULL, 0}
};

void R_init_zupmxode(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
