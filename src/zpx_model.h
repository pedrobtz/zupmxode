/* The compiled model as the solvers see it (src/zpx_model.c). */
#ifndef ZPX_MODEL_H
#define ZPX_MODEL_H

#define R_NO_REMAP
#include <Rinternals.h>

#include "core/ir.h"
#include "zpx_lower.h"

/* The struct behind the C API's opaque zupmxode_model. */
typedef struct zupmxode_model {
  zpx_lowered lowered;
  zpx_ir dims; /* counts only: its code and consts pointers are not kept */
  /* A copy of the IR, owned here, so an arena on the interpreter tier can
   * lower its own context (zumir D10: one context per thread). */
  int *code;
  double *consts;
  size_t ncode, nconsts;
  /* FNV-1a over the IR's bytes: a subject prepared for this model carries
   * it, so one prepared for another model is refused. */
  unsigned long long fingerprint;
} zpx_model;

/* The model behind an external pointer, or NULL. */
zpx_model *zpx_model_get(SEXP x);

#endif
