/* The compiled model as the solvers see it (src/zpx_model.c). */
#ifndef ZPX_MODEL_H
#define ZPX_MODEL_H

#define R_NO_REMAP
#include <Rinternals.h>

#include "core/ir.h"
#include "zpx_lower.h"

typedef struct zpx_model {
  zpx_lowered lowered;
  zpx_ir dims; /* counts only: its code and consts pointers are not kept */
} zpx_model;

/* The model behind an external pointer, or NULL. */
zpx_model *zpx_model_get(SEXP x);

#endif
