/* Lowering a verified IR to code through zumir (src/zpx_lower.c). */
#ifndef ZPX_LOWER_H
#define ZPX_LOWER_H

#include <stddef.h>

#include "core/ir.h"

enum { ZPX_TIER_NATIVE = 1, ZPX_TIER_INTERP = 2 };

enum {
  ZPX_LOWER_OK = 0,
  ZPX_LOWER_E_API,     /* zumir's table is not loaded */
  ZPX_LOWER_E_NO_TIER, /* the tier asked for does not run here */
  ZPX_LOWER_E_MIR      /* MIR reported an error; *message has it */
};

typedef void (*zpx_block_fn)(double t, const double *y, double *out, const double *par,
                             double *work);

struct zumir_ctx;

/* A lowered model: the context owning the code, and one function per
 * block, each a plain C function pointer on either tier. */
typedef struct zpx_lowered {
  struct zumir_ctx *z;
  zpx_block_fn fn[ZPX_IR_NBLOCKS];
  int tier;
} zpx_lowered;

/* Loads zumir's table; NULL, or why it cannot be loaded. Main thread. */
const char *zpx_lower_api(void);
/* Whether zumir reports that `tier` compiles and runs here. */
int zpx_lower_tier_available(int tier);
/* Bytes of scratch zpx_lower() needs for ir. */
size_t zpx_lower_scratch_size(const zpx_ir *ir);
/* Lowers ir on `tier` into *out. On failure *out may hold a context that
 * zpx_lower_free() must free; *message is MIR's message or NULL, and
 * *where_block the block that failed or -1. */
int zpx_lower(const zpx_ir *ir, int tier, zpx_lowered *out, void *scratch, const char **message,
              int *where_block);
/* Frees the code; NULL or twice does nothing. */
void zpx_lower_free(zpx_lowered *l);

#endif
