/* Per-subject evaluation (design §4.6): the C API's prepared subjects,
 * arenas and zupmxode_subject_eval(), in src/zpx_eval.c; its R glue in
 * src/zpx_api.c. */
#ifndef ZPX_EVAL_H
#define ZPX_EVAL_H

#include "zupmxode.h"

#include "core/solve.h"
#include "zpx_model.h"

struct zupmxode_subject {
  int n, nobs;  /* records, and observation records among them */
  long entries; /* timeline entries the records make */
  int ny, npar, ncap; /* the dimensions of the model it was prepared for */
  unsigned long long fingerprint; /* and its IR's (zpx_model) */
  int ntheta, neta, ncov;
  int *kind, *index; /* per parameter: 0 theta, 1 eta, 2 covariate; its index there */
  zpx_record *rec;
  double *cov; /* n x ncov, record i's covariate c at cov[i + n * c] */
  int method;  /* the solver chosen when it was prepared */
  int allowed; /* bit k set: solver k may be asked for */
  zpx_solve_opts opts;
};

struct zupmxode_arena {
  const zpx_model *model;
  zpx_blocks blocks;
  zpx_lowered own; /* the interpreter tier's own context, if has_own */
  int has_own;
  int max_n, max_nobs;
  size_t solve_size;
  char *solve;
  zpx_record *rec;
  double *par, *out;
  int *orow;
};

/* zupmxode_subject_eval(), and what went wrong for the R glue: res is
 * the solve's report, or untouched when the call did not get that far. */
int zpx_subject_eval(const zupmxode_model *m, const zupmxode_subject *s, const double *theta,
                     const double *eta, double *pred, double *var,
                     const zupmxode_eval_opts *opts, zupmxode_arena *a, zpx_subject_result *res);

/* An arena for the n subjects; borrow: use the model's own code even on
 * the interpreter tier (one thread, the model's). NULL on failure. Main
 * thread. */
zupmxode_arena *zpx_arena_new(const zupmxode_model *m, const zupmxode_subject *const *s, int n,
                              int borrow);
void zpx_arena_free(zupmxode_arena *a);

void zpx_subject_free(zupmxode_subject *s);
const char *zpx_status_message(int status);
const zupmxode_api *zpx_api_get(int version);

#endif
