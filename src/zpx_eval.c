/* Per-subject evaluation (design §4.6, D12): one prepared subject solved at
 * one theta and eta, its predictions and residual variances written into
 * the caller's buffers. zpx_subject_eval() runs on any thread: it calls no
 * R function, touches no R object and takes its memory from the arena.
 * The arena is made on R's main thread, where the interpreter tier's
 * context is lowered. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "zpx_eval.h"

static size_t align8(size_t x) { return (x + 7) & ~(size_t) 7; }

static int method_ok(const zupmxode_subject *s, int method) {
  return method >= 0 && method <= 3 && (s->allowed >> method & 1);
}

const char *zpx_status_message(int status) {
  switch (status) {
  case ZUPMXODE_OK: return "ok";
  case ZUPMXODE_SOLVER_FAILED: return "the solver failed; NaN from where it failed";
  case ZUPMXODE_SS_UNCONVERGED: return "a steady state did not converge";
  case ZUPMXODE_E_VERSION: return "unknown options version";
  case ZUPMXODE_E_ARGUMENT:
    return "a NULL pointer, or an output, error model or method the model and subject do "
           "not have, or a subject prepared for another model";
  case ZUPMXODE_E_BUFFER: return "the buffers are shorter than the subject's observations";
  case ZUPMXODE_E_ARENA: return "the arena is for another model, or too small";
  case ZUPMXODE_E_EVENTS: return "a record breaks an event rule";
  case ZUPMXODE_E_MODEL: return "the closed form's rate constants are not positive";
  default: return "unknown status";
  }
}

int zpx_subject_eval(const zupmxode_model *m, const zupmxode_subject *s, const double *theta,
                     const double *eta, double *pred, double *var,
                     const zupmxode_eval_opts *opts, zupmxode_arena *a, zpx_subject_result *res) {
  if (m == NULL || s == NULL || opts == NULL || a == NULL || pred == NULL || var == NULL)
    return ZUPMXODE_E_ARGUMENT;
  if (opts->version != ZUPMXODE_EVAL_OPTS_VERSION) return ZUPMXODE_E_VERSION;
  int ny = m->dims.ny, npar = m->dims.npar, ncap = m->dims.block[ZPX_POST].nout;
  if (s->fingerprint != m->fingerprint || s->ny != ny || s->npar != npar || s->ncap != ncap)
    return ZUPMXODE_E_ARGUMENT;
  if ((s->ntheta > 0 && theta == NULL) || (s->neta > 0 && eta == NULL))
    return ZUPMXODE_E_ARGUMENT;
  if (opts->output < 0 || opts->output >= ny + ncap) return ZUPMXODE_E_ARGUMENT;
  if (opts->error_model < ZUPMXODE_ERR_ADD || opts->error_model > ZUPMXODE_ERR_LNORM ||
      !(opts->add >= 0) || !(opts->prop >= 0) || !isfinite(opts->add) || !isfinite(opts->prop))
    return ZUPMXODE_E_ARGUMENT;
  int method = opts->method == ZUPMXODE_METHOD_DEFAULT ? s->method : opts->method;
  if (!method_ok(s, method)) return ZUPMXODE_E_ARGUMENT;
  if (opts->n_pred < (size_t) s->nobs) return ZUPMXODE_E_BUFFER;
  if (a->model != m || s->n > a->max_n || s->nobs > a->max_nobs ||
      zpx_solve_arena_size(&a->blocks, s->n, s->entries) > a->solve_size)
    return ZUPMXODE_E_ARENA;

  int n = s->n, nobs = s->nobs;
  memcpy(a->rec, s->rec, (size_t) n * sizeof(zpx_record));
  for (int p = 0; p < npar; p++) {
    double *col = a->par + (size_t) n * p;
    int k = s->index[p];
    switch (s->kind[p]) {
    case 0: for (int i = 0; i < n; i++) col[i] = theta[k]; break;
    case 1: for (int i = 0; i < n; i++) col[i] = eta[k]; break;
    default: memcpy(col, s->cov + (size_t) n * k, (size_t) n * sizeof(double)); break;
    }
  }
  long ncol = 1 + ny + ncap;
  for (long i = 0; i < (long) nobs * ncol; i++) a->out[i] = NAN;
  for (int i = 0; i < nobs; i++) a->orow[i] = -1;
  zpx_solve_opts o = s->opts;
  o.method = method;
  if (opts->rtol > 0) o.rtol = opts->rtol;
  if (opts->atol > 0) o.atol = opts->atol;
  zpx_subject_result r;
  zpx_solve_subject(&a->blocks, &o, a->rec, n, a->par, n, a->out, nobs, a->orow, 0, a->solve,
                    a->solve_size, &r);
  if (res != NULL) *res = r;
  if (r.error == ZPX_SOLVE_E_EVENTS) return ZUPMXODE_E_EVENTS;
  if (r.error != 0) return ZUPMXODE_E_MODEL;

  const double *f = a->out + (size_t) nobs * (size_t) (1 + opts->output);
  double add2 = opts->add * opts->add, prop2 = opts->prop * opts->prop;
  for (int i = 0; i < nobs; i++) {
    pred[i] = f[i];
    switch (opts->error_model) {
    case ZUPMXODE_ERR_ADD:
    case ZUPMXODE_ERR_LNORM: var[i] = isnan(f[i]) ? f[i] : add2; break;
    case ZUPMXODE_ERR_PROP: var[i] = prop2 * (f[i] * f[i]); break;
    default: var[i] = add2 + prop2 * (f[i] * f[i]); break;
    }
  }
  if (r.failed_row >= 0) return ZUPMXODE_SOLVER_FAILED;
  if (r.ss_unconverged > 0) return ZUPMXODE_SS_UNCONVERGED;
  return ZUPMXODE_OK;
}

static void set_blocks(zpx_blocks *b, const zpx_lowered *l, const zpx_model *m) {
  for (int k = 0; k < ZPX_IR_NBLOCKS; k++) b->fn[k] = (zpx_block) l->fn[k];
  b->ny = m->dims.ny;
  b->npar = m->dims.npar;
  b->nwork = m->dims.nwork;
  b->ncap = m->dims.block[ZPX_POST].nout;
}

zupmxode_arena *zpx_arena_new(const zupmxode_model *m, const zupmxode_subject *const *s, int n,
                              int borrow) {
  if (m == NULL || s == NULL) return NULL;
  zupmxode_arena *a = calloc(1, sizeof(*a));
  if (a == NULL) return NULL;
  a->model = m;
  set_blocks(&a->blocks, &m->lowered, m);
  int npar = m->dims.npar, ny = m->dims.ny, ncap = m->dims.block[ZPX_POST].nout;
  int any = 0;
  for (int i = 0; i < n; i++) {
    const zupmxode_subject *x = s[i];
    if (x == NULL || x->fingerprint != m->fingerprint) continue;
    any = 1;
    if (x->n > a->max_n) a->max_n = x->n;
    if (x->nobs > a->max_nobs) a->max_nobs = x->nobs;
    size_t z = zpx_solve_arena_size(&a->blocks, x->n, x->entries);
    if (z > a->solve_size) a->solve_size = z;
  }
  if (!any) {
    free(a);
    return NULL;
  }
  if (m->lowered.tier == ZPX_TIER_INTERP && !borrow) {
    zpx_ir ir;
    int bl, in;
    if (zpx_ir_verify(m->code, m->ncode, m->consts, m->nconsts, &ir, &bl, &in) != ZPX_IR_OK) {
      free(a);
      return NULL;
    }
    void *scratch = malloc(zpx_lower_scratch_size(&ir) + 1);
    const char *msg;
    int where;
    int st = scratch == NULL ? -1 : zpx_lower(&ir, ZPX_TIER_INTERP, &a->own, scratch, &msg, &where);
    free(scratch);
    if (st != ZPX_LOWER_OK) {
      zpx_lower_free(&a->own);
      free(a);
      return NULL;
    }
    a->has_own = 1;
    set_blocks(&a->blocks, &a->own, m);
  }
  size_t solve = align8(a->solve_size + 8);
  size_t rec = align8((size_t) a->max_n * sizeof(zpx_record) + 8);
  size_t par = align8((size_t) a->max_n * (size_t) npar * sizeof(double) + 8);
  size_t out = align8((size_t) a->max_nobs * (size_t) (1 + ny + ncap) * sizeof(double) + 8);
  size_t orow = align8((size_t) a->max_nobs * sizeof(int) + 8);
  char *mem = malloc(solve + rec + par + out + orow);
  if (mem == NULL) {
    zpx_arena_free(a);
    return NULL;
  }
  a->solve = mem;
  a->rec = (zpx_record *) (mem + solve);
  a->par = (double *) (mem + solve + rec);
  a->out = (double *) (mem + solve + rec + par);
  a->orow = (int *) (mem + solve + rec + par + out);
  return a;
}

void zpx_arena_free(zupmxode_arena *a) {
  if (a == NULL) return;
  if (a->has_own) zpx_lower_free(&a->own);
  free(a->solve);
  free(a);
}

void zpx_subject_free(zupmxode_subject *s) {
  if (s == NULL) return;
  free(s->kind);
  free(s->index);
  free(s->rec);
  free(s->cov);
  free(s);
}
