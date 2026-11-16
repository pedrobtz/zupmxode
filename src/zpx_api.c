/* The C API's R glue (design §4.6): prepared subjects as external
 * pointers, pmx_ode_subject_eval()'s .Call, the function table registered
 * for other packages, and the tests' callers. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "zpx_eval.h"

static SEXP subject_tag(void) {
  static SEXP tag = NULL;
  if (tag == NULL) tag = Rf_install("zupmxode_subject");
  return tag;
}

static void subject_finalize(SEXP x) {
  zupmxode_subject *s = R_ExternalPtrAddr(x);
  R_ClearExternalPtr(x);
  zpx_subject_free(s);
}

static zupmxode_subject *subject_get(SEXP x) {
  if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != subject_tag()) return NULL;
  return R_ExternalPtrAddr(x);
}

static const zupmxode_model *model_from_r(void *x) { return zpx_model_get((SEXP) x); }
static const zupmxode_subject *subject_from_r(void *x) { return subject_get((SEXP) x); }

static int api_eval(const zupmxode_model *m, const zupmxode_subject *s, const double *theta,
                    const double *eta, double *pred, double *var, const zupmxode_eval_opts *opts,
                    zupmxode_arena *a) {
  return zpx_subject_eval(m, s, theta, eta, pred, var, opts, a, NULL);
}

static zupmxode_arena *api_arena_new(const zupmxode_model *m, const zupmxode_subject *const *s,
                                     int n) {
  return zpx_arena_new(m, s, n, 0);
}

static int api_nobs(const zupmxode_subject *s) { return s == NULL ? -1 : s->nobs; }

static const zupmxode_api api = {ZUPMXODE_API_VERSION, api_eval,        api_arena_new,
                                 zpx_arena_free,       api_nobs,        model_from_r,
                                 subject_from_r,       zpx_status_message};

/* The table for a caller built against header version `version`: NULL
 * when that header is newer than this build. */
const zupmxode_api *zpx_api_get(int version) {
  return version >= 1 && version <= ZUPMXODE_API_VERSION ? &api : NULL;
}

/* A subject from its records (the columns), covariates (n x ncov),
 * parameter map (kind, index per parameter), the theta and eta counts, the
 * solver and the solvers allowed (a bit mask), the solve options c(rtol,
 * atol, ss_tol, ss_cap) and the closed form's spec. The external pointer
 * and its finalizer exist before the memory does. */
SEXP zpx_subject_new(SEXP model, SEXP time, SEXP evid, SEXP amt, SEXP cmt, SEXP rate,
                     SEXP dur, SEXP ii, SEXP addl, SEXP ss, SEXP cov, SEXP kind, SEXP index,
                     SEXP counts, SEXP s_method, SEXP s_opts, SEXP s_lc) {
  zpx_model *m = zpx_model_get(model);
  if (m == NULL) Rf_error("the model has no code: it was saved and not restored");
  int n = (int) XLENGTH(time), npar = m->dims.npar;
  if (XLENGTH(kind) != npar || XLENGTH(index) != npar || TYPEOF(kind) != INTSXP ||
      TYPEOF(index) != INTSXP)
    Rf_error("bad parameter map");
  if (TYPEOF(cov) != REALSXP || Rf_nrows(cov) != n) Rf_error("bad covariates");
  if (TYPEOF(counts) != INTSXP || XLENGTH(counts) != 2) Rf_error("bad counts");
  if (TYPEOF(s_method) != INTSXP || XLENGTH(s_method) != 2) Rf_error("bad method");
  if (TYPEOF(s_opts) != REALSXP || XLENGTH(s_opts) != 4) Rf_error("bad options");
  if (TYPEOF(s_lc) != INTSXP || XLENGTH(s_lc) != 10) Rf_error("bad linCmt spec");
  int ncov = Rf_ncols(cov);
  SEXP x = PROTECT(R_MakeExternalPtr(NULL, subject_tag(), R_NilValue));
  R_RegisterCFinalizerEx(x, subject_finalize, TRUE);
  zupmxode_subject *s = calloc(1, sizeof(*s));
  if (s == NULL) Rf_error("cannot allocate a subject");
  R_SetExternalPtrAddr(x, s);
  s->kind = malloc((size_t) npar * sizeof(int) + 1);
  s->index = malloc((size_t) npar * sizeof(int) + 1);
  s->rec = malloc((size_t) n * sizeof(zpx_record) + 1);
  s->cov = malloc((size_t) n * (size_t) ncov * sizeof(double) + 1);
  if (s->kind == NULL || s->index == NULL || s->rec == NULL || s->cov == NULL)
    Rf_error("cannot allocate a subject");
  s->n = n;
  s->ny = m->dims.ny;
  s->npar = npar;
  s->ncap = m->dims.block[ZPX_POST].nout;
  s->fingerprint = m->fingerprint;
  s->ntheta = INTEGER(counts)[0];
  s->neta = INTEGER(counts)[1];
  s->ncov = ncov;
  memcpy(s->kind, INTEGER(kind), (size_t) npar * sizeof(int));
  memcpy(s->index, INTEGER(index), (size_t) npar * sizeof(int));
  memcpy(s->cov, REAL(cov), (size_t) n * (size_t) ncov * sizeof(double));
  for (int i = 0; i < n; i++) {
    zpx_record *r = &s->rec[i];
    r->time = REAL(time)[i];
    r->amt = REAL(amt)[i];
    r->rate = REAL(rate)[i];
    r->dur = REAL(dur)[i];
    r->ii = REAL(ii)[i];
    r->evid = INTEGER(evid)[i];
    r->cmt = INTEGER(cmt)[i];
    r->addl = INTEGER(addl)[i];
    r->ss = INTEGER(ss)[i];
    /* Placeholders, as for the count in src/zpx_solve.c: pre fills them
     * on each call. */
    r->lag = 0;
    r->f = 1;
    r->mrate = 1;
    r->mdur = 1;
    s->nobs += r->evid == 0;
  }
  int error, where;
  long c = zpx_sched_count(s->rec, n, s->ny, &error, &where);
  s->entries = c < 0 ? 0 : c;
  s->method = INTEGER(s_method)[0];
  s->allowed = INTEGER(s_method)[1];
  memset(&s->opts, 0, sizeof(s->opts));
  s->opts.rtol = REAL(s_opts)[0];
  s->opts.atol = REAL(s_opts)[1];
  s->opts.ss_tol = REAL(s_opts)[2];
  s->opts.ss_cap = (long) REAL(s_opts)[3];
  memcpy(s->opts.lc, INTEGER(s_lc), sizeof(s->opts.lc));
  UNPROTECT(1);
  return x;
}

SEXP zpx_subject_valid(SEXP x) { return Rf_ScalarLogical(subject_get(x) != NULL); }

static void arena_finalize(SEXP x) {
  zupmxode_arena *a = R_ExternalPtrAddr(x);
  R_ClearExternalPtr(x);
  zpx_arena_free(a);
}

/* An arena owned by an external pointer, so an R error cannot strand it. */
static zupmxode_arena *held_arena(SEXP held, const zpx_model *m, const zupmxode_subject *const *s,
                                  int n, int borrow) {
  R_RegisterCFinalizerEx(held, arena_finalize, TRUE);
  zupmxode_arena *a = zpx_arena_new(m, s, n, borrow);
  if (a == NULL) Rf_error("cannot make an arena for the subject");
  R_SetExternalPtrAddr(held, a);
  return a;
}

static void opts_from(SEXP s_opts, zupmxode_eval_opts *o, size_t n_pred) {
  if (TYPEOF(s_opts) != REALSXP || XLENGTH(s_opts) != 7) Rf_error("bad options");
  const double *v = REAL(s_opts);
  memset(o, 0, sizeof(*o));
  o->version = ZUPMXODE_EVAL_OPTS_VERSION;
  o->output = (int) v[0];
  o->error_model = (int) v[1];
  o->add = v[2];
  o->prop = v[3];
  o->method = (int) v[4];
  o->rtol = v[5];
  o->atol = v[6];
  o->n_pred = n_pred;
}

static const double *maybe_real(SEXP x) { return XLENGTH(x) > 0 ? REAL(x) : NULL; }

/* pmx_ode_subject_eval(): list(status, pred, var, error_row, message,
 * failed_at, failed_row). opts = c(output, error model, add, prop, method,
 * rtol, atol). */
SEXP zpx_subject_eval_r(SEXP model, SEXP subject, SEXP theta, SEXP eta, SEXP s_opts) {
  zpx_model *m = zpx_model_get(model);
  zupmxode_subject *s = subject_get(subject);
  if (m == NULL || s == NULL) Rf_error("the model or subject was saved and not restored");
  if (TYPEOF(theta) != REALSXP || TYPEOF(eta) != REALSXP) Rf_error("bad theta or eta");
  if (XLENGTH(theta) != s->ntheta || XLENGTH(eta) != s->neta) Rf_error("bad theta or eta");
  zupmxode_eval_opts o;
  opts_from(s_opts, &o, (size_t) s->nobs);
  const char *names[] = {"status", "pred", "var", "error_row", "message", "failed_at",
                         "failed_row", ""};
  SEXP ans = PROTECT(Rf_mkNamed(VECSXP, names));
  SEXP pred = Rf_allocVector(REALSXP, s->nobs);
  SET_VECTOR_ELT(ans, 1, pred);
  SEXP var = Rf_allocVector(REALSXP, s->nobs);
  SET_VECTOR_ELT(ans, 2, var);
  if (s->fingerprint != m->fingerprint) {
    /* Another model's subject: no arena can be made for it, and the call
     * itself says so. */
    SET_VECTOR_ELT(ans, 0, Rf_ScalarInteger(ZUPMXODE_E_ARGUMENT));
    UNPROTECT(1);
    return ans;
  }
  SEXP held = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  const zupmxode_subject *one[1] = {s};
  zupmxode_arena *a = held_arena(held, m, one, 1, 1);
  zpx_subject_result r;
  memset(&r, 0, sizeof(r));
  r.failed_at = NAN;
  r.failed_row = -1;
  int st = zpx_subject_eval(m, s, maybe_real(theta), maybe_real(eta), REAL(pred), REAL(var), &o,
                            a, &r);
  R_ClearExternalPtr(held);
  zpx_arena_free(a);
  SET_VECTOR_ELT(ans, 0, Rf_ScalarInteger(st));
  SET_VECTOR_ELT(ans, 3, Rf_ScalarInteger(r.error_row));
  SET_VECTOR_ELT(ans, 4, Rf_mkString(r.message != NULL ? r.message : ""));
  SET_VECTOR_ELT(ans, 5, Rf_ScalarReal(r.failed_at));
  SET_VECTOR_ELT(ans, 6, Rf_ScalarInteger(r.failed_row));
  UNPROTECT(2);
  return ans;
}

/* ---- The tests' callers ------------------------------------------------ */

/* Every subject of `subjects` (a list of external pointers) evaluated
 * through the registered table, from `threads` threads, one arena each:
 * list(pred, status), pred one vector per subject. */
SEXP zpx_test_eval_threads(SEXP model, SEXP subjects, SEXP theta, SEXP eta, SEXP s_opts,
                           SEXP s_threads) {
  const zupmxode_api *(*get)(int) =
      (const zupmxode_api *(*)(int)) R_GetCCallable("zupmxode", "zupmxode_api");
  const zupmxode_api *A = get(ZUPMXODE_API_VERSION);
  if (A == NULL) Rf_error("no C API table");
  const zupmxode_model *m = A->model_from_r(model);
  if (m == NULL) Rf_error("not a model");
  int ns = (int) XLENGTH(subjects), threads = Rf_asInteger(s_threads);
  if (threads < 1) threads = 1;
#ifndef _OPENMP
  threads = 1;
#endif
  const zupmxode_subject **s =
      (const zupmxode_subject **) R_alloc((size_t) ns + 1, sizeof(zupmxode_subject *));
  for (int i = 0; i < ns; i++) {
    s[i] = A->subject_from_r(VECTOR_ELT(subjects, i));
    if (s[i] == NULL) Rf_error("not a subject");
  }
  const char *names[] = {"pred", "status", ""};
  SEXP ans = PROTECT(Rf_mkNamed(VECSXP, names));
  SEXP preds = Rf_allocVector(VECSXP, ns);
  SET_VECTOR_ELT(ans, 0, preds);
  SEXP status = Rf_allocVector(INTSXP, ns);
  SET_VECTOR_ELT(ans, 1, status);
  double **pp = (double **) R_alloc((size_t) ns + 1, sizeof(double *));
  double **vv = (double **) R_alloc((size_t) ns + 1, sizeof(double *));
  for (int i = 0; i < ns; i++) {
    SET_VECTOR_ELT(preds, i, Rf_allocVector(REALSXP, A->subject_nobs(s[i])));
    pp[i] = REAL(VECTOR_ELT(preds, i));
    vv[i] = (double *) R_alloc((size_t) A->subject_nobs(s[i]) + 1, sizeof(double));
  }
  SEXP held = PROTECT(Rf_allocVector(VECSXP, threads));
  zupmxode_arena **arenas = (zupmxode_arena **) R_alloc((size_t) threads, sizeof(void *));
  for (int t = 0; t < threads; t++) {
    SEXP h = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue);
    SET_VECTOR_ELT(held, t, h);
    R_RegisterCFinalizerEx(h, arena_finalize, TRUE);
    arenas[t] = A->arena_new(m, s, ns);
    if (arenas[t] == NULL) Rf_error("cannot make an arena");
    R_SetExternalPtrAddr(h, arenas[t]);
  }
  zupmxode_eval_opts o;
  opts_from(s_opts, &o, 0);
  const double *th = maybe_real(theta), *et = maybe_real(eta);
  int *st = INTEGER(status);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int i = 0; i < ns; i++) {
    int t = 0;
#ifdef _OPENMP
    t = omp_get_thread_num();
#endif
    zupmxode_eval_opts oi = o;
    oi.n_pred = (size_t) A->subject_nobs(s[i]);
    st[i] = A->subject_eval(m, s[i], th, et, pp[i], vv[i], &oi, arenas[t]);
  }
  UNPROTECT(2);
  return ans;
}

/* One call that should fail, by `what`: its status, with no R error. */
SEXP zpx_test_c_api(SEXP model, SEXP subject, SEXP small, SEXP s_what) {
  const zupmxode_api *A = zpx_api_get(ZUPMXODE_API_VERSION);
  const zupmxode_model *m = A->model_from_r(model);
  const zupmxode_subject *s = A->subject_from_r(subject);
  const zupmxode_subject *sm = A->subject_from_r(small);
  if (m == NULL || s == NULL || sm == NULL) Rf_error("not a model or subject");
  const char *what = CHAR(STRING_ELT(s_what, 0));
  int nobs = A->subject_nobs(s);
  double *pred = (double *) R_alloc((size_t) nobs + 1, sizeof(double));
  double *var = (double *) R_alloc((size_t) nobs + 1, sizeof(double));
  double *theta = (double *) R_alloc((size_t) s->ntheta + 1, sizeof(double));
  double *eta = (double *) R_alloc((size_t) s->neta + 1, sizeof(double));
  for (int i = 0; i < s->ntheta; i++) theta[i] = 0.1;
  for (int i = 0; i < s->neta; i++) eta[i] = 0;
  SEXP held = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  const zupmxode_subject *arena_for = strcmp(what, "arena") == 0 ? sm : s;
  zupmxode_arena *a = held_arena(held, m, &arena_for, 1, 0);
  zupmxode_eval_opts o;
  memset(&o, 0, sizeof(o));
  o.version = ZUPMXODE_EVAL_OPTS_VERSION;
  o.error_model = ZUPMXODE_ERR_ADD;
  o.add = 1;
  o.method = ZUPMXODE_METHOD_DEFAULT;
  o.n_pred = (size_t) nobs;
  const zupmxode_model *mm = m;
  if (strcmp(what, "version") == 0) o.version = 99;
  if (strcmp(what, "version0") == 0) o.version = 0;
  if (strcmp(what, "buffer") == 0) o.n_pred = (size_t) nobs - 1;
  if (strcmp(what, "output") == 0) o.output = 1 << 20;
  if (strcmp(what, "output_neg") == 0) o.output = -1;
  if (strcmp(what, "error_model") == 0) o.error_model = 9;
  if (strcmp(what, "sd") == 0) o.add = -1;
  if (strcmp(what, "sd_nan") == 0) o.prop = NAN;
  if (strcmp(what, "method") == 0) o.method = ZUPMXODE_METHOD_CLOSED;
  if (strcmp(what, "method_range") == 0) o.method = 7;
  int st;
  if (strcmp(what, "null") == 0) {
    st = A->subject_eval(mm, s, theta, eta, NULL, var, &o, a);
  } else if (strcmp(what, "other_model") == 0) {
    /* `small` was prepared for another model. */
    st = A->subject_eval(mm, sm, theta, eta, pred, var, &o, a);
  } else if (strcmp(what, "null_theta") == 0) {
    st = A->subject_eval(mm, s, NULL, eta, pred, var, &o, a);
  } else {
    st = A->subject_eval(mm, s, theta, eta, pred, var, &o, a);
  }
  UNPROTECT(1);
  return Rf_ScalarInteger(st);
}

SEXP zpx_status_message_r(SEXP s) { return Rf_mkString(zpx_status_message(Rf_asInteger(s))); }

void zpx_api_register(void) {
  R_RegisterCCallable("zupmxode", "zupmxode_api", (DL_FUNC) zpx_api_get);
}
