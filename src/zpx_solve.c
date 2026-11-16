/* The parallel runtime (design §4.4): every subject of a solve, split
 * across OpenMP threads, each subject solved by src/core/solve.c with an
 * arena allocated here, before the parallel region. Inside the region
 * nothing touches an R object or the R heap; R_CheckUserInterrupt() runs
 * between chunks of subjects, on the main thread, outside it. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/solve.h"
#include "zpx_lower.h"
#include "zpx_model.h"

/* The processors OpenMP sees, or 0 for a build without OpenMP. */
SEXP zpx_omp_max(void) {
#ifdef _OPENMP
  return Rf_ScalarInteger(omp_get_num_procs());
#else
  return Rf_ScalarInteger(0);
#endif
}

/* The interpreter tier's extra contexts, one per thread past the first
 * (zumir D10), owned by an external pointer so an interrupt between
 * chunks cannot strand them. */
typedef struct copies {
  int n;
  zpx_lowered *l;
} copies;

static void copies_finalize(SEXP x) {
  copies *c = R_ExternalPtrAddr(x);
  R_ClearExternalPtr(x);
  if (c == NULL) return;
  for (int i = 0; i < c->n; i++) zpx_lower_free(&c->l[i]);
  free(c->l);
  free(c);
}

static void set_blocks(zpx_blocks *b, const zpx_lowered *l, const zpx_model *m) {
  for (int k = 0; k < ZPX_IR_NBLOCKS; k++) b->fn[k] = (zpx_block) l->fn[k];
  b->ny = m->dims.ny;
  b->npar = m->dims.npar;
  b->nwork = m->dims.nwork;
  b->ncap = m->dims.block[ZPX_POST].nout;
}

/* x: the model; code, consts: its IR (lowered again per thread on the
 * interpreter tier); the record columns of every subject, subject s's
 * records at [start[s], start[s + 1]); par: nrec x npar parameters, the
 * subject's slice of each column; method: per subject; s_opts = c(rtol,
 * atol, ss_tol, ss_cap); lc: the closed form's spec; threads. */
SEXP zpx_solve_all(SEXP x, SEXP code, SEXP consts, SEXP time, SEXP evid, SEXP amt, SEXP cmt,
                   SEXP rate, SEXP dur, SEXP ii, SEXP addl, SEXP ss, SEXP start, SEXP par,
                   SEXP method, SEXP s_opts, SEXP s_lc, SEXP s_threads) {
  zpx_model *m = zpx_model_get(x);
  if (m == NULL) Rf_error("the model has no code: it was saved and not restored");
  R_xlen_t nrec = XLENGTH(time);
  int nsub = (int) XLENGTH(start) - 1;
  int ny = m->dims.ny, npar = m->dims.npar, ncap = m->dims.block[ZPX_POST].nout;
  if (TYPEOF(par) != REALSXP || XLENGTH(par) != nrec * npar) Rf_error("bad parameters");
  if (TYPEOF(s_lc) != INTSXP || XLENGTH(s_lc) != 10) Rf_error("bad linCmt spec");
  if (TYPEOF(s_opts) != REALSXP || XLENGTH(s_opts) != 4) Rf_error("bad options");
  if (TYPEOF(method) != INTSXP || XLENGTH(method) != nsub) Rf_error("bad methods");
  const int *st = INTEGER(start);
  int threads = Rf_asInteger(s_threads);
  if (threads < 1) threads = 1;
#ifndef _OPENMP
  threads = 1;
#endif
  if (threads > nsub && nsub > 0) threads = nsub;

  /* The records, in one array; each subject's slice is its own. */
  zpx_record *rec = (zpx_record *) R_alloc((size_t) nrec + 1, sizeof(zpx_record));
  for (R_xlen_t i = 0; i < nrec; i++) {
    rec[i].time = REAL(time)[i];
    rec[i].amt = REAL(amt)[i];
    rec[i].rate = REAL(rate)[i];
    rec[i].dur = REAL(dur)[i];
    rec[i].ii = REAL(ii)[i];
    rec[i].evid = INTEGER(evid)[i];
    rec[i].cmt = INTEGER(cmt)[i];
    rec[i].addl = INTEGER(addl)[i];
    rec[i].ss = INTEGER(ss)[i];
    /* Placeholders for the count: the lags and rates pre gives later do
     * not change how many entries a record makes. */
    rec[i].lag = 0;
    rec[i].f = 1;
    rec[i].mrate = 1;
    rec[i].mdur = 1;
  }

  /* A serial pass: each subject's observations (its rows of the output)
   * and its entries (the arena it needs). A count that fails is the
   * subject's error, reported by its solve. */
  long *row0 = (long *) R_alloc((size_t) nsub + 1, sizeof(long));
  long nobs = 0, max_entries = 0;
  int max_n = 0;
  for (int s = 0; s < nsub; s++) {
    row0[s] = nobs;
    int n = st[s + 1] - st[s];
    if (n > max_n) max_n = n;
    for (int i = st[s]; i < st[s + 1]; i++) nobs += rec[i].evid == 0;
    int error, where;
    long c = zpx_sched_count(rec + st[s], n, ny, &error, &where);
    if (c > max_entries) max_entries = c;
  }
  row0[nsub] = nobs;

  zpx_blocks *blocks = (zpx_blocks *) R_alloc((size_t) threads, sizeof(zpx_blocks));
  set_blocks(&blocks[0], &m->lowered, m);
  SEXP held = PROTECT(R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(held, copies_finalize, TRUE);
  if (m->lowered.tier == ZPX_TIER_INTERP && threads > 1) {
    zpx_ir ir;
    int bl, in;
    if (zpx_ir_verify(INTEGER(code), (size_t) XLENGTH(code), REAL(consts),
                      (size_t) XLENGTH(consts), &ir, &bl, &in) != ZPX_IR_OK)
      Rf_error("the model's IR does not verify");
    copies *cp = calloc(1, sizeof(copies));
    if (cp == NULL) Rf_error("cannot allocate the interpreter's contexts");
    R_SetExternalPtrAddr(held, cp);
    cp->l = calloc((size_t) threads - 1, sizeof(zpx_lowered));
    if (cp->l == NULL) Rf_error("cannot allocate the interpreter's contexts");
    void *scratch = R_alloc(zpx_lower_scratch_size(&ir) + 1, 1);
    for (int t = 1; t < threads; t++) {
      const char *msg;
      int where;
      cp->n = t;
      if (zpx_lower(&ir, ZPX_TIER_INTERP, &cp->l[t - 1], scratch, &msg, &where) != ZPX_LOWER_OK)
        Rf_error("lowering a thread's copy of the model failed%s%s", msg ? ": " : "",
                 msg ? msg : "");
      set_blocks(&blocks[t], &cp->l[t - 1], m);
    }
  } else {
    for (int t = 1; t < threads; t++) blocks[t] = blocks[0];
  }

  zpx_solve_opts o;
  memset(&o, 0, sizeof(o));
  o.rtol = REAL(s_opts)[0];
  o.atol = REAL(s_opts)[1];
  o.ss_tol = REAL(s_opts)[2];
  o.ss_cap = (long) REAL(s_opts)[3];
  memcpy(o.lc, INTEGER(s_lc), sizeof(o.lc));

  size_t arena_size = zpx_solve_arena_size(&blocks[0], max_n, max_entries);
  char **arenas = (char **) R_alloc((size_t) threads, sizeof(char *));
  for (int t = 0; t < threads; t++) arenas[t] = R_alloc(arena_size + 8, 1);

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) nobs, 1 + ny + ncap));
  SEXP orow = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nobs));
  zpx_subject_result *res = (zpx_subject_result *) R_alloc((size_t) nsub + 1, sizeof(zpx_subject_result));
  double *o_out = REAL(out);
  int *o_row = INTEGER(orow);
  for (R_xlen_t i = 0; i < XLENGTH(out); i++) o_out[i] = NA_REAL;
  for (long i = 0; i < nobs; i++) o_row[i] = -1;
  const double *p = REAL(par);
  const int *meth = INTEGER(method);

  /* Chunks of subjects: the parallel region inside, the interrupt check
   * between them. */
  int chunk = threads * 8;
  for (int s0 = 0; s0 < nsub; s0 += chunk) {
    int s1 = s0 + chunk < nsub ? s0 + chunk : nsub;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (int s = s0; s < s1; s++) {
      int t = 0;
#ifdef _OPENMP
      t = omp_get_thread_num();
#endif
      int n = st[s + 1] - st[s];
      zpx_solve_opts os = o;
      os.method = meth[s];
      /* The subject's parameters: its slice of each column, read where
       * the full matrix holds them. */
      zpx_solve_subject(&blocks[t], &os, rec + st[s], n, p + st[s], (long) nrec, o_out, nobs, o_row, row0[s],
                        arenas[t], arena_size, &res[s]);
    }
    R_CheckUserInterrupt();
  }

  const char *names[] = {"out", "row", "error", "message", "error_row", "failed_at",
                         "failed_row", "ss_unconverged", "ss_row", ""};
  SEXP ans = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(ans, 0, out);
  SET_VECTOR_ELT(ans, 1, orow);
  SEXP er = Rf_allocVector(INTSXP, nsub);
  SET_VECTOR_ELT(ans, 2, er);
  SEXP msg = Rf_allocVector(STRSXP, nsub);
  SET_VECTOR_ELT(ans, 3, msg);
  SEXP erow = Rf_allocVector(INTSXP, nsub);
  SET_VECTOR_ELT(ans, 4, erow);
  SEXP fat = Rf_allocVector(REALSXP, nsub);
  SET_VECTOR_ELT(ans, 5, fat);
  SEXP frow = Rf_allocVector(INTSXP, nsub);
  SET_VECTOR_ELT(ans, 6, frow);
  SEXP ssu = Rf_allocVector(INTSXP, nsub);
  SET_VECTOR_ELT(ans, 7, ssu);
  SEXP ssr = Rf_allocVector(INTSXP, nsub);
  SET_VECTOR_ELT(ans, 8, ssr);
  for (int s = 0; s < nsub; s++) {
    INTEGER(er)[s] = res[s].error;
    SET_STRING_ELT(msg, s, Rf_mkChar(res[s].message != NULL ? res[s].message : ""));
    INTEGER(erow)[s] = res[s].error_row;
    REAL(fat)[s] = res[s].failed_at;
    INTEGER(frow)[s] = res[s].failed_row;
    INTEGER(ssu)[s] = res[s].ss_unconverged;
    INTEGER(ssr)[s] = res[s].ss_row;
  }
  UNPROTECT(4);
  return ans;
}
