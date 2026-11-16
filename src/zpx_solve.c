/* The closed-form solve of one subject (design §4.3, §4.5, D21): the
 * scheduler's timeline walked with the state carried exactly between its
 * entries (src/core/linear.c), the model's ivl and post blocks giving the
 * rate constants and the outputs. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "core/linear.h"
#include "core/sched.h"
#include "zpx_model.h"

static SEXP fail(const char *kind, const char *message, int row) {
  const char *names[] = {"error", "message", "row", ""};
  SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 0, Rf_mkString(kind));
  SET_VECTOR_ELT(res, 1, Rf_mkString(message));
  SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(row));
  UNPROTECT(1);
  return res;
}

/* The rate constants from the captures, by the indices R passes:
 * lc = c(ncmt, depot, micro, then 7 capture indices, -1 for none) for
 * ka and either (cl, v, q, vp, q2, vp2) or (k, v, k12, k21, k13, k31). */
static int rate_constants(const int *lc, const double *cap, double *ka, double *k10,
                          double *k12, double *k21, double *k13, double *k31) {
  int ncmt = lc[0], micro = lc[2];
  const int *ix = lc + 3;
#define CAP(j) (ix[j] >= 0 ? cap[ix[j]] : NAN)
  *ka = CAP(0);
  if (micro) {
    *k10 = CAP(1);
    *k12 = CAP(3);
    *k21 = CAP(4);
    *k13 = CAP(5);
    *k31 = CAP(6);
  } else {
    double cl = CAP(1), v = CAP(2);
    *k10 = cl / v;
    *k12 = ncmt >= 2 ? CAP(3) / v : 0;
    *k21 = ncmt >= 2 ? CAP(3) / CAP(4) : 0;
    *k13 = ncmt >= 3 ? CAP(5) / v : 0;
    *k31 = ncmt >= 3 ? CAP(5) / CAP(6) : 0;
  }
#undef CAP
  return ncmt;
}

/* x: the model; records as zpx_schedule_r takes them; par: nrec x npar
 * (column-major) parameters at each record; lc: see rate_constants(). */
SEXP zpx_solve_closed(SEXP x, SEXP time, SEXP evid, SEXP amt, SEXP cmt, SEXP rate, SEXP dur,
                      SEXP ii, SEXP addl, SEXP ss, SEXP lag, SEXP f, SEXP mrate, SEXP mdur,
                      SEXP par, SEXP s_lc, SEXP init) {
  zpx_model *m = zpx_model_get(x);
  if (m == NULL) Rf_error("the model has no code: it was saved and not restored");
  int n = (int) XLENGTH(time);
  int ny = m->dims.ny, npar = m->dims.npar, nwork = m->dims.nwork;
  int ncap = m->dims.block[ZPX_POST].nout;
  if (TYPEOF(par) != REALSXP || XLENGTH(par) != (R_xlen_t) n * npar) Rf_error("bad parameters");
  if (TYPEOF(s_lc) != INTSXP || XLENGTH(s_lc) != 10) Rf_error("bad linCmt spec");
  if (TYPEOF(init) != REALSXP || XLENGTH(init) != ny) Rf_error("bad initial amounts");
  const int *lc = INTEGER(s_lc);
  if (lc[0] + lc[1] != ny) Rf_error("the model is not a closed form of its compartments");

  zpx_record *rec = (zpx_record *) R_alloc((size_t) n + 1, sizeof(zpx_record));
  for (int i = 0; i < n; i++) {
    rec[i].time = REAL(time)[i];
    rec[i].amt = REAL(amt)[i];
    rec[i].rate = REAL(rate)[i];
    rec[i].dur = REAL(dur)[i];
    rec[i].ii = REAL(ii)[i];
    rec[i].lag = REAL(lag)[i];
    rec[i].f = REAL(f)[i];
    rec[i].mrate = REAL(mrate)[i];
    rec[i].mdur = REAL(mdur)[i];
    rec[i].evid = INTEGER(evid)[i];
    rec[i].cmt = INTEGER(cmt)[i];
    rec[i].addl = INTEGER(addl)[i];
    rec[i].ss = INTEGER(ss)[i];
  }
  int error, where;
  long count = zpx_sched_count(rec, n, ny, &error, &where);
  if (count < 0) return fail("events", zpx_sched_error_message(error), where);
  zpx_entry *ent = (zpx_entry *) R_alloc((size_t) count + 1, sizeof(zpx_entry));
  count = zpx_sched_build(rec, n, ny, ent, &error, &where);
  if (count < 0) return fail("events", zpx_sched_error_message(error), where);

  long nobs = 0, nlink = 0;
  for (long i = 0; i < count; i++) {
    if (ent[i].kind == ZPX_E_OBS) nobs++;
    if (ent[i].link + 1 > nlink) nlink = ent[i].link + 1;
    if (ent[i].kind == ZPX_E_OFF) return fail("solver", "a compartment switched off (EVID 2 with a negative CMT) is not supported by the closed-form solver", ent[i].rec);
  }
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) nobs, 1 + ny + ncap));
  SEXP orow = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nobs));
  double *o = REAL(out);
  double *work = (double *) R_alloc((size_t) nwork + 1, sizeof(double));
  double *cap = (double *) R_alloc((size_t) ncap + 1, sizeof(double));
  double *active = (double *) R_alloc((size_t) nlink + 1, sizeof(double));
  double *pcur = (double *) R_alloc((size_t) npar + 1, sizeof(double));
  double y[ZPX_LIN_MAX] = {0}, r[ZPX_LIN_MAX] = {0}, zero[ZPX_LIN_MAX] = {0};
  for (long l = 0; l < nlink; l++) active[l] = 0;
  zpx_lin L;
  int cur = -1;
  double t = count > 0 ? ent[0].time : 0;
  memcpy(y, REAL(init), sizeof(double) * (size_t) ny);

  /* Parameters, work and rate constants at record i. */
#define SET_RECORD(i)                                                                 \
  do {                                                                                \
    cur = (i);                                                                        \
    for (int p = 0; p < npar; p++) pcur[p] = REAL(par)[(R_xlen_t) p * n + cur];        \
    m->lowered.fn[ZPX_IVL](t, zero, work, pcur, work);                                \
    m->lowered.fn[ZPX_POST](t, zero, cap, pcur, work);                                \
    double ka, k10, k12, k21, k13, k31;                                               \
    int nc = rate_constants(lc, cap, &ka, &k10, &k12, &k21, &k13, &k31);              \
    if (zpx_lin_setup(&L, nc, lc[1], ka, k10, k12, k21, k13, k31) != ZPX_LIN_OK) {     \
      UNPROTECT(2);                                                                   \
      return fail("solver", "a rate constant of the closed form is not a finite positive number", cur); \
    }                                                                                 \
  } while (0)

  SET_RECORD(0);
  long k = 0;
  for (long i = 0; i < count; i++) {
    const zpx_entry *e = &ent[i];
    if (e->time > t) {
      zpx_lin_step(&L, e->time - t, y, r);
      t = e->time;
    }
    if (e->phase == 2 && e->rec != cur) SET_RECORD(e->rec);
    int c = e->cmt - 1;
    switch (e->kind) {
    case ZPX_E_OBS:
      m->lowered.fn[ZPX_POST](t, y, cap, pcur, work);
      o[k] = t;
      for (int j = 0; j < ny; j++) o[k + (R_xlen_t) nobs * (1 + j)] = y[j];
      for (int j = 0; j < ncap; j++) o[k + (R_xlen_t) nobs * (1 + ny + j)] = cap[j];
      INTEGER(orow)[k] = e->rec;
      k++;
      break;
    case ZPX_E_BOLUS: y[c] += e->amount; break;
    case ZPX_E_INF_START:
      active[e->link] = e->rate;
      r[c] += e->rate;
      break;
    case ZPX_E_INF_END:
      /* An infusion a reset stopped has nothing left to end. */
      if (active[e->link] != 0) {
        r[c] -= active[e->link];
        active[e->link] = 0;
        if (fabs(r[c]) < 1e-12 * e->rate) r[c] = 0;
      }
      break;
    case ZPX_E_RESET:
      for (int j = 0; j < ny; j++) y[j] = r[j] = 0;
      for (long l = 0; l < nlink; l++) active[l] = 0;
      break;
    case ZPX_E_SS: {
      double s[ZPX_LIN_MAX] = {0}, v[ZPX_LIN_MAX] = {0};
      if (e->ss == 1) {
        for (int j = 0; j < ny; j++) y[j] = r[j] = 0;
        for (long l = 0; l < nlink; l++) active[l] = 0;
      }
      if (e->rate > 0 && e->ii == 0) {
        v[c] = e->rate;
        zpx_lin_ss_constant(&L, v, s);
      } else if (e->rate > 0) {
        v[c] = e->rate;
        zpx_lin_ss_infusion(&L, e->ii, e->dur, v, s);
      } else {
        v[c] = e->amount;
        zpx_lin_ss_bolus(&L, e->ii, v, s);
      }
      for (int j = 0; j < ny; j++) y[j] += s[j];
      break;
    }
    default: break;
    }
  }
#undef SET_RECORD
  const char *names[] = {"out", "row", ""};
  SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 0, out);
  SET_VECTOR_ELT(res, 1, orow);
  UNPROTECT(3);
  return res;
}

/* ---- The general solvers (design §4.3, D5): LSODA and Dormand-Prince ---- */

#include "core/dopri.h"
#include "vendor/liblsoda/lsoda.h"

typedef struct general {
  zpx_model *m;
  const double *par;
  double *work;
  const double *rate; /* infusion rate into each compartment */
  const int *off;     /* compartments switched off */
  int ny;
} general;

/* dy/dt: the model's rhs block, plus the infusions running, with switched-
 * off compartments held where they are (at 0). */
static void general_rhs(double t, const double *y, double *dydt, void *data) {
  general *g = data;
  g->m->lowered.fn[ZPX_RHS](t, y, dydt, g->par, g->work);
  for (int i = 0; i < g->ny; i++) {
    dydt[i] += g->rate[i];
    if (g->off[i]) dydt[i] = 0;
  }
}

static int lsoda_rhs(double t, double *y, double *dydt, void *data) {
  general_rhs(t, y, dydt, data);
  return 0;
}

/* Integrates y from t0 to t1, on `method` (1 LSODA, 2 Dormand-Prince).
 * LSODA restarts at every call: each interval begins at a discontinuity
 * or carries on after an observation, and a restart costs little against
 * the accuracy it keeps. Returns 0, or a solver error. */
static int finite_state(const double *y, int n) {
  for (int i = 0; i < n; i++)
    if (!isfinite(y[i])) return 0;
  return 1;
}

/* *reached is where the integration got to: t1, or where it failed. A
 * state that is no longer finite is a failure too: LSODA can carry a
 * solution through a blow-up without reporting one. */
static int integrate(int method, general *g, double rtol, double atol, double t0, double t1,
                     double *y, double *dwork, zpx_dopri *D, double *reached) {
  *reached = t1;
  if (!(t1 > t0)) return 0;
  if (method == 2) {
    D->n = g->ny;
    D->rtol = rtol;
    D->atol = atol;
    int st = zpx_dopri_solve(D, general_rhs, g, t0, t1, y, NULL, 0, NULL, dwork);
    if (st != 0) *reached = D->t;
    if (st == 0 && !finite_state(y, g->ny)) st = ZPX_DOPRI_E_NONFINITE;
    return st;
  }
  double *rtv = dwork + ZPX_DOPRI_WORK(g->ny), *atv = rtv + g->ny;
  for (int i = 0; i < g->ny; i++) {
    rtv[i] = rtol;
    atv[i] = atol;
  }
  struct lsoda_opt_t opt = {0};
  opt.ixpr = 0;
  opt.rtol = rtv;
  opt.atol = atv;
  opt.itask = 1;
  opt.mxstep = 100000;
  struct lsoda_context_t ctx = {0};
  ctx.function = lsoda_rhs;
  ctx.data = g;
  ctx.neq = g->ny;
  ctx.state = 1;
  if (!lsoda_prepare(&ctx, &opt)) return 1;
  double t = t0;
  lsoda(&ctx, y, &t, t1);
  /* LSODA can return a success short of tout, or carry on through a
   * blow-up; both are failures here. */
  int ok = ctx.state > 0 && t >= t1 && finite_state(y, g->ny);
  lsoda_free(&ctx);
  if (!ok) *reached = t;
  return ok ? 0 : 1;
}

/* Like zpx_solve_closed(), for any model: s_opts = c(method, rtol, atol,
 * ss_tol, ss_cap). Returns list(out, row, warnings): a solver failure
 * leaves the subject's outputs NA from where it happened, with a warning;
 * a steady state that does not converge within ss_cap intervals, a
 * warning. */
SEXP zpx_solve_general(SEXP x, SEXP time, SEXP evid, SEXP amt, SEXP cmt, SEXP rate, SEXP dur,
                       SEXP ii, SEXP addl, SEXP ss, SEXP lag, SEXP f, SEXP mrate, SEXP mdur,
                       SEXP par, SEXP init, SEXP s_opts) {
  zpx_model *m = zpx_model_get(x);
  if (m == NULL) Rf_error("the model has no code: it was saved and not restored");
  int n = (int) XLENGTH(time);
  int ny = m->dims.ny, npar = m->dims.npar, nwork = m->dims.nwork;
  int ncap = m->dims.block[ZPX_POST].nout;
  if (TYPEOF(par) != REALSXP || XLENGTH(par) != (R_xlen_t) n * npar) Rf_error("bad parameters");
  if (TYPEOF(init) != REALSXP || XLENGTH(init) != ny) Rf_error("bad initial amounts");
  if (TYPEOF(s_opts) != REALSXP || XLENGTH(s_opts) != 5) Rf_error("bad options");
  int method = (int) REAL(s_opts)[0];
  double rtol = REAL(s_opts)[1], atol = REAL(s_opts)[2], ss_tol = REAL(s_opts)[3];
  long ss_cap = (long) REAL(s_opts)[4];

  zpx_record *rec = (zpx_record *) R_alloc((size_t) n + 1, sizeof(zpx_record));
  for (int i = 0; i < n; i++) {
    rec[i].time = REAL(time)[i];
    rec[i].amt = REAL(amt)[i];
    rec[i].rate = REAL(rate)[i];
    rec[i].dur = REAL(dur)[i];
    rec[i].ii = REAL(ii)[i];
    rec[i].lag = REAL(lag)[i];
    rec[i].f = REAL(f)[i];
    rec[i].mrate = REAL(mrate)[i];
    rec[i].mdur = REAL(mdur)[i];
    rec[i].evid = INTEGER(evid)[i];
    rec[i].cmt = INTEGER(cmt)[i];
    rec[i].addl = INTEGER(addl)[i];
    rec[i].ss = INTEGER(ss)[i];
  }
  int error, where;
  long count = zpx_sched_count(rec, n, ny, &error, &where);
  if (count < 0) return fail("events", zpx_sched_error_message(error), where);
  zpx_entry *ent = (zpx_entry *) R_alloc((size_t) count + 1, sizeof(zpx_entry));
  count = zpx_sched_build(rec, n, ny, ent, &error, &where);
  if (count < 0) return fail("events", zpx_sched_error_message(error), where);

  long nobs = 0, nlink = 0;
  for (long i = 0; i < count; i++) {
    if (ent[i].kind == ZPX_E_OBS) nobs++;
    if (ent[i].link + 1 > nlink) nlink = ent[i].link + 1;
  }
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) nobs, 1 + ny + ncap));
  SEXP orow = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nobs));
  double *o = REAL(out);
  double *work = (double *) R_alloc((size_t) nwork + 1, sizeof(double));
  double *cap = (double *) R_alloc((size_t) ncap + 1, sizeof(double));
  double *active = (double *) R_alloc((size_t) nlink + 1, sizeof(double));
  double *pcur = (double *) R_alloc((size_t) npar + 1, sizeof(double));
  double *y = (double *) R_alloc((size_t) ny + 1, sizeof(double));
  double *r = (double *) R_alloc((size_t) ny + 1, sizeof(double));
  double *s = (double *) R_alloc((size_t) ny + 1, sizeof(double));
  double *prev = (double *) R_alloc((size_t) ny + 1, sizeof(double));
  double *zero = (double *) R_alloc((size_t) ny + 1, sizeof(double));
  double *dwork = (double *) R_alloc((size_t) ZPX_DOPRI_WORK(ny) + 2 * (size_t) ny + 1, sizeof(double));
  int *off = (int *) R_alloc((size_t) ny + 1, sizeof(int));
  for (int j = 0; j < ny; j++) {
    r[j] = zero[j] = 0;
    off[j] = 0;
  }
  for (long l = 0; l < nlink; l++) active[l] = 0;
  memcpy(y, REAL(init), sizeof(double) * (size_t) ny);
  general g = {m, pcur, work, r, off, ny};
  zpx_dopri D = {ny, rtol, atol, 1000000, 0, 0};
  int cur = -1, failed = 0, failed_row = -1, ss_warn = 0, ss_warn_row = -1;
  double t = count > 0 ? ent[0].time : 0, failed_at = NAN;

#define SET_RECORD(i)                                                          \
  do {                                                                         \
    cur = (i);                                                                 \
    for (int p = 0; p < npar; p++) pcur[p] = REAL(par)[(R_xlen_t) p * n + cur]; \
    m->lowered.fn[ZPX_IVL](t, zero, work, pcur, work);                         \
    D.h = 0;                                                                   \
  } while (0)

  /* Dormand-Prince's dense output: a run of observations with nothing
   * between them that changes the system is one integration, to the last
   * of them, the others read from the continuous extension. */
  double *dense_t = (double *) R_alloc((size_t) count + 1, sizeof(double));
  double *dense_y = (double *) R_alloc(((size_t) count + 1) * (size_t) ny, sizeof(double));
  long dense_from = -1, dense_to = -1;

  SET_RECORD(0);
  long k = 0;
  for (long i = 0; i < count; i++) {
    const zpx_entry *e = &ent[i];
    int dense_hit = method == 2 && i > dense_from && i <= dense_to;
    if (!failed && e->time > t && !dense_hit) {
      long last = i;
      if (method == 2 && e->kind == ZPX_E_OBS) {
        /* The run: observations of one record's parameters, in time. */
        while (last + 1 < count && ent[last + 1].kind == ZPX_E_OBS &&
               ent[last + 1].rec == ent[last].rec + 1 && ent[last + 1].phase == 2)
          last++;
        /* Every observation is its own record, so a parameter change can
         * come with each; the run is cut where the parameters differ. */
        long cut = i;
        while (cut + 1 <= last) {
          int same = 1;
          for (int p = 0; p < npar && same; p++)
            same = REAL(par)[(R_xlen_t) p * n + ent[cut + 1].rec] ==
                   REAL(par)[(R_xlen_t) p * n + cur];
          if (!same) break;
          cut++;
        }
        last = cut;
      }
      int st;
      double reached;
      if (last > i) {
        int nd = (int) (last - i);
        for (int q = 0; q < nd; q++) dense_t[q] = ent[i + q].time;
        D.n = ny;
        D.rtol = rtol;
        D.atol = atol;
        st = zpx_dopri_solve(&D, general_rhs, &g, t, ent[last].time, y, dense_t, nd, dense_y,
                             dwork);
        if (st == 0 && !finite_state(y, ny)) st = ZPX_DOPRI_E_NONFINITE;
        reached = st != 0 ? D.t : ent[last].time;
        dense_from = i;
        /* After a failure only the outputs written before it are kept. */
        dense_to = st != 0 ? i + D.nout_done : last;
      } else {
        st = integrate(method, &g, rtol, atol, t, e->time, y, dwork, &D, &reached);
      }
      if (st != 0) {
        failed = 1;
        failed_at = reached;
        failed_row = e->rec;
      }
    }
    t = e->time;
    if (e->phase == 2 && e->rec != cur) SET_RECORD(e->rec);
    if (e->kind == ZPX_E_OBS && method == 2 && i >= dense_from && i < dense_to) {
      /* An observation inside the run: its state from the dense output. */
      double *yd = dense_y + (i - dense_from) * (long) ny;
      m->lowered.fn[ZPX_POST](t, yd, cap, pcur, work);
      o[k] = t;
      for (int j = 0; j < ny; j++) o[k + (R_xlen_t) nobs * (1 + j)] = yd[j];
      for (int j = 0; j < ncap; j++) o[k + (R_xlen_t) nobs * (1 + ny + j)] = cap[j];
      INTEGER(orow)[k] = e->rec;
      k++;
      continue;
    }
    int c = e->cmt - 1;
    switch (e->kind) {
    case ZPX_E_OBS:
      o[k] = t;
      if (failed) {
        for (int j = 0; j < ny + ncap; j++) o[k + (R_xlen_t) nobs * (1 + j)] = NA_REAL;
      } else {
        m->lowered.fn[ZPX_POST](t, y, cap, pcur, work);
        for (int j = 0; j < ny; j++) o[k + (R_xlen_t) nobs * (1 + j)] = y[j];
        for (int j = 0; j < ncap; j++) o[k + (R_xlen_t) nobs * (1 + ny + j)] = cap[j];
      }
      INTEGER(orow)[k] = e->rec;
      k++;
      break;
    case ZPX_E_BOLUS:
      y[c] += e->amount;
      off[c] = 0; /* a dose turns its compartment on (PREDPP) */
      D.h = 0;
      break;
    case ZPX_E_INF_START:
      active[e->link] = e->rate;
      r[c] += e->rate;
      off[c] = 0;
      D.h = 0;
      break;
    case ZPX_E_INF_END:
      if (active[e->link] != 0) {
        r[c] -= active[e->link];
        active[e->link] = 0;
        if (fabs(r[c]) < 1e-12 * e->rate) r[c] = 0;
      }
      D.h = 0;
      break;
    case ZPX_E_RESET:
      for (int j = 0; j < ny; j++) {
        y[j] = r[j] = 0;
        off[j] = 0;
      }
      for (long l = 0; l < nlink; l++) active[l] = 0;
      D.h = 0;
      break;
    case ZPX_E_OFF:
      y[c] = 0;
      off[c] = 1;
      D.h = 0;
      break;
    case ZPX_E_ON: off[c] = 0; break;
    case ZPX_E_SS: {
      /* Repeated dosing from empty, until the state before a dose moves by
       * less than ss_tol (relative, with atol); at most ss_cap intervals. */
      for (int j = 0; j < ny; j++) s[j] = 0;
      double *saved_r = (double *) R_alloc((size_t) ny + 1, sizeof(double));
      memcpy(saved_r, r, sizeof(double) * (size_t) ny);
      int converged = 0;
      double ss_reached;
      for (long it = 0; it < ss_cap && !failed; it++) {
        memcpy(prev, s, sizeof(double) * (size_t) ny);
        for (int j = 0; j < ny; j++) r[j] = 0;
        if (e->rate > 0 && e->ii == 0) {
          /* A constant infusion: 24 time units at a time until nothing
           * moves, under the same test and cap. */
          r[c] = e->rate;
          if (integrate(method, &g, rtol, atol, 0, 24, s, dwork, &D, &ss_reached) != 0) failed = 1;
        } else if (e->rate > 0) {
          r[c] = e->rate;
          if (integrate(method, &g, rtol, atol, 0, e->dur, s, dwork, &D, &ss_reached) != 0) failed = 1;
          r[c] = 0;
          if (!failed && integrate(method, &g, rtol, atol, e->dur, e->ii, s, dwork, &D, &ss_reached) != 0)
            failed = 1;
        } else {
          s[c] += e->amount;
          if (integrate(method, &g, rtol, atol, 0, e->ii, s, dwork, &D, &ss_reached) != 0) failed = 1;
        }
        double change = 0;
        for (int j = 0; j < ny; j++) {
          double d = fabs(s[j] - prev[j]) / (atol + ss_tol * fmax(fabs(s[j]), fabs(prev[j])));
          if (d > change) change = d;
        }
        if (it > 0 && change <= 1) {
          converged = 1;
          break;
        }
      }
      memcpy(r, saved_r, sizeof(double) * (size_t) ny);
      if (failed) {
        failed_at = t;
        failed_row = e->rec;
        break;
      }
      if (!converged) {
        ss_warn++;
        ss_warn_row = e->rec;
      }
      if (e->ss == 1) {
        for (int j = 0; j < ny; j++) {
          y[j] = r[j] = 0;
          off[j] = 0;
        }
        for (long l = 0; l < nlink; l++) active[l] = 0;
      }
      for (int j = 0; j < ny; j++) y[j] += s[j];
      D.h = 0;
      break;
    }
    default: break;
    }
  }
#undef SET_RECORD
  const char *names[] = {"out", "row", "failed_at", "failed_row", "ss_unconverged", "ss_row", ""};
  SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 0, out);
  SET_VECTOR_ELT(res, 1, orow);
  SET_VECTOR_ELT(res, 2, Rf_ScalarReal(failed ? failed_at : NA_REAL));
  SET_VECTOR_ELT(res, 3, Rf_ScalarInteger(failed ? failed_row : -1));
  SET_VECTOR_ELT(res, 4, Rf_ScalarInteger(ss_warn));
  SET_VECTOR_ELT(res, 5, Rf_ScalarInteger(ss_warn_row));
  UNPROTECT(3);
  return res;
}

/* ---- The linear solver (design §4.3): the matrix exponential ------------ */

#include "core/expm.h"

#define ZPX_LINEAR_MAX_STATES (ZPX_EXPM_MAX - 1)

typedef struct linear {
  int n;
  double a[ZPX_LINEAR_MAX_STATES * ZPX_LINEAR_MAX_STATES]; /* A, row-major */
  double c[ZPX_LINEAR_MAX_STATES];                         /* the constant input */
  double mwork[8 * ZPX_EXPM_MAX * ZPX_EXPM_MAX];
  double M[ZPX_EXPM_MAX * ZPX_EXPM_MAX], X[ZPX_EXPM_MAX * ZPX_EXPM_MAX];
} linear;

/* A and c of dy/dt = A y + c from the model's rhs block: c at y = 0, each
 * column of A as rhs(e_j) - c. Exact for a linear model, whose rhs is
 * affine in y, up to rounding. */
static void linear_probe(linear *L, zpx_model *m, double t, const double *par, double *work,
                         double *ybuf, double *dbuf) {
  int n = L->n;
  for (int i = 0; i < n; i++) ybuf[i] = 0;
  m->lowered.fn[ZPX_RHS](t, ybuf, L->c, par, work);
  for (int j = 0; j < n; j++) {
    ybuf[j] = 1;
    m->lowered.fn[ZPX_RHS](t, ybuf, dbuf, par, work);
    for (int i = 0; i < n; i++) L->a[i * n + j] = dbuf[i] - L->c[i];
    ybuf[j] = 0;
  }
}

/* exp of [[A dt, (c + r) dt], [0, 0]] into L->X: E(dt) in its top-left
 * block, G(dt)(c + r) in its last column. r may be NULL; with_c says
 * whether the constant input is included. */
static int linear_exp(linear *L, double dt, const double *r, int with_c) {
  int n = L->n, m = n + 1;
  memset(L->M, 0, sizeof(double) * (size_t) (m * m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) L->M[i * m + j] = L->a[i * n + j] * dt;
    L->M[i * m + n] = ((with_c ? L->c[i] : 0) + (r != NULL ? r[i] : 0)) * dt;
  }
  return zpx_expm(m, L->M, L->X, L->mwork);
}

/* y <- E(dt) y + G(dt)(c + r). */
static int linear_step(linear *L, double dt, double *y, const double *r) {
  int n = L->n, m = n + 1;
  if (!(dt > 0)) return 0;
  if (linear_exp(L, dt, r, 1) != 0) return 1;
  double out[ZPX_LINEAR_MAX_STATES];
  for (int i = 0; i < n; i++) {
    double s = L->X[i * m + n];
    for (int j = 0; j < n; j++) s += L->X[i * m + j] * y[j];
    out[i] = s;
  }
  memcpy(y, out, sizeof(double) * (size_t) n);
  return 0;
}

/* Solves (I - E(ii)) y = b for y, E from L->X after linear_exp(ii). */
static int linear_ss_solve(linear *L, double *b, double *y) {
  int n = L->n, m = n + 1;
  double M[ZPX_LINEAR_MAX_STATES * ZPX_LINEAR_MAX_STATES];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) M[i * n + j] = (i == j ? 1.0 : 0.0) - L->X[i * m + j];
  /* Gaussian elimination with partial pivoting on M y = b. */
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int r = c + 1; r < n; r++)
      if (fabs(M[r * n + c]) > fabs(M[p * n + c])) p = r;
    if (M[p * n + c] == 0) return 1;
    if (p != c) {
      for (int k = 0; k < n; k++) {
        double t = M[c * n + k];
        M[c * n + k] = M[p * n + k];
        M[p * n + k] = t;
      }
      double t = b[c];
      b[c] = b[p];
      b[p] = t;
    }
    for (int r = c + 1; r < n; r++) {
      double f = M[r * n + c] / M[c * n + c];
      for (int k = c; k < n; k++) M[r * n + k] -= f * M[c * n + k];
      b[r] -= f * b[c];
    }
  }
  for (int c = n - 1; c >= 0; c--) {
    double s = b[c];
    for (int k = c + 1; k < n; k++) s -= M[c * n + k] * y[k];
    y[c] = s / M[c * n + c];
  }
  return 0;
}

/* Steady state just before a dose (design §4.2), with the constant input
 * c running throughout: a bolus b every ii, an infusion r for dur every
 * ii, or (ii = 0) a constant infusion r. */
static int linear_ss(linear *L, int kind, double ii, double dur, const double *v, double *y) {
  int n = L->n, m = n + 1;
  double b[ZPX_LINEAR_MAX_STATES] = {0};
  if (kind == 2) {
    /* A y + c + r = 0. */
    for (int i = 0; i < n; i++) b[i] = -(L->c[i] + v[i]);
    /* Reuse the solver on (I - E) with E = I + A: I - E = -A. */
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) L->X[i * m + j] = (i == j ? 1.0 : 0.0) + L->a[i * n + j];
    for (int i = 0; i < n; i++) b[i] = -b[i];
    return linear_ss_solve(L, b, y);
  }
  if (kind == 0) {
    /* y = E (y + v) + G c  =>  (I - E) y = E v + G c. */
    if (linear_exp(L, ii, NULL, 1) != 0) return 1;
    for (int i = 0; i < n; i++) {
      double s = L->X[i * m + n];
      for (int j = 0; j < n; j++) s += L->X[i * m + j] * v[j];
      b[i] = s;
    }
    return linear_ss_solve(L, b, y);
  }
  /* Infusion: from a trough y, dur with r then ii - dur without:
   * y = E(ii) y + [E(ii - dur) G(dur)(c + r) + G(ii - dur) c]. */
  double z[ZPX_LINEAR_MAX_STATES] = {0};
  if (linear_step(L, dur, z, v) != 0) return 1;
  double zero[ZPX_LINEAR_MAX_STATES] = {0};
  if (linear_step(L, ii - dur, z, zero) != 0) return 1;
  memcpy(b, z, sizeof(double) * (size_t) n);
  if (linear_exp(L, ii, NULL, 0) != 0) return 1;
  return linear_ss_solve(L, b, y);
}

SEXP zpx_solve_linear(SEXP x, SEXP time, SEXP evid, SEXP amt, SEXP cmt, SEXP rate, SEXP dur,
                      SEXP ii, SEXP addl, SEXP ss, SEXP lag, SEXP f, SEXP mrate, SEXP mdur,
                      SEXP par, SEXP init) {
  zpx_model *m = zpx_model_get(x);
  if (m == NULL) Rf_error("the model has no code: it was saved and not restored");
  int n = (int) XLENGTH(time);
  int ny = m->dims.ny, npar = m->dims.npar, nwork = m->dims.nwork;
  int ncap = m->dims.block[ZPX_POST].nout;
  if (ny > ZPX_LINEAR_MAX_STATES) Rf_error("too many states for the linear solver");
  if (TYPEOF(par) != REALSXP || XLENGTH(par) != (R_xlen_t) n * npar) Rf_error("bad parameters");
  if (TYPEOF(init) != REALSXP || XLENGTH(init) != ny) Rf_error("bad initial amounts");
  zpx_record *rec = (zpx_record *) R_alloc((size_t) n + 1, sizeof(zpx_record));
  for (int i = 0; i < n; i++) {
    rec[i].time = REAL(time)[i];
    rec[i].amt = REAL(amt)[i];
    rec[i].rate = REAL(rate)[i];
    rec[i].dur = REAL(dur)[i];
    rec[i].ii = REAL(ii)[i];
    rec[i].lag = REAL(lag)[i];
    rec[i].f = REAL(f)[i];
    rec[i].mrate = REAL(mrate)[i];
    rec[i].mdur = REAL(mdur)[i];
    rec[i].evid = INTEGER(evid)[i];
    rec[i].cmt = INTEGER(cmt)[i];
    rec[i].addl = INTEGER(addl)[i];
    rec[i].ss = INTEGER(ss)[i];
  }
  int error, where;
  long count = zpx_sched_count(rec, n, ny, &error, &where);
  if (count < 0) return fail("events", zpx_sched_error_message(error), where);
  zpx_entry *ent = (zpx_entry *) R_alloc((size_t) count + 1, sizeof(zpx_entry));
  count = zpx_sched_build(rec, n, ny, ent, &error, &where);
  if (count < 0) return fail("events", zpx_sched_error_message(error), where);
  long nobs = 0, nlink = 0;
  for (long i = 0; i < count; i++) {
    if (ent[i].kind == ZPX_E_OBS) nobs++;
    if (ent[i].link + 1 > nlink) nlink = ent[i].link + 1;
    if (ent[i].kind == ZPX_E_OFF)
      return fail("solver", "a compartment switched off is not supported by the linear solver", ent[i].rec);
  }
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) nobs, 1 + ny + ncap));
  SEXP orow = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nobs));
  double *o = REAL(out);
  double *work = (double *) R_alloc((size_t) nwork + 1, sizeof(double));
  double *cap = (double *) R_alloc((size_t) ncap + 1, sizeof(double));
  double *active = (double *) R_alloc((size_t) nlink + 1, sizeof(double));
  double *pcur = (double *) R_alloc((size_t) npar + 1, sizeof(double));
  linear *L = (linear *) R_alloc(1, sizeof(linear));
  L->n = ny;
  double y[ZPX_LINEAR_MAX_STATES] = {0}, r[ZPX_LINEAR_MAX_STATES] = {0},
         zero[ZPX_LINEAR_MAX_STATES] = {0}, ybuf[ZPX_LINEAR_MAX_STATES],
         dbuf[ZPX_LINEAR_MAX_STATES];
  for (long l = 0; l < nlink; l++) active[l] = 0;
  memcpy(y, REAL(init), sizeof(double) * (size_t) ny);
  int cur = -1;
  double t = count > 0 ? ent[0].time : 0;
#define SET_RECORD(i)                                                          \
  do {                                                                         \
    cur = (i);                                                                 \
    for (int p = 0; p < npar; p++) pcur[p] = REAL(par)[(R_xlen_t) p * n + cur]; \
    m->lowered.fn[ZPX_IVL](t, zero, work, pcur, work);                         \
    linear_probe(L, m, t, pcur, work, ybuf, dbuf);                             \
  } while (0)
  SET_RECORD(0);
  long k = 0;
  for (long i = 0; i < count; i++) {
    const zpx_entry *e = &ent[i];
    if (e->time > t) {
      if (linear_step(L, e->time - t, y, r) != 0) {
        UNPROTECT(2);
        return fail("solver", "the linear system's matrix is not finite", e->rec);
      }
      t = e->time;
    }
    if (e->phase == 2 && e->rec != cur) SET_RECORD(e->rec);
    int c = e->cmt - 1;
    switch (e->kind) {
    case ZPX_E_OBS:
      m->lowered.fn[ZPX_POST](t, y, cap, pcur, work);
      o[k] = t;
      for (int j = 0; j < ny; j++) o[k + (R_xlen_t) nobs * (1 + j)] = y[j];
      for (int j = 0; j < ncap; j++) o[k + (R_xlen_t) nobs * (1 + ny + j)] = cap[j];
      INTEGER(orow)[k] = e->rec;
      k++;
      break;
    case ZPX_E_BOLUS: y[c] += e->amount; break;
    case ZPX_E_INF_START:
      active[e->link] = e->rate;
      r[c] += e->rate;
      break;
    case ZPX_E_INF_END:
      if (active[e->link] != 0) {
        r[c] -= active[e->link];
        active[e->link] = 0;
        if (fabs(r[c]) < 1e-12 * e->rate) r[c] = 0;
      }
      break;
    case ZPX_E_RESET:
      for (int j = 0; j < ny; j++) y[j] = r[j] = 0;
      for (long l = 0; l < nlink; l++) active[l] = 0;
      break;
    case ZPX_E_SS: {
      double s[ZPX_LINEAR_MAX_STATES] = {0}, v[ZPX_LINEAR_MAX_STATES] = {0};
      if (e->ss == 1) {
        for (int j = 0; j < ny; j++) y[j] = r[j] = 0;
        for (long l = 0; l < nlink; l++) active[l] = 0;
      }
      int kind = e->rate > 0 ? (e->ii == 0 ? 2 : 1) : 0;
      v[c] = kind == 0 ? e->amount : e->rate;
      /* SS = 2 adds the dose's own steady state to what is there, which
       * already carries the constant input's share: without c. */
      double saved_c[ZPX_LINEAR_MAX_STATES];
      memcpy(saved_c, L->c, sizeof(double) * (size_t) ny);
      if (e->ss == 2)
        for (int j = 0; j < ny; j++) L->c[j] = 0;
      int bad = linear_ss(L, kind, e->ii, e->dur, v, s);
      memcpy(L->c, saved_c, sizeof(double) * (size_t) ny);
      if (bad) {
        UNPROTECT(2);
        return fail("solver", "the linear system has no steady state", e->rec);
      }
      for (int j = 0; j < ny; j++) y[j] = (e->ss == 1 ? 0 : y[j]) + s[j];
      break;
    }
    default: break;
    }
  }
#undef SET_RECORD
  const char *names[] = {"out", "row", ""};
  SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(res, 0, out);
  SET_VECTOR_ELT(res, 1, orow);
  UNPROTECT(3);
  return res;
}
