/* One subject's solve (src/core/solve.h). Builds without R. */
#include "solve.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "dopri.h"
#include "expm.h"
#include "linear.h"
#include "../vendor/liblsoda/lsoda.h"

/* ---- The arena ---------------------------------------------------------- */

typedef struct lin_sys { /* the matrix exponential's system */
  int n;
  double a[ZPX_SOLVE_LINEAR_MAX * ZPX_SOLVE_LINEAR_MAX];
  double c[ZPX_SOLVE_LINEAR_MAX];
  double mwork[8 * ZPX_EXPM_MAX * ZPX_EXPM_MAX];
  double M[ZPX_EXPM_MAX * ZPX_EXPM_MAX], X[ZPX_EXPM_MAX * ZPX_EXPM_MAX];
} lin_sys;

typedef struct arena {
  zpx_entry *ent;
  double *y, *r, *s, *prev, *zero, *work, *cap, *pcur, *dwork, *active, *dense_t, *dense_y,
      *preout, *ybuf, *dbuf;
  int *off;
  lin_sys *ls;
} arena;

static size_t align8(size_t x) { return (x + 7) & ~(size_t) 7; }

/* Lays the arena out, returning the bytes it takes (base NULL: count). */
static size_t layout(const zpx_blocks *b, int n, long entries, char *base, arena *A) {
  size_t at = 0;
  int ny = b->ny;
#define TAKE(field, type, count)                                \
  do {                                                          \
    if (base != NULL) A->field = (type *) (void *) (base + at); \
    at = align8(at + sizeof(type) * (size_t) (count));          \
  } while (0)
  /* What does not depend on the entry count first, so that its offsets
   * are the same however many entries there turn out to be: pre runs
   * before the count is known and leaves the initial amounts in y. */
  TAKE(y, double, ny + 1);
  TAKE(r, double, ny + 1);
  TAKE(s, double, ny + 1);
  TAKE(prev, double, ny + 1);
  TAKE(zero, double, ny + 1);
  TAKE(ybuf, double, ny + 1);
  TAKE(dbuf, double, ny + 1);
  TAKE(work, double, b->nwork + 1);
  TAKE(cap, double, b->ncap + 1);
  TAKE(pcur, double, b->npar + 1);
  TAKE(dwork, double, ZPX_DOPRI_WORK(ny) + 2 * ny + 1);
  TAKE(preout, double, 5 * ny + 1);
  TAKE(off, int, ny + 1);
  TAKE(ls, lin_sys, 1);
  TAKE(ent, zpx_entry, entries + 1);
  TAKE(active, double, entries + 1);
  TAKE(dense_t, double, entries + 1);
  TAKE(dense_y, double, (entries + 1) * ny);
#undef TAKE
  (void) n;
  return at;
}

size_t zpx_solve_arena_size(const zpx_blocks *b, int n, long entries) {
  return layout(b, n, entries, NULL, NULL);
}

/* ---- The closed form ---------------------------------------------------- */

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

/* ---- The matrix exponential --------------------------------------------- */

static void lin_probe(lin_sys *L, const zpx_blocks *b, double t, const double *par, double *work,
                      double *ybuf, double *dbuf) {
  int n = L->n;
  for (int i = 0; i < n; i++) ybuf[i] = 0;
  b->fn[ZPX_RHS](t, ybuf, L->c, par, work);
  for (int j = 0; j < n; j++) {
    ybuf[j] = 1;
    b->fn[ZPX_RHS](t, ybuf, dbuf, par, work);
    for (int i = 0; i < n; i++) L->a[i * n + j] = dbuf[i] - L->c[i];
    ybuf[j] = 0;
  }
}

static int lin_exp(lin_sys *L, double dt, const double *r, int with_c) {
  int n = L->n, m = n + 1;
  memset(L->M, 0, sizeof(double) * (size_t) (m * m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) L->M[i * m + j] = L->a[i * n + j] * dt;
    L->M[i * m + n] = ((with_c ? L->c[i] : 0) + (r != NULL ? r[i] : 0)) * dt;
  }
  return zpx_expm(m, L->M, L->X, L->mwork);
}

static int lin_step(lin_sys *L, double dt, double *y, const double *r) {
  int n = L->n, m = n + 1;
  if (!(dt > 0)) return 0;
  if (lin_exp(L, dt, r, 1) != 0) return 1;
  double out[ZPX_SOLVE_LINEAR_MAX];
  for (int i = 0; i < n; i++) {
    double s = L->X[i * m + n];
    for (int j = 0; j < n; j++) s += L->X[i * m + j] * y[j];
    out[i] = s;
  }
  memcpy(y, out, sizeof(double) * (size_t) n);
  return 0;
}

/* Solves (I - E) y = b, E the top-left block of L->X. */
static int lin_ss_solve(lin_sys *L, double *b, double *y) {
  int n = L->n, m = n + 1;
  double M[ZPX_SOLVE_LINEAR_MAX * ZPX_SOLVE_LINEAR_MAX];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) M[i * n + j] = (i == j ? 1.0 : 0.0) - L->X[i * m + j];
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

/* kind 0 bolus v every ii; 1 infusion v for dur every ii; 2 constant v. */
static int lin_ss(lin_sys *L, int kind, double ii, double dur, const double *v, double *y) {
  int n = L->n, m = n + 1;
  double b[ZPX_SOLVE_LINEAR_MAX] = {0};
  if (kind == 2) {
    /* A y + c + v = 0, as (I - E) y = c + v with E = I + A. */
    for (int i = 0; i < n; i++) {
      b[i] = L->c[i] + v[i];
      for (int j = 0; j < n; j++) L->X[i * m + j] = (i == j ? 1.0 : 0.0) + L->a[i * n + j];
    }
    return lin_ss_solve(L, b, y);
  }
  if (kind == 0) {
    if (lin_exp(L, ii, NULL, 1) != 0) return 1;
    for (int i = 0; i < n; i++) {
      double s = L->X[i * m + n];
      for (int j = 0; j < n; j++) s += L->X[i * m + j] * v[j];
      b[i] = s;
    }
    return lin_ss_solve(L, b, y);
  }
  double z[ZPX_SOLVE_LINEAR_MAX] = {0}, zero[ZPX_SOLVE_LINEAR_MAX] = {0};
  if (lin_step(L, dur, z, v) != 0) return 1;
  if (lin_step(L, ii - dur, z, zero) != 0) return 1;
  memcpy(b, z, sizeof(double) * (size_t) n);
  if (lin_exp(L, ii, NULL, 0) != 0) return 1;
  return lin_ss_solve(L, b, y);
}

/* ---- The general solvers ------------------------------------------------ */

typedef struct general {
  const zpx_blocks *b;
  const double *par;
  double *work;
  const double *rate;
  const int *off;
  int ny;
} general;

static void general_rhs(double t, const double *y, double *dydt, void *data) {
  general *g = data;
  g->b->fn[ZPX_RHS](t, y, dydt, g->par, g->work);
  for (int i = 0; i < g->ny; i++) {
    dydt[i] += g->rate[i];
    if (g->off[i]) dydt[i] = 0;
  }
}

static int lsoda_rhs(double t, double *y, double *dydt, void *data) {
  general_rhs(t, y, dydt, data);
  return 0;
}

static int finite_state(const double *y, int n) {
  for (int i = 0; i < n; i++)
    if (!isfinite(y[i])) return 0;
  return 1;
}

/* integrate()'s return when LSODA cannot allocate its workspace: an
 * error of the subject, not a solver failure to report as NA. */
#define INTEGRATE_NO_MEMORY 99
static const char *no_memory = "cannot allocate LSODA's workspace";

static int integrate(int method, general *g, double rtol, double atol, double t0, double t1,
                     double *y, double *dwork, zpx_dopri *D, double *reached) {
  *reached = t1;
  if (!(t1 > t0)) return 0;
  if (method == ZPX_M_DOPRI) {
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
  struct lsoda_opt_t opt;
  memset(&opt, 0, sizeof(opt));
  opt.ixpr = 0;
  opt.rtol = rtv;
  opt.atol = atv;
  opt.itask = 1;
  opt.mxstep = 100000;
  struct lsoda_context_t ctx;
  memset(&ctx, 0, sizeof(ctx));
  ctx.function = lsoda_rhs;
  ctx.data = g;
  ctx.neq = g->ny;
  ctx.state = 1;
  if (!lsoda_prepare(&ctx, &opt)) {
    /* The options are always valid here: a failed prepare is a failed
     * allocation (tools/vendor/liblsoda.sh, patch 12). */
    if (ctx.common != NULL) lsoda_free(&ctx);
    return INTEGRATE_NO_MEMORY;
  }
  double t = t0;
  lsoda(&ctx, y, &t, t1);
  /* LSODA can return a success short of tout, or carry on through a
   * blow-up; both are failures here. */
  int ok = ctx.state > 0 && t >= t1 && finite_state(y, g->ny);
  lsoda_free(&ctx);
  if (!ok) *reached = t;
  return ok ? 0 : 1;
}

/* ---- The walk ----------------------------------------------------------- */

static void reset_state(int ny, double *y, double *r, int *off, double *active, long nlink) {
  for (int j = 0; j < ny; j++) {
    y[j] = r[j] = 0;
    off[j] = 0;
  }
  for (long l = 0; l < nlink; l++) active[l] = 0;
}

void zpx_solve_subject(const zpx_blocks *b, const zpx_solve_opts *o, zpx_record *rec, int n,
                       const double *par, long par_ld, double *out, long ld, int *orow, long row0,
                       void *arena_mem, size_t arena_size, zpx_subject_result *res) {
  int ny = b->ny, npar = b->npar, ncap = b->ncap;
  int method = o->method;
  arena A;
  memset(res, 0, sizeof(*res));
  res->failed_at = NAN;
  res->failed_row = -1;
  res->ss_row = -1;
  res->error_row = -1;

  /* The pre block at each dose record (and the first, for the initial
   * amounts), with that record's parameters. It needs only the part of
   * the arena laid out before the entries. */
  int error, where;
  layout(b, n, 0, (char *) arena_mem, &A);
  for (int j = 0; j <= ny; j++) A.zero[j] = 0;
  for (int j = 0; j < ny; j++) A.y[j] = 0;
  for (int i = 0; i < n; i++) {
    rec[i].lag = 0;
    rec[i].f = 1;
    rec[i].mrate = 0;
    rec[i].mdur = 0;
    int dose = rec[i].evid == 1 || rec[i].evid == 4;
    if (i != 0 && !dose) continue;
    for (int p = 0; p < npar; p++) A.pcur[p] = par[i + par_ld * p];
    b->fn[ZPX_PRE](rec[i].time, A.zero, A.preout, A.pcur, A.work);
    if (i == 0)
      for (int j = 0; j < ny; j++) A.y[j] = A.preout[4 * ny + j];
    int c = rec[i].cmt - 1;
    if (dose && c >= 0 && c < ny) {
      rec[i].lag = A.preout[c];
      rec[i].f = A.preout[ny + c];
      rec[i].mrate = A.preout[2 * ny + c];
      rec[i].mdur = A.preout[3 * ny + c];
    }
  }
  long count = zpx_sched_count(rec, n, ny, &error, &where);
  if (count < 0) {
    res->error = ZPX_SOLVE_E_EVENTS;
    res->message = zpx_sched_error_message(error);
    res->error_row = where;
    return;
  }
  if (layout(b, n, count, NULL, NULL) > arena_size) {
    res->error = ZPX_SOLVE_E_SOLVER;
    res->message = "the subject's arena is too small (a bug in zupmxode)";
    return;
  }
  layout(b, n, count, (char *) arena_mem, &A);
  count = zpx_sched_build(rec, n, ny, A.ent, &error, &where);
  if (count < 0) {
    res->error = ZPX_SOLVE_E_EVENTS;
    res->message = zpx_sched_error_message(error);
    res->error_row = where;
    return;
  }

  long nlink = 0;
  for (long i = 0; i < count; i++) {
    if (A.ent[i].link + 1 > nlink) nlink = A.ent[i].link + 1;
    if (A.ent[i].kind == ZPX_E_OFF && method != ZPX_M_LSODA && method != ZPX_M_DOPRI) {
      res->error = ZPX_SOLVE_E_SOLVER;
      res->message = "a compartment switched off (EVID 2 with a negative CMT) needs lsoda or dopri5";
      res->error_row = A.ent[i].rec;
      return;
    }
  }
  if (method == ZPX_M_EXPM && ny > ZPX_SOLVE_LINEAR_MAX) {
    res->error = ZPX_SOLVE_E_SOLVER;
    res->message = "too many compartments for the matrix exponential";
    return;
  }
  for (int j = 0; j < ny; j++) {
    A.r[j] = 0;
    A.off[j] = 0;
  }
  for (long l = 0; l < nlink; l++) A.active[l] = 0;

  zpx_lin L;
  A.ls->n = ny;
  general g = {b, A.pcur, A.work, A.r, A.off, ny};
  zpx_dopri D = {ny, o->rtol, o->atol, 1000000, 0, 0, 0, 0};
  int cur = -1, failed = 0;
  double t = count > 0 ? A.ent[0].time : 0;
  long dense_from = -1, dense_to = -1;

#define SET_RECORD(i)                                                                     \
  do {                                                                                    \
    cur = (i);                                                                            \
    for (int p = 0; p < npar; p++) A.pcur[p] = par[cur + par_ld * p];                   \
    b->fn[ZPX_IVL](t, A.zero, A.work, A.pcur, A.work);                                    \
    if (method == ZPX_M_CLOSED) {                                                         \
      b->fn[ZPX_POST](t, A.zero, A.cap, A.pcur, A.work);                                  \
      double ka, k10, k12, k21, k13, k31;                                                 \
      int nc = rate_constants(o->lc, A.cap, &ka, &k10, &k12, &k21, &k13, &k31);           \
      if (zpx_lin_setup(&L, nc, o->lc[1], ka, k10, k12, k21, k13, k31) != ZPX_LIN_OK) {   \
        res->error = ZPX_SOLVE_E_SOLVER;                                                  \
        res->message = "a rate constant of the closed form is not a finite positive number"; \
        res->error_row = cur;                                                             \
        return;                                                                           \
      }                                                                                   \
    } else if (method == ZPX_M_EXPM) {                                                    \
      lin_probe(A.ls, b, t, A.pcur, A.work, A.ybuf, A.dbuf);                              \
    } else {                                                                              \
      D.h = 0;                                                                            \
    }                                                                                     \
  } while (0)

  SET_RECORD(0);
  long k = row0;
  for (long i = 0; i < count; i++) {
    const zpx_entry *e = &A.ent[i];
    int dense_hit = method == ZPX_M_DOPRI && i > dense_from && i <= dense_to;
    if (!failed && e->time > t && !dense_hit) {
      int st = 0;
      double reached = e->time;
      if (method == ZPX_M_CLOSED) {
        zpx_lin_step(&L, e->time - t, A.y, A.r);
      } else if (method == ZPX_M_EXPM) {
        if (lin_step(A.ls, e->time - t, A.y, A.r) != 0) {
          res->error = ZPX_SOLVE_E_SOLVER;
          res->message = "the linear system's matrix is not finite";
          res->error_row = e->rec;
          return;
        }
      } else {
        long last = i;
        if (method == ZPX_M_DOPRI && e->kind == ZPX_E_OBS) {
          /* A run of observations, nothing between them that changes the
           * system: one integration, read by dense output. */
          while (last + 1 < count && A.ent[last + 1].kind == ZPX_E_OBS &&
                 A.ent[last + 1].rec == A.ent[last].rec + 1 && A.ent[last + 1].phase == 2)
            last++;
          long cut = i;
          while (cut + 1 <= last) {
            int same = 1;
            for (int p = 0; p < npar && same; p++)
              same = par[A.ent[cut + 1].rec + par_ld * p] == par[cur + par_ld * p];
            if (!same) break;
            cut++;
          }
          last = cut;
        }
        if (last > i) {
          int nd = (int) (last - i);
          for (int q = 0; q < nd; q++) A.dense_t[q] = A.ent[i + q].time;
          D.n = ny;
          D.rtol = o->rtol;
          D.atol = o->atol;
          st = zpx_dopri_solve(&D, general_rhs, &g, t, A.ent[last].time, A.y, A.dense_t, nd,
                               A.dense_y, A.dwork);
          if (st == 0 && !finite_state(A.y, ny)) st = ZPX_DOPRI_E_NONFINITE;
          reached = st != 0 ? D.t : A.ent[last].time;
          dense_from = i;
          dense_to = st != 0 ? i + D.nout_done : last;
        } else {
          st = integrate(method, &g, o->rtol, o->atol, t, e->time, A.y, A.dwork, &D, &reached);
        }
      }
      if (st == INTEGRATE_NO_MEMORY) {
        res->error = ZPX_SOLVE_E_SOLVER;
        res->message = no_memory;
        res->error_row = e->rec;
        return;
      }
      if (st != 0) {
        failed = 1;
        res->failed_at = reached;
        res->failed_row = e->rec;
      }
    }
    t = e->time;
    if (e->phase == 2 && e->rec != cur) SET_RECORD(e->rec);
    int c = e->cmt - 1;
    if (e->kind == ZPX_E_OBS && method == ZPX_M_DOPRI && i >= dense_from && i < dense_to) {
      double *yd = A.dense_y + (i - dense_from) * (long) ny;
      b->fn[ZPX_POST](t, yd, A.cap, A.pcur, A.work);
      out[k] = t;
      for (int j = 0; j < ny; j++) out[k + ld * (1 + j)] = yd[j];
      for (int j = 0; j < ncap; j++) out[k + ld * (1 + ny + j)] = A.cap[j];
      orow[k] = e->rec;
      k++;
      continue;
    }
    switch (e->kind) {
    case ZPX_E_OBS:
      out[k] = t;
      if (failed) {
        for (int j = 0; j < ny + ncap; j++) out[k + ld * (1 + j)] = NAN;
      } else {
        b->fn[ZPX_POST](t, A.y, A.cap, A.pcur, A.work);
        for (int j = 0; j < ny; j++) out[k + ld * (1 + j)] = A.y[j];
        for (int j = 0; j < ncap; j++) out[k + ld * (1 + ny + j)] = A.cap[j];
      }
      orow[k] = e->rec;
      k++;
      break;
    case ZPX_E_BOLUS:
      A.y[c] += e->amount;
      A.off[c] = 0; /* a dose turns its compartment on (PREDPP) */
      D.h = 0;
      break;
    case ZPX_E_INF_START:
      A.active[e->link] = e->rate;
      A.r[c] += e->rate;
      A.off[c] = 0;
      D.h = 0;
      break;
    case ZPX_E_INF_END:
      /* An infusion a reset stopped has nothing left to end. */
      if (A.active[e->link] != 0) {
        A.r[c] -= A.active[e->link];
        A.active[e->link] = 0;
        if (fabs(A.r[c]) < 1e-12 * e->rate) A.r[c] = 0;
      }
      D.h = 0;
      break;
    case ZPX_E_RESET:
      reset_state(ny, A.y, A.r, A.off, A.active, nlink);
      D.h = 0;
      break;
    case ZPX_E_OFF:
      A.y[c] = 0;
      A.off[c] = 1;
      D.h = 0;
      break;
    case ZPX_E_ON: A.off[c] = 0; break;
    case ZPX_E_SS: {
      if (failed) break;
      for (int j = 0; j < ny; j++) A.s[j] = A.dbuf[j] = 0;
      double *v = A.dbuf;
      if (method == ZPX_M_CLOSED || method == ZPX_M_EXPM) {
        if (e->ss == 1) reset_state(ny, A.y, A.r, A.off, A.active, nlink);
        int kind = e->rate > 0 ? (e->ii == 0 ? 2 : 1) : 0;
        v[c] = kind == 0 ? e->amount : e->rate;
        if (method == ZPX_M_CLOSED) {
          if (kind == 2) zpx_lin_ss_constant(&L, v, A.s);
          else if (kind == 1) zpx_lin_ss_infusion(&L, e->ii, e->dur, v, A.s);
          else zpx_lin_ss_bolus(&L, e->ii, v, A.s);
        } else {
          /* SS = 2 adds the dose's own steady state to what is there,
           * which already carries the constant input's share: no c. */
          double saved[ZPX_SOLVE_LINEAR_MAX];
          memcpy(saved, A.ls->c, sizeof(double) * (size_t) ny);
          if (e->ss == 2)
            for (int j = 0; j < ny; j++) A.ls->c[j] = 0;
          int bad = lin_ss(A.ls, kind, e->ii, e->dur, v, A.s);
          memcpy(A.ls->c, saved, sizeof(double) * (size_t) ny);
          if (bad) {
            res->error = ZPX_SOLVE_E_SOLVER;
            res->message = "the linear system has no steady state";
            res->error_row = e->rec;
            return;
          }
        }
        for (int j = 0; j < ny; j++) A.y[j] += A.s[j];
        break;
      }
      /* Repeated dosing from empty, until the state before a dose moves
       * by less than ss_tol; at most ss_cap intervals. */
      /* The infusions running now are put back after: ybuf holds them
       * (it is the matrix exponential's alone, not this path's). */
      memcpy(A.ybuf, A.r, sizeof(double) * (size_t) ny);
      int converged = 0;
      double ss_reached;
      int ist = 0;
      for (long it = 0; it < o->ss_cap && !failed; it++) {
        memcpy(A.prev, A.s, sizeof(double) * (size_t) ny);
        for (int j = 0; j < ny; j++) A.r[j] = 0;
        if (e->rate > 0 && e->ii == 0) {
          A.r[c] = e->rate;
          if ((ist = integrate(method, &g, o->rtol, o->atol, 0, 24, A.s, A.dwork, &D, &ss_reached)) != 0)
            failed = ist;
        } else if (e->rate > 0) {
          A.r[c] = e->rate;
          if ((ist = integrate(method, &g, o->rtol, o->atol, 0, e->dur, A.s, A.dwork, &D, &ss_reached)) != 0)
            failed = ist;
          A.r[c] = 0;
          if (!failed &&
              (ist = integrate(method, &g, o->rtol, o->atol, e->dur, e->ii, A.s, A.dwork, &D, &ss_reached)) != 0)
            failed = ist;
        } else {
          A.s[c] += e->amount;
          if ((ist = integrate(method, &g, o->rtol, o->atol, 0, e->ii, A.s, A.dwork, &D, &ss_reached)) != 0)
            failed = ist;
        }
        double change = 0;
        for (int j = 0; j < ny; j++) {
          double d = fabs(A.s[j] - A.prev[j]) /
                     (o->atol + o->ss_tol * fmax(fabs(A.s[j]), fabs(A.prev[j])));
          if (d > change) change = d;
        }
        if (it > 0 && change <= 1) {
          converged = 1;
          break;
        }
      }
      memcpy(A.r, A.ybuf, sizeof(double) * (size_t) ny);
      if (failed == INTEGRATE_NO_MEMORY) {
        res->error = ZPX_SOLVE_E_SOLVER;
        res->message = no_memory;
        res->error_row = e->rec;
        return;
      }
      if (failed) {
        res->failed_at = t;
        res->failed_row = e->rec;
        break;
      }
      if (!converged) {
        res->ss_unconverged++;
        res->ss_row = e->rec;
      }
      if (e->ss == 1) reset_state(ny, A.y, A.r, A.off, A.active, nlink);
      for (int j = 0; j < ny; j++) A.y[j] += A.s[j];
      D.h = 0;
      break;
    }
    default: break;
    }
  }
#undef SET_RECORD
}
