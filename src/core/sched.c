/* The event scheduler (src/core/sched.h). Builds without R. */
#include "sched.h"

#include <math.h>
#include <stdlib.h>

const char *zpx_sched_error_message(int error) {
  switch (error) {
  case ZPX_SCHED_OK: return "ok";
  case ZPX_SCHED_E_TIME: return "TIME is missing, not finite, or decreases";
  case ZPX_SCHED_E_EVID: return "EVID must be 0, 1, 2, 3 or 4";
  case ZPX_SCHED_E_AMT:
    return "AMT must be a finite amount >= 0 on a dose (EVID 1 or 4) and 0 on any other record";
  case ZPX_SCHED_E_CMT: return "CMT is not a compartment of the model";
  case ZPX_SCHED_E_RATE: return "RATE must be >= 0, -1 or -2 on a dose and 0 on any other record";
  case ZPX_SCHED_E_MRATE: return "RATE is -1 but the model's rate() for the compartment is not > 0";
  case ZPX_SCHED_E_MDUR: return "RATE is -2 but the model's dur() for the compartment is not > 0";
  case ZPX_SCHED_E_DUR: return "DUR must be >= 0, and 0 when RATE is given";
  case ZPX_SCHED_E_ADDL: return "ADDL must be >= 0, and needs II > 0";
  case ZPX_SCHED_E_II: return "II must be a finite interval >= 0";
  case ZPX_SCHED_E_SS:
    return "SS must be 0, 1 or 2, only on a dose, and needs II > 0 unless the dose is an infusion";
  case ZPX_SCHED_E_LAGF:
    return "the model's lag and bioavailability for the compartment must be finite and >= 0";
  case ZPX_SCHED_E_TOO_MANY: return "the records make more timeline entries than the limit";
  case ZPX_SCHED_E_NO_MEMORY: return "not enough memory for the timeline";
  case ZPX_SCHED_E_INTERNAL: return "the scheduler's two passes disagree (a bug in zupmxode)";
  }
  return "unknown scheduler error";
}

static int is_dose(int evid) { return evid == 1 || evid == 4; }

/* Whether a dose record is an infusion, and its rate and duration once
 * bioavailability is applied: a given RATE keeps its rate (F changes the
 * duration, as in PREDPP); a given DUR, or the model's dur(), keeps its
 * duration (F changes the rate). */
static int infusion(const zpx_record *r, double *rate, double *dur) {
  double a = r->amt * r->f;
  if (r->rate > 0) {
    *rate = r->rate;
    *dur = a / r->rate;
    return 1;
  }
  if (r->rate == -1) {
    *rate = r->mrate;
    *dur = a / r->mrate;
    return 1;
  }
  if (r->rate == -2) {
    *dur = r->mdur;
    *rate = a / r->mdur;
    return 1;
  }
  if (r->dur > 0) {
    *dur = r->dur;
    *rate = a / r->dur;
    return 1;
  }
  return 0;
}

/* Checks one record; the error, or 0. */
static int check(const zpx_record *r, double prev_time, int ncmt) {
  if (!isfinite(r->time) || r->time < prev_time) return ZPX_SCHED_E_TIME;
  if (r->evid < 0 || r->evid > 4) return ZPX_SCHED_E_EVID;
  int dose = is_dose(r->evid);
  if (!isfinite(r->amt) || r->amt < 0 || (!dose && r->amt != 0)) return ZPX_SCHED_E_AMT;
  if (dose && (r->cmt < 1 || r->cmt > ncmt)) return ZPX_SCHED_E_CMT;
  if (r->evid == 0 && (r->cmt < 0 || r->cmt > ncmt)) return ZPX_SCHED_E_CMT;
  if (r->evid == 2 && (r->cmt < -ncmt || r->cmt > ncmt)) return ZPX_SCHED_E_CMT;
  if (!(isfinite(r->rate) && (r->rate >= 0 || r->rate == -1 || r->rate == -2)) ||
      (!dose && r->rate != 0))
    return ZPX_SCHED_E_RATE;
  if (!isfinite(r->dur) || r->dur < 0 || (r->dur > 0 && r->rate != 0) || (!dose && r->dur != 0))
    return ZPX_SCHED_E_DUR;
  if (!isfinite(r->ii) || r->ii < 0) return ZPX_SCHED_E_II;
  if (r->addl < 0 || (r->addl > 0 && !(r->ii > 0))) return ZPX_SCHED_E_ADDL;
  if (dose) {
    if (!isfinite(r->lag) || r->lag < 0 || !isfinite(r->f) || r->f < 0) return ZPX_SCHED_E_LAGF;
    if (r->rate == -1 && !(isfinite(r->mrate) && r->mrate > 0)) return ZPX_SCHED_E_MRATE;
    if (r->rate == -2 && !(isfinite(r->mdur) && r->mdur > 0)) return ZPX_SCHED_E_MDUR;
  }
  if (r->ss < 0 || r->ss > 2 || (r->ss > 0 && !dose)) return ZPX_SCHED_E_SS;
  if (r->ss > 0) {
    double rate, dur;
    int inf = infusion(r, &rate, &dur);
    if (!(r->ii > 0) && !inf) return ZPX_SCHED_E_SS;
  }
  return ZPX_SCHED_OK;
}

/* Entries one checked record makes. */
static long entries(const zpx_record *r) {
  if (!is_dose(r->evid)) return 1;
  double rate, dur;
  int inf = infusion(r, &rate, &dur);
  /* A steady-state infusion with II 0 runs on: no end. */
  int endless = inf && r->ss > 0 && r->ii == 0;
  long per = inf ? (endless ? 1 : 2) : 1;
  long n = (r->evid == 4) + (r->ss > 0) + per;
  return n + (long) r->addl * per;
}

long zpx_sched_count(const zpx_record *rec, int n, int ncmt, int *error, int *where) {
  long total = 0;
  *error = ZPX_SCHED_OK;
  *where = -1;
  for (int i = 0; i < n; i++) {
    int e = check(&rec[i], i > 0 ? rec[i - 1].time : -INFINITY, ncmt);
    if (e == ZPX_SCHED_OK && rec[i].addl > ZPX_SCHED_MAX_ENTRIES) e = ZPX_SCHED_E_TOO_MANY;
    if (e != ZPX_SCHED_OK) {
      *error = e;
      *where = i;
      return -1;
    }
    total += entries(&rec[i]);
    if (total > ZPX_SCHED_MAX_ENTRIES) {
      *error = ZPX_SCHED_E_TOO_MANY;
      *where = i;
      return -1;
    }
  }
  return total;
}

static int cmp_entry(const void *a, const void *b) {
  const zpx_entry *x = a, *y = b;
  if (x->time != y->time) return x->time < y->time ? -1 : 1;
  if (x->phase != y->phase) return x->phase < y->phase ? -1 : 1;
  return x->seq < y->seq ? -1 : x->seq > y->seq;
}

long zpx_sched_build(const zpx_record *rec, int n, int ncmt, zpx_entry *out, int *error,
                     int *where) {
  long count = zpx_sched_count(rec, n, ncmt, error, where);
  if (count < 0) return -1;
  long k = 0, gen = 0, links = 0;
  for (int i = 0; i < n; i++) {
    const zpx_record *r = &rec[i];
    /* Record entries keep data order: four slots per record (reset, steady
     * state, dose, other). */
    long base = (long) i * 4;
    zpx_entry e = {0};
    e.link = -1;
    e.rec = i;
    e.time = r->time;
    e.phase = 2;
    if (!is_dose(r->evid)) {
      e.kind = r->evid == 0 ? ZPX_E_OBS
               : r->evid == 3 ? ZPX_E_RESET
               : r->cmt < 0 ? ZPX_E_OFF
               : r->cmt > 0 ? ZPX_E_ON
                            : ZPX_E_MARK;
      /* A reset and a mark are about no one compartment. */
      e.cmt = (e.kind == ZPX_E_RESET || e.kind == ZPX_E_MARK) ? 0 : (r->cmt < 0 ? -r->cmt : r->cmt);
      e.seq = base + 3;
      out[k++] = e;
      continue;
    }
    if (r->evid == 4) {
      e.kind = ZPX_E_RESET;
      e.seq = base;
      out[k++] = e;
    }
    double rate = 0, dur = 0;
    int inf = infusion(r, &rate, &dur);
    int endless = inf && r->ss > 0 && r->ii == 0;
    double amount = r->amt * r->f;
    for (long d = 0; d <= r->addl; d++) {
      zpx_entry s = {0};
      s.link = -1;
      s.rec = i;
      s.cmt = r->cmt;
      /* k * II from the record's time, not summed, so the doses do not
       * drift. */
      s.time = r->time + (double) d * r->ii + r->lag;
      if (d == 0 && r->lag == 0) {
        s.phase = 2;
        s.seq = base + 2;
      } else {
        s.phase = 1;
        s.seq = gen++;
      }
      if (d == 0 && r->ss > 0) {
        zpx_entry ss = s;
        ss.kind = ZPX_E_SS;
        ss.ss = r->ss;
        ss.amount = amount;
        ss.rate = inf ? rate : 0;
        ss.dur = inf ? dur : 0;
        ss.ii = r->ii;
        /* Before the dose it precedes, in either phase. */
        ss.seq = s.phase == 2 ? base + 1 : gen++;
        if (s.phase == 1) s.seq = gen++;
        out[k++] = ss;
      }
      if (inf) {
        s.kind = ZPX_E_INF_START;
        s.link = links++;
        s.rate = rate;
        s.dur = dur;
        out[k++] = s;
        if (!endless) {
          zpx_entry end = s;
          end.kind = ZPX_E_INF_END;
          end.time = s.time + dur;
          end.phase = 0;
          end.seq = gen++;
          out[k++] = end;
        }
      } else {
        s.kind = ZPX_E_BOLUS;
        s.amount = amount;
        out[k++] = s;
      }
    }
  }
  if (k != count) {
    /* The two passes disagree: a bug here, never the caller's input. */
    *error = ZPX_SCHED_E_INTERNAL;
    *where = -1;
    return -1;
  }
  qsort(out, (size_t) k, sizeof(zpx_entry), cmp_entry);
  return k;
}
