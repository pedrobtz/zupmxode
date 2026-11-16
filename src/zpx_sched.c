/* The scheduler's .Call entry (src/core/sched.c). */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <stdlib.h>

#include "core/sched.h"

static double num(SEXP x, R_xlen_t i) { return REAL(x)[i]; }
static int ival(SEXP x, R_xlen_t i) { return INTEGER(x)[i]; }

SEXP zpx_schedule_r(SEXP time, SEXP evid, SEXP amt, SEXP cmt, SEXP rate, SEXP dur, SEXP ii,
                    SEXP addl, SEXP ss, SEXP lag, SEXP f, SEXP mrate, SEXP mdur, SEXP s_ncmt) {
  R_xlen_t n = XLENGTH(time);
  SEXP dbl[] = {time, amt, rate, dur, ii, lag, f, mrate, mdur};
  SEXP ints[] = {evid, cmt, addl, ss};
  for (int k = 0; k < 9; k++)
    if (TYPEOF(dbl[k]) != REALSXP || XLENGTH(dbl[k]) != n) Rf_error("bad double column %d", k);
  for (int k = 0; k < 4; k++)
    if (TYPEOF(ints[k]) != INTSXP || XLENGTH(ints[k]) != n) Rf_error("bad integer column %d", k);
  if (n > 100000000) Rf_error("too many records");
  zpx_record *rec = (zpx_record *) R_alloc((size_t) n + 1, sizeof(zpx_record));
  for (R_xlen_t i = 0; i < n; i++) {
    rec[i].time = num(time, i);
    rec[i].amt = num(amt, i);
    rec[i].rate = num(rate, i);
    rec[i].dur = num(dur, i);
    rec[i].ii = num(ii, i);
    rec[i].lag = num(lag, i);
    rec[i].f = num(f, i);
    rec[i].mrate = num(mrate, i);
    rec[i].mdur = num(mdur, i);
    rec[i].evid = ival(evid, i);
    rec[i].cmt = ival(cmt, i);
    rec[i].addl = ival(addl, i);
    rec[i].ss = ival(ss, i);
  }
  int ncmt = Rf_asInteger(s_ncmt), error, where;
  long count = zpx_sched_count(rec, (int) n, ncmt, &error, &where);
  zpx_entry *out = NULL;
  if (count >= 0) {
    out = (zpx_entry *) R_alloc((size_t) count + 1, sizeof(zpx_entry));
    count = zpx_sched_build(rec, (int) n, ncmt, out, &error, &where);
  }
  if (count < 0) {
    const char *names[] = {"error", "row", ""};
    SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, Rf_mkString(zpx_sched_error_message(error)));
    SET_VECTOR_ELT(res, 1, Rf_ScalarInteger(where));
    UNPROTECT(1);
    return res;
  }
  const char *names[] = {"time", "kind", "cmt", "amount", "rate", "dur", "ii", "ss", "rec", ""};
  SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
  SEXPTYPE types[] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP};
  for (int k = 0; k < 9; k++) SET_VECTOR_ELT(res, k, Rf_allocVector(types[k], (R_xlen_t) count));
  for (long i = 0; i < count; i++) {
    REAL(VECTOR_ELT(res, 0))[i] = out[i].time;
    INTEGER(VECTOR_ELT(res, 1))[i] = out[i].kind;
    INTEGER(VECTOR_ELT(res, 2))[i] = out[i].cmt;
    REAL(VECTOR_ELT(res, 3))[i] = out[i].amount;
    REAL(VECTOR_ELT(res, 4))[i] = out[i].rate;
    REAL(VECTOR_ELT(res, 5))[i] = out[i].dur;
    REAL(VECTOR_ELT(res, 6))[i] = out[i].ii;
    INTEGER(VECTOR_ELT(res, 7))[i] = out[i].ss;
    INTEGER(VECTOR_ELT(res, 8))[i] = out[i].rec;
  }
  UNPROTECT(1);
  return res;
}
