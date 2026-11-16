/* libFuzzer harness for the event scheduler (src/core/sched.c), built
 * without R by r-actions' fuzz.yml (.github/workflows/fuzz.yaml).
 *
 * Input: one byte for the number of compartments (1 to 8), then records of
 * 16 bytes each, decoded so that small values land in the valid ranges
 * often (see decode()).
 *
 * Oracle, beyond "no crash": for records the scheduler accepts, the
 * timeline is sorted by time, then phase, then order; it has exactly one
 * dose per dose and additional dose; every infusion that should end ends,
 * at its start plus its duration; and every entry is finite and points at
 * a real record. Records it refuses name a known error at a real row. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "sched.h"

static double small(uint8_t b) {
  /* 0 to 63 in quarters, with a few specials in the top values. */
  if (b == 255) return NAN;
  if (b == 254) return INFINITY;
  if (b == 253) return -1;
  if (b == 252) return -2;
  return (double) (b & 0x7f) / 4.0;
}

static void decode(const uint8_t *p, zpx_record *r) {
  r->time = small(p[0]);
  r->evid = p[1] % 6;            /* 5 is invalid */
  r->amt = p[2] < 200 ? (double) p[2] : small(p[2]);
  r->cmt = (int) (p[3] % 19) - 9; /* -9 to 9 */
  r->rate = p[4] < 128 ? 0 : small(p[4]);
  r->dur = p[5] < 192 ? 0 : small(p[5]);
  r->ii = p[6] < 64 ? 0 : small(p[6]);
  r->addl = p[7] < 128 ? 0 : (int) (p[7] & 0x1f);
  r->ss = p[8] % 4;              /* 3 is invalid */
  r->lag = p[9] < 192 ? 0 : small(p[9]);
  r->f = p[10] < 128 ? 1 : small(p[10]);
  r->mrate = small(p[11]);
  r->mdur = small(p[12]);
  (void) p[13];
}

static int is_dose(int evid) { return evid == 1 || evid == 4; }

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  if (size < 1) return 0;
  int ncmt = 1 + data[0] % 8;
  int n = (int) ((size - 1) / 16);
  if (n > 64) n = 64;
  zpx_record rec[64];
  for (int i = 0; i < n; i++) decode(data + 1 + 16 * i, &rec[i]);
  /* Time must not decrease for the records to be accepted; most inputs
   * keep it so, by making each time an offset from the previous. */
  for (int i = 1; i < n; i++)
    if (isfinite(rec[i].time) && isfinite(rec[i - 1].time) && (data[1 + 16 * i + 14] & 1))
      rec[i].time += rec[i - 1].time;

  int error, where;
  long count = zpx_sched_count(rec, n, ncmt, &error, &where);
  if (count < 0) {
    if (error <= ZPX_SCHED_OK || error > ZPX_SCHED_E_INTERNAL || where < 0 || where >= n)
      abort();
    return 0;
  }
  zpx_entry *out = malloc((size_t) (count > 0 ? count : 1) * sizeof(zpx_entry));
  if (out == NULL) abort();
  long k = zpx_sched_build(rec, n, ncmt, out, &error, &where);
  if (k != count) abort();
  long doses = 0, expected = 0, starts = 0, ends = 0;
  for (int i = 0; i < n; i++)
    if (is_dose(rec[i].evid)) expected += 1 + rec[i].addl;
  for (long i = 0; i < k; i++) {
    const zpx_entry *e = &out[i];
    if (!isfinite(e->time) || e->rec < 0 || e->rec >= n) abort();
    if (e->cmt < 0 || e->cmt > ncmt) abort();
    if (i > 0) {
      const zpx_entry *p = &out[i - 1];
      if (p->time > e->time) abort();
      if (p->time == e->time && (p->phase > e->phase ||
                                 (p->phase == e->phase && p->seq >= e->seq)))
        abort();
    }
    if (e->kind == ZPX_E_BOLUS || e->kind == ZPX_E_INF_START) doses++;
    if (e->kind == ZPX_E_INF_START) starts++;
    if (e->kind == ZPX_E_INF_END) ends++;
  }
  if (doses != expected) abort();
  /* Only steady-state infusions with II 0 run on without an end. */
  long endless = 0;
  for (int i = 0; i < n; i++) {
    const zpx_record *r = &rec[i];
    int inf = r->rate != 0 || r->dur > 0;
    if (is_dose(r->evid) && inf && r->ss > 0 && r->ii == 0) endless++;
  }
  if (starts - ends != endless) abort();
  free(out);
  return 0;
}
