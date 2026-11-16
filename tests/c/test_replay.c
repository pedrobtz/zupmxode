/* Corpus cases replayed through src/core/ alone: the scheduler, every
 * solver and the reference evaluator for the model's blocks, compared with
 * pmx_ode_solve()'s output (tests/c/replay/, written by
 * tools/make-c-replay.R; the R suite ties that output to rxode2 and
 * mrgsolve). The mutation job runs this binary, so a mutant that changes a
 * number in the core is caught here.
 *
 * Not bit for bit: the cases are written on one platform and replayed on
 * every other, and the C library's exp(), log() and pow() differ between
 * them in the last bit, as does the compiler's use of fused multiply-add
 * (platforms.yaml, Windows and macOS arm64, Stage 12). The closed form and
 * the matrix exponential carry that ulp through, so they are held to
 * 1e-12 relative. LSODA and Dormand-Prince choose their steps from the
 * values, so an ulp changes the steps and the result moves within the
 * solver's error control (1e-9 seen on arm64); they are held to 1e-8, a
 * hundredth of their rtol. On one platform the replay is exact. Arguments: the case files; with
 * none, every file in ZPX_REPLAY_DIR (mull-runner passes none). */
#define _POSIX_C_SOURCE 200809L
#include <dirent.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ir.h"
#include "sched.h"
#include "solve.h"

static zpx_ir g_ir;
static double *g_regs;

#define BLOCK(k)                                                                         \
  static void block##k(double t, const double *y, double *out, const double *par,       \
                       double *work) {                                                  \
    zpx_ir_eval(&g_ir, k, t, y, out, par, work, g_regs);                                \
  }
BLOCK(0)
BLOCK(1)
BLOCK(2)
BLOCK(3)
BLOCK(4)

static double hex_double(const char *h) {
  uint64_t u = (uint64_t) strtoull(h, NULL, 16);
  double d;
  memcpy(&d, &u, sizeof d);
  return d;
}

static int read_doubles(FILE *f, double *x, long n) {
  char h[32];
  for (long i = 0; i < n; i++) {
    if (fscanf(f, "%31s", h) != 1) return 0;
    x[i] = hex_double(h);
  }
  return 1;
}

static int expect(FILE *f, const char *word) {
  char w[64];
  return fscanf(f, "%63s", w) == 1 && strcmp(w, word) == 0;
}

/* 0 when the case replays within the tolerance; 1 when it differs; 2 when the
 * file cannot be read. */
static int replay(const char *path) {
  FILE *f = fopen(path, "r");
  if (f == NULL) return 2;
  char name[256];
  long ncode, nconst, npar_total, nobs, ncol;
  int ncap, method, n, lc[10], ss_cap;
  double tol[3];
  int ok = expect(f, "case") && fscanf(f, "%255s", name) == 1 && expect(f, "code") &&
           fscanf(f, "%ld", &ncode) == 1;
  int32_t *code = ok ? malloc((size_t) ncode * sizeof(int32_t) + 1) : NULL;
  for (long i = 0; ok && i < ncode; i++) ok = fscanf(f, "%d", &code[i]) == 1;
  ok = ok && expect(f, "consts") && fscanf(f, "%ld", &nconst) == 1;
  double *consts = ok ? malloc((size_t) nconst * sizeof(double) + 1) : NULL;
  ok = ok && read_doubles(f, consts, nconst) && expect(f, "ncap") &&
       fscanf(f, "%d", &ncap) == 1 && expect(f, "method") && fscanf(f, "%d", &method) == 1 &&
       expect(f, "opts") && read_doubles(f, tol, 3) && fscanf(f, "%d", &ss_cap) == 1 &&
       expect(f, "lc");
  for (int i = 0; ok && i < 10; i++) ok = fscanf(f, "%d", &lc[i]) == 1;
  ok = ok && expect(f, "records") && fscanf(f, "%d", &n) == 1;
  zpx_record *rec = ok ? calloc((size_t) n + 1, sizeof(zpx_record)) : NULL;
  for (int i = 0; ok && i < n; i++) {
    double v[5];
    ok = read_doubles(f, v, 5) && fscanf(f, "%d %d %d %d", &rec[i].evid, &rec[i].cmt,
                                         &rec[i].addl, &rec[i].ss) == 4;
    rec[i].time = v[0];
    rec[i].amt = v[1];
    rec[i].rate = v[2];
    rec[i].dur = v[3];
    rec[i].ii = v[4];
    rec[i].f = 1;
    rec[i].mrate = 1;
    rec[i].mdur = 1;
  }
  ok = ok && expect(f, "par") && fscanf(f, "%ld", &npar_total) == 1;
  double *par = ok ? malloc((size_t) npar_total * sizeof(double) + 1) : NULL;
  ok = ok && read_doubles(f, par, npar_total) && expect(f, "out") &&
       fscanf(f, "%ld %ld", &nobs, &ncol) == 2;
  double *want = ok ? malloc((size_t) (nobs * ncol) * sizeof(double) + 1) : NULL;
  ok = ok && read_doubles(f, want, nobs * ncol);
  fclose(f);
  int block, insn;
  ok = ok && zpx_ir_verify(code, (size_t) ncode, consts, (size_t) nconst, &g_ir, &block,
                           &insn) == ZPX_IR_OK;
  if (!ok) {
    fprintf(stderr, "FAIL: %s: cannot read the case\n", path);
    return 2;
  }
  int maxreg = 1;
  for (int b = 0; b < ZPX_IR_NBLOCKS; b++)
    if (g_ir.block[b].nreg > maxreg) maxreg = g_ir.block[b].nreg;
  g_regs = malloc((size_t) maxreg * sizeof(double));

  zpx_blocks bl = {{block0, block1, block2, block3, block4},
                   g_ir.ny, g_ir.npar, g_ir.nwork, ncap};
  zpx_solve_opts o;
  memset(&o, 0, sizeof o);
  o.method = method;
  o.rtol = tol[0];
  o.atol = tol[1];
  o.ss_tol = tol[2];
  o.ss_cap = ss_cap;
  memcpy(o.lc, lc, sizeof o.lc);
  int error, where;
  long count = zpx_sched_count(rec, n, g_ir.ny, &error, &where);
  size_t size = zpx_solve_arena_size(&bl, n, count < 0 ? 0 : count);
  void *arena = malloc(size + 8);
  double *got = malloc((size_t) (nobs * (ncol + 1)) * sizeof(double) + 1);
  int *orow = malloc((size_t) nobs * sizeof(int) + 1);
  for (long i = 0; i < nobs * (ncol + 1); i++) got[i] = NAN;
  zpx_subject_result res;
  zpx_solve_subject(&bl, &o, rec, n, par, n, got, nobs, orow, 0, arena, size, &res);
  /* out holds time, the compartments and the outputs; the file holds them
   * without the id. Both NaN, of any payload, is the same; infinities
   * must match exactly. */
  int differ = res.error != 0;
  for (long k = 0; !differ && k < nobs * ncol; k++) {
    double a = got[k], b = want[k];
    double rel = method == ZPX_M_LSODA || method == ZPX_M_DOPRI ? 1e-8 : 1e-12;
    int same = (isnan(a) && isnan(b)) || a == b ||
               (isfinite(a) && isfinite(b) && fabs(a - b) <= rel * fmax(fabs(a), fabs(b)));
    if (!same) {
      fprintf(stderr, "FAIL: %s: row %ld column %ld: %.17g, want %.17g\n", name,
              k % nobs + 1, k / nobs + 1, a, b);
      differ = 1;
    }
  }
  if (res.error != 0) fprintf(stderr, "FAIL: %s: %s\n", name, res.message);
  free(code);
  free(consts);
  free(rec);
  free(par);
  free(want);
  free(got);
  free(orow);
  free(arena);
  free(g_regs);
  return differ;
}

#ifndef ZPX_REPLAY_DIR
#define ZPX_REPLAY_DIR "replay"
#endif

int main(int argc, char **argv) {
  int failed = 0, cases = 0;
  for (int i = 1; i < argc; i++, cases++) failed |= replay(argv[i]) != 0;
  if (argc < 2) {
    DIR *d = opendir(ZPX_REPLAY_DIR);
    struct dirent *e;
    char path[4096];
    while (d != NULL && (e = readdir(d)) != NULL) {
      size_t len = strlen(e->d_name);
      if (len < 5 || strcmp(e->d_name + len - 4, ".txt") != 0) continue;
      snprintf(path, sizeof path, "%s/%s", ZPX_REPLAY_DIR, e->d_name);
      failed |= replay(path) != 0;
      cases++;
    }
    if (d != NULL) closedir(d);
  }
  /* Too few cases is a failure: a build that lost them proves nothing. */
  if (cases < 19) {
    fprintf(stderr, "FAIL: %d case files, expected at least 19\n", cases);
    failed = 1;
  }
  if (!failed) printf("replay: ok (%d cases; 1e-12 relative, 1e-8 for the adaptive solvers)\n", cases);
  return failed;
}
