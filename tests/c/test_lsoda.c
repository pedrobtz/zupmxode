/* liblsoda as vendored (src/vendor/liblsoda) is re-entrant per context
 * (design D5): two contexts stepped in turn give what each gives alone.
 * Built without R. Exit status 0 on success. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "lsoda.h"

static int decay(double t, double *y, double *ydot, void *data) {
  (void) t;
  const double *k = data;
  ydot[0] = -k[0] * y[0];
  ydot[1] = k[0] * y[0] - k[1] * y[1];
  return 0;
}

/* Solves from 0 to each of n times into out, stepping ctx along. */
static int run(struct lsoda_context_t *ctx, double *y, double *t, double tout) {
  lsoda(ctx, y, t, tout);
  return ctx->state > 0;
}

static int setup(struct lsoda_context_t *ctx, struct lsoda_opt_t *opt, double *rtol,
                 double *atol, double *k) {
  memset(ctx, 0, sizeof(*ctx));
  memset(opt, 0, sizeof(*opt));
  opt->itask = 1;
  opt->rtol = rtol;
  opt->atol = atol;
  ctx->function = decay;
  ctx->data = k;
  ctx->neq = 2;
  ctx->state = 1;
  return lsoda_prepare(ctx, opt);
}

int main(void) {
  double rtol[2] = {1e-10, 1e-10}, atol[2] = {1e-12, 1e-12};
  double ka[2] = {0.8, 0.1}, kb[2] = {2.5, 0.03};
  double alone_a[10][2], alone_b[10][2];
  struct lsoda_context_t a, b;
  struct lsoda_opt_t oa, ob;
  int failed = 0;

  for (int pass = 0; pass < 2; pass++) {
    double ya[2] = {100, 0}, yb[2] = {50, 0}, ta = 0, tb = 0;
    if (!setup(&a, &oa, rtol, atol, ka) || !setup(&b, &ob, rtol, atol, kb)) return 2;
    for (int i = 0; i < 10; i++) {
      double tout = i + 1;
      if (pass == 0) {
        /* Alone: a to the end, then b. */
        if (!run(&a, ya, &ta, tout)) failed = 1;
        memcpy(alone_a[i], ya, sizeof(ya));
      } else {
        /* In turn. */
        if (!run(&a, ya, &ta, tout) || !run(&b, yb, &tb, tout)) failed = 1;
        if (memcmp(ya, alone_a[i], sizeof(ya)) != 0 || memcmp(yb, alone_b[i], sizeof(yb)) != 0) {
          fprintf(stderr, "FAIL: interleaved contexts differ at t = %g\n", tout);
          failed = 1;
        }
      }
    }
    if (pass == 0) {
      for (int i = 0; i < 10; i++) {
        if (!run(&b, yb, &tb, i + 1.0)) failed = 1;
        memcpy(alone_b[i], yb, sizeof(yb));
      }
      /* Against the exact solution, to know the runs are right. */
      double t = 10, e = 100 * exp(-0.8 * t);
      if (fabs(alone_a[9][0] - e) > 1e-8 * e) {
        fprintf(stderr, "FAIL: lsoda off the exact solution\n");
        failed = 1;
      }
    }
    lsoda_free(&a);
    lsoda_free(&b);
  }
  if (!failed) printf("lsoda: ok (two contexts in turn, bit for bit)\n");
  return failed;
}
