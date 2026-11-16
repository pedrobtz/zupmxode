/* Dormand-Prince 5(4) with dense output (src/core/dopri.h). */
#include "dopri.h"

#include <math.h>
#include <string.h>

static const double c2 = 1.0 / 5, c3 = 3.0 / 10, c4 = 4.0 / 5, c5 = 8.0 / 9;
static const double a21 = 1.0 / 5;
static const double a31 = 3.0 / 40, a32 = 9.0 / 40;
static const double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
static const double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561,
                    a54 = -212.0 / 729;
static const double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247,
                    a64 = 49.0 / 176, a65 = -5103.0 / 18656;
static const double a71 = 35.0 / 384, a73 = 500.0 / 1113, a74 = 125.0 / 192,
                    a75 = -2187.0 / 6784, a76 = 11.0 / 84;
static const double e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920,
                    e5 = -17253.0 / 339200, e6 = 22.0 / 525, e7 = -1.0 / 40;
static const double d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
                    d4 = -10690763975.0 / 1880347072.0, d5 = 701980252875.0 / 199316789632.0,
                    d6 = -1453857185.0 / 822651844.0, d7 = 69997945.0 / 29380423.0;

/* The weighted RMS norm of x against the scale atol + rtol max(|a|, |b|). */
static double wnorm(int n, const double *x, const double *a, const double *b, double rtol,
                    double atol) {
  double s = 0;
  for (int i = 0; i < n; i++) {
    double sc = atol + rtol * fmax(fabs(a[i]), fabs(b[i]));
    double q = x[i] / sc;
    s += q * q;
  }
  return sqrt(s / (n > 0 ? n : 1));
}

/* Hairer's initial step: from the size of y and f, then one Euler step. */
static double initial_step(int n, zpx_rhs f, void *data, double t, const double *y,
                           const double *k1, double rtol, double atol, double span, double *tmp,
                           double *f1) {
  double d0 = wnorm(n, y, y, y, rtol, atol), d1n = wnorm(n, k1, y, y, rtol, atol);
  double h0 = (d0 < 1e-5 || d1n < 1e-5) ? 1e-6 : 0.01 * d0 / d1n;
  if (h0 > span) h0 = span;
  for (int i = 0; i < n; i++) tmp[i] = y[i] + h0 * k1[i];
  f(t + h0, tmp, f1, data);
  for (int i = 0; i < n; i++) tmp[i] = (f1[i] - k1[i]) / h0;
  double d2 = wnorm(n, tmp, y, y, rtol, atol);
  double m = fmax(d1n, d2);
  double h1 = m <= 1e-15 ? fmax(1e-6, h0 * 1e-3) : pow(0.01 / m, 1.0 / 5);
  double h = fmin(100 * h0, h1);
  return h > span ? span : h;
}

int zpx_dopri_solve(zpx_dopri *D, zpx_rhs f, void *data, double t0, double t1, double *y,
                    const double *tout, int nout, double *yout, double *work) {
  int n = D->n;
  double *k1 = work, *k2 = work + n, *k3 = work + 2 * n, *k4 = work + 3 * n, *k5 = work + 4 * n,
         *k6 = work + 5 * n, *k7 = work + 6 * n, *yt = work + 7 * n, *y1 = work + 8 * n,
         *err = work + 9 * n;
  double t = t0, span = t1 - t0;
  int next = 0;
  D->steps = 0;
  D->t = t0;
  D->nout_done = 0;
  if (!(span > 0)) return ZPX_DOPRI_OK;
  f(t, y, k1, data);
  double h = D->h > 0 ? fmin(D->h, span) : initial_step(n, f, data, t, y, k1, D->rtol, D->atol,
                                                        span, yt, k2);
  int last_rejected = 0;
  while (t < t1) {
    D->t = t;
    D->nout_done = next;
    if (D->steps >= D->max_steps) return ZPX_DOPRI_E_STEPS;
    int final = 0;
    if (t + h >= t1 || t + 1.01 * h >= t1) {
      h = t1 - t;
      final = 1;
    }
    if (!(h > fabs(t) * 1e-15) && !final) return ZPX_DOPRI_E_STEP_SIZE;
    for (int i = 0; i < n; i++) yt[i] = y[i] + h * a21 * k1[i];
    f(t + c2 * h, yt, k2, data);
    for (int i = 0; i < n; i++) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
    f(t + c3 * h, yt, k3, data);
    for (int i = 0; i < n; i++) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    f(t + c4 * h, yt, k4, data);
    for (int i = 0; i < n; i++)
      yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    f(t + c5 * h, yt, k5, data);
    for (int i = 0; i < n; i++)
      yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    double tn = final ? t1 : t + h;
    f(tn, yt, k6, data);
    for (int i = 0; i < n; i++)
      y1[i] = y[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    f(tn, y1, k7, data);
    D->steps++;
    for (int i = 0; i < n; i++)
      err[i] = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
    double en = wnorm(n, err, y, y1, D->rtol, D->atol);
    for (int i = 0; i < n; i++)
      if (!isfinite(y1[i])) en = INFINITY;
    double fac = en == 0 ? 10 : 0.9 * pow(en, -1.0 / 5);
    if (fac > 10) fac = 10;
    if (fac < 0.2) fac = 0.2;
    if (en <= 1.0) {
      /* Accepted: dense output for every requested time in (t, tn]. */
      while (next < nout && tout[next] <= tn) {
        double th = (tout[next] - t) / h;
        double *o = yout + (long) next * n;
        for (int i = 0; i < n; i++) {
          double ydiff = y1[i] - y[i];
          double bspl = h * k1[i] - ydiff;
          double r4 = ydiff - h * k7[i] - bspl;
          double r5 = h * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] +
                           d7 * k7[i]);
          o[i] = y[i] + th * (ydiff + (1 - th) * (bspl + th * (r4 + (1 - th) * r5)));
        }
        next++;
      }
      memcpy(y, y1, sizeof(double) * (size_t) n);
      memcpy(k1, k7, sizeof(double) * (size_t) n); /* first same as last */
      t = tn;
      if (last_rejected && fac > 1) fac = 1;
      last_rejected = 0;
      if (!final) h *= fac;
    } else {
      if (!isfinite(en) && h < fabs(t) * 1e-12) return ZPX_DOPRI_E_NONFINITE;
      h *= fac;
      last_rejected = 1;
    }
  }
  D->h = h;
  D->t = t1;
  for (; next < nout; next++) memcpy(yout + (long) next * n, y, sizeof(double) * (size_t) n);
  D->nout_done = nout;
  return ZPX_DOPRI_OK;
}
