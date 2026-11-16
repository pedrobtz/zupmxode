/* The closed-form solver (src/core/linear.h). Builds without R. */
#include "linear.h"

#include <math.h>
#include <string.h>

/* ---- Small dense algebra ------------------------------------------------ */

/* Solves M x = b for n <= ZPX_LIN_MAX + 1 by Gaussian elimination with
 * partial pivoting; M and b are overwritten. Returns 0 when M is singular. */
static int solve(int n, double M[][ZPX_LIN_MAX + 1], double *b) {
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int r = c + 1; r < n; r++)
      if (fabs(M[r][c]) > fabs(M[p][c])) p = r;
    if (M[p][c] == 0.0) return 0;
    if (p != c) {
      for (int k = 0; k < n; k++) {
        double t = M[c][k];
        M[c][k] = M[p][k];
        M[p][k] = t;
      }
      double t = b[c];
      b[c] = b[p];
      b[p] = t;
    }
    for (int r = c + 1; r < n; r++) {
      double f = M[r][c] / M[c][c];
      for (int k = c; k < n; k++) M[r][k] -= f * M[c][k];
      b[r] -= f * b[c];
    }
  }
  for (int c = n - 1; c >= 0; c--) {
    double s = b[c];
    for (int k = c + 1; k < n; k++) s -= M[c][k] * b[k];
    b[c] = s / M[c][c];
  }
  return 1;
}

/* A unit null vector of (A - mu I), n <= 4, by elimination with full
 * pivoting: the column left without a pivot is set to 1. Returns 0 if the
 * null space is not one-dimensional to working precision. */
static int null_vector(int n, double A[ZPX_LIN_MAX][ZPX_LIN_MAX], double mu, double *x) {
  double M[ZPX_LIN_MAX][ZPX_LIN_MAX];
  int col[ZPX_LIN_MAX];
  double scale = 0;
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) {
      M[r][c] = A[r][c] - (r == c ? mu : 0.0);
      if (fabs(M[r][c]) > scale) scale = fabs(M[r][c]);
    }
  for (int c = 0; c < n; c++) col[c] = c;
  int rank = 0;
  for (int k = 0; k < n; k++) {
    int pr = -1, pc = -1;
    double best = 0;
    for (int r = k; r < n; r++)
      for (int c = k; c < n; c++)
        if (fabs(M[r][c]) > best) {
          best = fabs(M[r][c]);
          pr = r;
          pc = c;
        }
    if (best <= 1e-12 * scale) break;
    if (pr != k)
      for (int c = 0; c < n; c++) {
        double t = M[k][c];
        M[k][c] = M[pr][c];
        M[pr][c] = t;
      }
    if (pc != k) {
      for (int r = 0; r < n; r++) {
        double t = M[r][k];
        M[r][k] = M[r][pc];
        M[r][pc] = t;
      }
      int t = col[k];
      col[k] = col[pc];
      col[pc] = t;
    }
    for (int r = k + 1; r < n; r++) {
      double f = M[r][k] / M[k][k];
      for (int c = k; c < n; c++) M[r][c] -= f * M[k][c];
    }
    rank++;
  }
  if (rank != n - 1) return 0;
  double y[ZPX_LIN_MAX];
  y[n - 1] = 1.0;
  for (int k = n - 2; k >= 0; k--) {
    double s = 0;
    for (int c = k + 1; c < n; c++) s -= M[k][c] * y[c];
    y[k] = s / M[k][k];
  }
  double norm = 0;
  for (int k = 0; k < n; k++) norm += y[k] * y[k];
  norm = sqrt(norm);
  for (int k = 0; k < n; k++) x[col[k]] = y[k] / norm;
  return 1;
}

/* ---- Eigenvalues -------------------------------------------------------- */

/* The roots of x^3 + a2 x^2 + a1 x + a0 for the three-compartment
 * disposition, which are real and negative: Viete's trigonometric form,
 * each refined by Newton's method on the polynomial. */
static void cubic_roots(double a2, double a1, double a0, double *r) {
  double q = (3 * a1 - a2 * a2) / 9;
  double s = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
  double m = sqrt(-q);
  double c = s / (m * m * m);
  if (c > 1) c = 1;
  if (c < -1) c = -1;
  double th = acos(c);
  const double pi = 3.14159265358979323846;
  for (int k = 0; k < 3; k++) {
    double x = 2 * m * cos((th + 2 * pi * k) / 3) - a2 / 3;
    for (int it = 0; it < 4; it++) {
      double f = ((x + a2) * x + a1) * x + a0;
      double d = (3 * x + 2 * a2) * x + a1;
      if (d == 0) break;
      x -= f / d;
    }
    r[k] = x;
  }
}

static int disposition_eigen(int ncmt, double k10, double k12, double k21, double k13,
                             double k31, double *lam) {
  if (ncmt == 1) {
    lam[0] = -k10;
    return 1;
  }
  if (ncmt == 2) {
    double b = k10 + k12 + k21, c = k10 * k21;
    double disc = b * b - 4 * c;
    if (disc < 0) disc = 0;
    /* The larger root without cancellation, the smaller from the product. */
    double big = -(b + sqrt(disc)) / 2;
    lam[0] = big;
    lam[1] = c / big;
    return 1;
  }
  double a2 = k10 + k12 + k13 + k21 + k31;
  double a1 = k10 * k21 + k10 * k31 + k12 * k31 + k13 * k21 + k21 * k31;
  double a0 = k10 * k21 * k31;
  cubic_roots(a2, a1, a0, lam);
  return 1;
}

/* ---- Setup ---------------------------------------------------------------- */

static int ok_rate(double k, int may_be_zero) {
  return isfinite(k) && (k > 0 || (may_be_zero && k == 0));
}

int zpx_lin_setup(zpx_lin *L, int ncmt, int depot, double ka, double k10, double k12,
                  double k21, double k13, double k31) {
  memset(L, 0, sizeof(*L));
  if (ncmt < 1 || ncmt > 3) return ZPX_LIN_E_PARAM;
  if (!ok_rate(k10, 1) || (depot && !ok_rate(ka, 0))) return ZPX_LIN_E_PARAM;
  if (ncmt >= 2 && (!ok_rate(k12, 0) || !ok_rate(k21, 0))) return ZPX_LIN_E_PARAM;
  if (ncmt >= 3 && (!ok_rate(k13, 0) || !ok_rate(k31, 0))) return ZPX_LIN_E_PARAM;
  int n = ncmt + (depot != 0);
  int c0 = depot != 0; /* the central compartment's state */
  L->n = n;
  L->depot = depot != 0;
  memset(L->a, 0, sizeof(L->a));
  if (depot) {
    L->a[0][0] = -ka;
    L->a[c0][0] = ka;
  }
  L->a[c0][c0] = -(k10 + (ncmt >= 2 ? k12 : 0) + (ncmt >= 3 ? k13 : 0));
  if (ncmt >= 2) {
    L->a[c0][c0 + 1] = k21;
    L->a[c0 + 1][c0] = k12;
    L->a[c0 + 1][c0 + 1] = -k21;
  }
  if (ncmt >= 3) {
    L->a[c0][c0 + 2] = k31;
    L->a[c0 + 2][c0] = k13;
    L->a[c0 + 2][c0 + 2] = -k31;
  }
  double lam[3];
  disposition_eigen(ncmt, k10, k12, k21, k13, k31, lam);
  int m = 0;
  if (depot) L->mu[m++] = -ka;
  for (int i = 0; i < ncmt; i++) L->mu[m++] = lam[i];
  /* Eigenvalues too close to separate (ka at an elimination rate, or a
   * zero elimination rate next to nothing): the eigenvector basis is
   * ill-conditioned there, so the matrix exponential takes over. */
  L->taylor = 0;
  for (int i = 0; i < n && !L->taylor; i++) {
    if (L->mu[i] == 0) L->taylor = 1;
    for (int j = i + 1; j < n; j++)
      if (fabs(L->mu[i] - L->mu[j]) <= 1e-6 * fmax(fabs(L->mu[i]), fabs(L->mu[j])))
        L->taylor = 1;
  }
  if (!L->taylor) {
    for (int i = 0; i < n; i++) {
      double x[ZPX_LIN_MAX];
      if (!null_vector(n, L->a, L->mu[i], x)) {
        L->taylor = 1;
        break;
      }
      for (int r = 0; r < n; r++) L->v[r][i] = x[r];
    }
  }
  if (!L->taylor) {
    /* w = v^-1, column by column. */
    for (int c = 0; c < n; c++) {
      double M[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1], b[ZPX_LIN_MAX + 1];
      for (int r = 0; r < n; r++) {
        for (int k = 0; k < n; k++) M[r][k] = L->v[r][k];
        b[r] = r == c ? 1.0 : 0.0;
      }
      if (!solve(n, M, b)) {
        L->taylor = 1;
        break;
      }
      for (int r = 0; r < n; r++) L->w[r][c] = b[r];
    }
  }
  return ZPX_LIN_OK;
}

/* ---- The matrix exponential, for the confluent cases --------------------- */

/* exp(B) for the (n+1) x (n+1) matrix B = [[A dt, r dt], [0, 0]], whose
 * top-right column is G(dt) r and top-left block E(dt): Taylor's series
 * after scaling by 2^-s, then squaring s times. */
static void expm_augmented(const zpx_lin *L, double dt, const double *r,
                           double X[][ZPX_LIN_MAX + 1]) {
  int n = L->n, m = n + 1;
  double B[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1] = {{0}};
  double norm = 0;
  for (int i = 0; i < n; i++) {
    double row = 0;
    for (int j = 0; j < n; j++) {
      B[i][j] = L->a[i][j] * dt;
      row += fabs(B[i][j]);
    }
    B[i][n] = (r != NULL ? r[i] : 0.0) * dt;
    row += fabs(B[i][n]);
    if (row > norm) norm = row;
  }
  int s = 0;
  while (norm > 0.25 && s < 1000) {
    norm /= 2;
    s++;
  }
  double scale = ldexp(1.0, -s);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) B[i][j] *= scale;
  double T[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1], P[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1];
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) X[i][j] = T[i][j] = i == j ? 1.0 : 0.0;
  for (int k = 1; k <= 24; k++) {
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        double acc = 0;
        for (int l = 0; l < m; l++) acc += T[i][l] * B[l][j];
        P[i][j] = acc / k;
      }
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        T[i][j] = P[i][j];
        X[i][j] += T[i][j];
      }
  }
  for (int q = 0; q < s; q++) {
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        double acc = 0;
        for (int l = 0; l < m; l++) acc += X[i][l] * X[l][j];
        P[i][j] = acc;
      }
    memcpy(X, P, sizeof(P));
  }
}

/* ---- Stepping and steady state ------------------------------------------ */

/* (e^x - 1) / x, exact near 0. */
static double phi1(double x) { return x == 0 ? 1.0 : expm1(x) / x; }

void zpx_lin_step(const zpx_lin *L, double dt, double *y, const double *r) {
  int n = L->n;
  double out[ZPX_LIN_MAX] = {0};
  if (dt <= 0) return;
  if (L->taylor) {
    double X[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1];
    expm_augmented(L, dt, r, X);
    for (int i = 0; i < n; i++) {
      double acc = X[i][n];
      for (int j = 0; j < n; j++) acc += X[i][j] * y[j];
      out[i] = acc;
    }
    memcpy(y, out, sizeof(double) * (size_t) n);
    return;
  }
  /* In the eigenbasis: z = w y, z_i <- e^{mu_i dt} z_i + dt phi1(mu_i dt) (w r)_i. */
  double z[ZPX_LIN_MAX];
  for (int i = 0; i < n; i++) {
    double zy = 0, zr = 0;
    for (int j = 0; j < n; j++) {
      zy += L->w[i][j] * y[j];
      if (r != NULL) zr += L->w[i][j] * r[j];
    }
    double x = L->mu[i] * dt;
    z[i] = exp(x) * zy + dt * phi1(x) * zr;
  }
  for (int i = 0; i < n; i++) {
    double acc = 0;
    for (int j = 0; j < n; j++) acc += L->v[i][j] * z[j];
    out[i] = acc;
  }
  memcpy(y, out, sizeof(double) * (size_t) n);
}

/* Solves (I - E(ii)) y = c for y, through the matrix exponential, in the
 * confluent cases. */
static void ss_solve(const zpx_lin *L, double ii, const double *c, double *y) {
  int n = L->n;
  double X[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1];
  expm_augmented(L, ii, NULL, X);
  double M[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1], b[ZPX_LIN_MAX + 1];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) M[i][j] = (i == j ? 1.0 : 0.0) - X[i][j];
    b[i] = c[i];
  }
  if (!solve(n, M, b))
    for (int i = 0; i < n; i++) b[i] = NAN;
  memcpy(y, b, sizeof(double) * (size_t) n);
}

void zpx_lin_ss_bolus(const zpx_lin *L, double ii, const double *b, double *y) {
  int n = L->n;
  if (L->taylor) {
    /* Trough y = (I - E)^-1 E b. */
    double c[ZPX_LIN_MAX];
    memcpy(c, b, sizeof(double) * (size_t) n);
    double zero[ZPX_LIN_MAX] = {0};
    zpx_lin_step(L, ii, c, zero);
    ss_solve(L, ii, c, y);
    return;
  }
  double z[ZPX_LIN_MAX];
  for (int i = 0; i < n; i++) {
    double zb = 0;
    for (int j = 0; j < n; j++) zb += L->w[i][j] * b[j];
    /* e^x / (1 - e^x) = 1 / (e^-x - 1). */
    z[i] = zb / expm1(-L->mu[i] * ii);
  }
  for (int i = 0; i < n; i++) {
    double acc = 0;
    for (int j = 0; j < n; j++) acc += L->v[i][j] * z[j];
    y[i] = acc;
  }
}

void zpx_lin_ss_infusion(const zpx_lin *L, double ii, double dur, const double *r, double *y) {
  int n = L->n;
  if (L->taylor) {
    /* Trough y = (I - E(ii))^-1 E(ii - dur) G(dur) r. */
    double c[ZPX_LIN_MAX] = {0}, zero[ZPX_LIN_MAX] = {0};
    zpx_lin_step(L, dur, c, r);
    zpx_lin_step(L, ii - dur, c, zero);
    ss_solve(L, ii, c, y);
    return;
  }
  double z[ZPX_LIN_MAX];
  for (int i = 0; i < n; i++) {
    double zr = 0;
    for (int j = 0; j < n; j++) zr += L->w[i][j] * r[j];
    double mu = L->mu[i];
    z[i] = exp(mu * (ii - dur)) * dur * phi1(mu * dur) * zr / -expm1(mu * ii);
  }
  for (int i = 0; i < n; i++) {
    double acc = 0;
    for (int j = 0; j < n; j++) acc += L->v[i][j] * z[j];
    y[i] = acc;
  }
}

void zpx_lin_ss_constant(const zpx_lin *L, const double *r, double *y) {
  int n = L->n;
  /* A y + r = 0. */
  double M[ZPX_LIN_MAX + 1][ZPX_LIN_MAX + 1], b[ZPX_LIN_MAX + 1];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) M[i][j] = L->a[i][j];
    b[i] = -r[i];
  }
  if (!solve(n, M, b))
    for (int i = 0; i < n; i++) b[i] = INFINITY;
  memcpy(y, b, sizeof(double) * (size_t) n);
}
