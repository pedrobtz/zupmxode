/* Pade 13 with scaling and squaring (src/core/expm.h). */
#include "expm.h"

#include <math.h>
#include <string.h>

/* C = A B, m x m, row-major. */
static void matmul(int m, const double *A, const double *B, double *C) {
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
      double s = 0;
      for (int k = 0; k < m; k++) s += A[i * m + k] * B[k * m + j];
      C[i * m + j] = s;
    }
}

static double norm1(int m, const double *A) {
  double best = 0;
  for (int j = 0; j < m; j++) {
    double s = 0;
    for (int i = 0; i < m; i++) s += fabs(A[i * m + j]);
    if (s > best) best = s;
  }
  return best;
}

/* Solves Q X = P for X, m x m, by elimination with partial pivoting; Q and
 * P are overwritten, X written to P. Returns 0 when Q is singular. */
static int solve_matrix(int m, double *Q, double *P) {
  for (int c = 0; c < m; c++) {
    int p = c;
    for (int r = c + 1; r < m; r++)
      if (fabs(Q[r * m + c]) > fabs(Q[p * m + c])) p = r;
    if (Q[p * m + c] == 0) return 0;
    if (p != c)
      for (int k = 0; k < m; k++) {
        double t = Q[c * m + k];
        Q[c * m + k] = Q[p * m + k];
        Q[p * m + k] = t;
        t = P[c * m + k];
        P[c * m + k] = P[p * m + k];
        P[p * m + k] = t;
      }
    for (int r = c + 1; r < m; r++) {
      double f = Q[r * m + c] / Q[c * m + c];
      if (f == 0) continue;
      for (int k = c; k < m; k++) Q[r * m + k] -= f * Q[c * m + k];
      for (int k = 0; k < m; k++) P[r * m + k] -= f * P[c * m + k];
    }
  }
  for (int c = m - 1; c >= 0; c--)
    for (int k = 0; k < m; k++) {
      double s = P[c * m + k];
      for (int j = c + 1; j < m; j++) s -= Q[c * m + j] * P[j * m + k];
      P[c * m + k] = s / Q[c * m + c];
    }
  return 1;
}

int zpx_expm(int m, const double *A, double *X, double *work) {
  static const double b[14] = {
    64764752532480000.0, 32382376266240000.0, 7771770303897600.0, 1187353796428800.0,
    129060195264000.0,   10559470521600.0,    670442572800.0,     33522128640.0,
    1323241920.0,        40840800.0,          960960.0,           16380.0,
    182.0,               1.0};
  const double theta13 = 5.371920351148152;
  int mm = m * m;
  double *As = work, *A2 = work + mm, *A4 = work + 2 * mm, *A6 = work + 3 * mm,
         *U = work + 4 * mm, *V = work + 5 * mm, *T = work + 6 * mm, *Q = work + 7 * mm;
  for (int i = 0; i < mm; i++)
    if (!isfinite(A[i])) return 1;
  double nrm = norm1(m, A);
  int s = 0;
  if (nrm > theta13) s = (int) ceil(log2(nrm / theta13));
  if (s < 0) s = 0;
  double scale = ldexp(1.0, -s);
  for (int i = 0; i < mm; i++) As[i] = A[i] * scale;
  matmul(m, As, As, A2);
  matmul(m, A2, A2, A4);
  matmul(m, A4, A2, A6);
  /* U = A [A6 (b13 A6 + b11 A4 + b9 A2) + b7 A6 + b5 A4 + b3 A2 + b1 I] */
  for (int i = 0; i < mm; i++) T[i] = b[13] * A6[i] + b[11] * A4[i] + b[9] * A2[i];
  matmul(m, A6, T, U);
  for (int i = 0; i < mm; i++) U[i] += b[7] * A6[i] + b[5] * A4[i] + b[3] * A2[i];
  for (int i = 0; i < m; i++) U[i * m + i] += b[1];
  memcpy(T, U, sizeof(double) * (size_t) mm);
  matmul(m, As, T, U);
  /* V = A6 (b12 A6 + b10 A4 + b8 A2) + b6 A6 + b4 A4 + b2 A2 + b0 I */
  for (int i = 0; i < mm; i++) T[i] = b[12] * A6[i] + b[10] * A4[i] + b[8] * A2[i];
  matmul(m, A6, T, V);
  for (int i = 0; i < mm; i++) V[i] += b[6] * A6[i] + b[4] * A4[i] + b[2] * A2[i];
  for (int i = 0; i < m; i++) V[i * m + i] += b[0];
  /* (V - U) X = (V + U) */
  for (int i = 0; i < mm; i++) {
    Q[i] = V[i] - U[i];
    X[i] = V[i] + U[i];
  }
  if (!solve_matrix(m, Q, X)) return 1;
  for (int k = 0; k < s; k++) {
    matmul(m, X, X, T);
    memcpy(X, T, sizeof(double) * (size_t) mm);
  }
  return 0;
}
