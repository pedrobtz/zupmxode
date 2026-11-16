/* The matrix exponential (design §4.3): Pade of degree 13 with scaling and
 * squaring, Higham (2005), "The scaling and squaring method for the matrix
 * exponential revisited", SIAM J. Matrix Anal. Appl. 26(4). For the linear
 * solver's systems of up to ZPX_EXPM_MAX - 1 states and their augmented
 * matrices. Builds without R. */
#ifndef ZUPMXODE_EXPM_H
#define ZUPMXODE_EXPM_H

/* The largest matrix: the linear solver's state limit (ZPX_LINEAR_MAX_STATES,
 * 20) plus the row of the augmented system. */
#define ZPX_EXPM_MAX 21

/* X = exp(A), both m x m with m <= ZPX_EXPM_MAX, row-major. work holds
 * 8 * m * m doubles. Returns 0, or 1 when A is not finite. */
int zpx_expm(int m, const double *A, double *X, double *work);

#endif
