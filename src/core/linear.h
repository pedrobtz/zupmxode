/* The closed-form solver (design §4.3, D21): the linear mammillary system
 * of linCmt() -- an optional depot, a central compartment and up to two
 * peripheral ones -- carried exactly from one event to the next. Builds
 * without R.
 *
 * With states y and constant infusion rates r over an interval of length
 * dt, dy/dt = A y + r has the solution
 *
 *   y(dt) = E(dt) y(0) + G(dt) r,   E = exp(A dt),  G = A^-1 (E - I),
 *
 * computed through A's eigenvalues, which are -ka for the depot and the
 * roots of the disposition's characteristic polynomial, in closed form for
 * one, two and three compartments. Steady state is the fixed point of one
 * dosing interval. */
#ifndef ZUPMXODE_LINEAR_H
#define ZUPMXODE_LINEAR_H

#define ZPX_LIN_MAX 4

typedef struct zpx_lin {
  int n;       /* states: ncmt, plus 1 with a depot */
  int depot;   /* whether state 0 is a depot */
  int taylor;  /* eigenvalues too close: the matrix exponential instead */
  double a[ZPX_LIN_MAX][ZPX_LIN_MAX];    /* the system matrix A */
  double mu[ZPX_LIN_MAX];                /* eigenvalues, all < 0 */
  double v[ZPX_LIN_MAX][ZPX_LIN_MAX];    /* eigenvectors, by column */
  double w[ZPX_LIN_MAX][ZPX_LIN_MAX];    /* v's inverse */
} zpx_lin;

/* Errors: 0 is success. */
enum {
  ZPX_LIN_OK = 0,
  ZPX_LIN_E_PARAM,     /* a rate constant negative, zero where it cannot be, or not finite */
  ZPX_LIN_E_SINGULAR   /* eigenvalues too close to separate (see linear.c) */
};

/* Sets L up for ncmt = 1 to 3 compartments, with a depot when depot is
 * nonzero: elimination k10, distribution k12, k21 (two or more), k13, k31
 * (three), absorption ka (with a depot). */
int zpx_lin_setup(zpx_lin *L, int ncmt, int depot, double ka, double k10, double k12,
                  double k21, double k13, double k31);

/* y <- E(dt) y + G(dt) r over an interval dt >= 0. */
void zpx_lin_step(const zpx_lin *L, double dt, double *y, const double *r);

/* The steady state just before a dose, for a dose repeated every ii > 0:
 * a bolus b (a vector of amounts), or an infusion of rate r for dur
 * (0 < dur <= ii). With ii = 0, the steady state of a constant infusion r.
 * Writes it to y. */
void zpx_lin_ss_bolus(const zpx_lin *L, double ii, const double *b, double *y);
void zpx_lin_ss_infusion(const zpx_lin *L, double ii, double dur, const double *r, double *y);
void zpx_lin_ss_constant(const zpx_lin *L, const double *r, double *y);

#endif
