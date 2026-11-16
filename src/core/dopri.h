/* Dormand-Prince 5(4) (design §4.3, D5), written for the package, with
 * dense output: Hairer, Norsett and Wanner, Solving Ordinary Differential
 * Equations I, II.5 (the coefficients and the continuous extension of
 * their dopri5). Builds without R. */
#ifndef ZUPMXODE_DOPRI_H
#define ZUPMXODE_DOPRI_H

typedef void (*zpx_rhs)(double t, const double *y, double *dydt, void *data);

typedef struct zpx_dopri {
  int n;
  double rtol, atol;
  long max_steps; /* per call */
  double h;       /* the next step size, carried between calls; 0 to choose */
  long steps;     /* taken in the last call */
  double t;       /* where the last call got to: t1, or where it failed */
  int nout_done;  /* dense outputs the last call wrote */
} zpx_dopri;

/* A failure leaves y at the last accepted step, D->t its time, and the
 * dense outputs before it written (D->nout_done). */
enum { ZPX_DOPRI_OK = 0, ZPX_DOPRI_E_STEPS, ZPX_DOPRI_E_STEP_SIZE, ZPX_DOPRI_E_NONFINITE };

/* Doubles of workspace zpx_dopri_solve() needs for n states. */
#define ZPX_DOPRI_WORK(n) (10 * (n))

/* Integrates y from t0 to t1 > t0. Each of the nout times in tout, in
 * increasing order within (t0, t1], gets the solution there by dense
 * output, written to yout[k * n ...]. y is the solution at t1 after. */
int zpx_dopri_solve(zpx_dopri *D, zpx_rhs f, void *data, double t0, double t1, double *y,
                    const double *tout, int nout, double *yout, double *work);

#endif
