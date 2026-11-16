/* One subject's solve (design §4.3, §4.4): the scheduler's timeline walked
 * with the closed form, the matrix exponential, LSODA or Dormand-Prince,
 * the model's blocks called through function pointers. Builds without R,
 * so the parallel runtime can run it on any thread: nothing here can touch
 * an R object or the R heap. Memory comes from a caller's arena; only
 * LSODA allocates, on the C heap, inside its own calls. */
#ifndef ZUPMXODE_SOLVE_H
#define ZUPMXODE_SOLVE_H

#include <stddef.h>

#include "ir.h"
#include "sched.h"

typedef void (*zpx_block)(double t, const double *y, double *out, const double *par,
                          double *work);

/* The model as the solver sees it. */
typedef struct zpx_blocks {
  zpx_block fn[ZPX_IR_NBLOCKS];
  int ny, npar, nwork, ncap;
} zpx_blocks;

enum { ZPX_M_CLOSED = 0, ZPX_M_LSODA = 1, ZPX_M_DOPRI = 2, ZPX_M_EXPM = 3 };

typedef struct zpx_solve_opts {
  int method;
  double rtol, atol, ss_tol;
  long ss_cap;
  /* linCmt()'s rate constants for the closed form: ncmt, depot, micro,
   * then 7 capture indices (-1 for none), as R/solve.R lincmt_spec(). */
  int lc[10];
} zpx_solve_opts;

/* What happened to one subject. */
typedef struct zpx_subject_result {
  int error;          /* 0, ZPX_SOLVE_E_EVENTS or ZPX_SOLVE_E_SOLVER */
  const char *message; /* a static string for the error */
  int error_row;      /* 0-based record of the subject */
  double failed_at;   /* a solver failure: the time, NaN if none */
  int failed_row;     /* its record, -1 if none */
  int ss_unconverged; /* steady states that hit the cap */
  int ss_row;         /* the last of them */
} zpx_subject_result;

enum { ZPX_SOLVE_OK = 0, ZPX_SOLVE_E_EVENTS = 1, ZPX_SOLVE_E_SOLVER = 2 };

/* The largest linear system the matrix exponential takes. */
#define ZPX_SOLVE_LINEAR_MAX 20

/* Bytes of arena a subject of n records making `entries` timeline
 * entries needs. */
size_t zpx_solve_arena_size(const zpx_blocks *b, int n, long entries);

/* Solves one subject. rec[0..n) are its records (their lag, f, mrate and
 * mdur are filled here from the pre block); par its parameters, record i's
 * parameter p at par[i + par_ld * p]. Writes, for
 * each of its observation records, the time, the compartments and the
 * outputs into out at rows row0, row0 + 1, ..., column c at
 * out[row + c * ld]; orow[row] is the record. The arena must hold
 * zpx_solve_arena_size() bytes. */
void zpx_solve_subject(const zpx_blocks *b, const zpx_solve_opts *o, zpx_record *rec, int n,
                       const double *par, long par_ld, double *out, long ld, int *orow, long row0,
                       void *arena, size_t arena_size, zpx_subject_result *res);

#endif
