/* zupmxode's C API (design §4.6, D9, D12): one subject's predictions and
 * residual variances at one parameter vector, written into the caller's
 * buffers, for estimators.
 *
 * Plain C99: this header needs no R header. A package that links to
 * zupmxode (LinkingTo: zupmxode, and zupmxode in Imports so its namespace
 * is loaded) includes R_ext/Rdynload.h first and calls
 * zupmxode_api_load(), which fetches the function table registered with
 * R_RegisterCCallable(). Everything else goes through that table.
 *
 * Threads: zupmxode_subject_eval() is re-entrant and keeps no global
 * state; two threads may evaluate subjects of the same model at once, each
 * with its own arena. It never calls R, never raises an R error and never
 * allocates on R's heap: it returns a status code. LSODA, the one solver
 * that allocates, does so on the C heap within the call. The other
 * entries (arena_new, model_from_r, subject_from_r) run on R's main
 * thread.
 *
 * Lifetimes: a model and its subjects are owned by their R objects. Keep
 * those alive (protected or reachable) for as long as an arena made for
 * them or a pointer taken from them is used. A subject holds copies of its
 * records and points into no R object. */
#ifndef ZUPMXODE_H
#define ZUPMXODE_H

#include <stddef.h>

/* The version of this header's table and of zupmxode_eval_opts. */
#define ZUPMXODE_API_VERSION 1
#define ZUPMXODE_EVAL_OPTS_VERSION 1

typedef struct zupmxode_model zupmxode_model;     /* a compiled model */
typedef struct zupmxode_subject zupmxode_subject; /* a prepared subject */
typedef struct zupmxode_arena zupmxode_arena;     /* one thread's workspace */

/* Status codes. Below ZUPMXODE_E_FIRST_ERROR, pred and var are written. */
enum zupmxode_status {
  ZUPMXODE_OK = 0,
  ZUPMXODE_SOLVER_FAILED = 1,  /* pred and var NaN from where it failed */
  ZUPMXODE_SS_UNCONVERGED = 2, /* a steady state hit its cap; the rest is solved */
  ZUPMXODE_E_FIRST_ERROR = 16,
  ZUPMXODE_E_VERSION = 16,  /* opts->version is not one this build knows */
  ZUPMXODE_E_ARGUMENT = 17, /* a NULL pointer, an output, error model or
                               method the model and subject do not have,
                               or a subject prepared for another model */
  ZUPMXODE_E_BUFFER = 18,   /* opts->n_pred is below the subject's observations */
  ZUPMXODE_E_ARENA = 19,    /* the arena is for another model, or too small */
  ZUPMXODE_E_EVENTS = 20,   /* a record breaks an event rule at these parameters */
  ZUPMXODE_E_MODEL = 21     /* the closed form's rate constants are not positive */
};

/* Residual error models (design D19), each sd a standard deviation,
 * computed as written: additive var = add^2; proportional
 * var = prop^2 * (f * f); combined var = add^2 + prop^2 * (f * f);
 * exponential (lnorm) var = add^2, on the log scale. */
enum zupmxode_error_model {
  ZUPMXODE_ERR_ADD = 1,
  ZUPMXODE_ERR_PROP = 2,
  ZUPMXODE_ERR_COMBINED = 3,
  ZUPMXODE_ERR_LNORM = 4
};

/* Solvers, for zupmxode_eval_opts.method. */
enum zupmxode_method {
  ZUPMXODE_METHOD_DEFAULT = -1, /* the subject's, chosen when it was prepared */
  ZUPMXODE_METHOD_CLOSED = 0,
  ZUPMXODE_METHOD_LSODA = 1,
  ZUPMXODE_METHOD_DOPRI5 = 2,
  ZUPMXODE_METHOD_EXPM = 3
};

/* Options of one call. Set version to ZUPMXODE_EVAL_OPTS_VERSION; a later
 * version adds fields at the end (sensitivity buffers, design §14) and
 * this build refuses a version it does not know. */
typedef struct zupmxode_eval_opts {
  int version;
  int output;       /* 0-based: the compartments, then the model's outputs */
  int error_model;  /* enum zupmxode_error_model */
  double add, prop; /* standard deviations; lnorm takes add */
  int method;       /* enum zupmxode_method */
  double rtol, atol; /* the general solvers' tolerances; <= 0: the subject's */
  size_t n_pred;    /* the length of pred and of var */
} zupmxode_eval_opts;

typedef struct zupmxode_api {
  int version; /* ZUPMXODE_API_VERSION of the build */
  /* pred and var get one value per observation record of s, in record
   * order. theta and eta are as the subject was prepared: its theta names
   * and its eta names, in that order. */
  int (*subject_eval)(const zupmxode_model *m, const zupmxode_subject *s, const double *theta,
                      const double *eta, double *pred, double *var,
                      const zupmxode_eval_opts *opts, zupmxode_arena *arena);
  /* An arena large enough for each of the n subjects, on R's main thread;
   * NULL when out of memory or given no subject of m. */
  zupmxode_arena *(*arena_new)(const zupmxode_model *m, const zupmxode_subject *const *s, int n);
  void (*arena_free)(zupmxode_arena *arena);
  int (*subject_nobs)(const zupmxode_subject *s);
  /* The model or subject behind an R object (a SEXP, passed as void *);
   * NULL if it is not one, or was saved and not prepared again. */
  const zupmxode_model *(*model_from_r)(void *sexp);
  const zupmxode_subject *(*subject_from_r)(void *sexp);
  /* A static description of a status code. */
  const char *(*status_message)(int status);
} zupmxode_api;

#ifdef R_EXT_DYNLOAD_H_
/* The table, or NULL when zupmxode's build is older than this header. */
static inline const zupmxode_api *zupmxode_api_load(void) {
  const zupmxode_api *(*get)(int) =
      (const zupmxode_api *(*)(int)) R_GetCCallable("zupmxode", "zupmxode_api");
  return get(ZUPMXODE_API_VERSION);
}
#endif

#endif
