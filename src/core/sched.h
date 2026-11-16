/* The event scheduler (design §4.2, §6.2): one subject's records become a
 * sorted timeline of what happens to the compartments and when. Builds
 * without R.
 *
 * A record is a row of the data: its time, EVID, AMT, compartment, RATE,
 * DUR, II, ADDL and SS, and what the model's `pre` block gave for its dose
 * compartment at that record (lag, bioavailability, modelled rate and
 * duration). The timeline holds one entry per thing that happens: an
 * observation, a bolus, the start and end of an infusion, a reset, a
 * steady-state dose, a compartment switched off or on, or a record with
 * nothing to do but mark a time (EVID 2, where covariates change).
 *
 * Order (decided 2026-11-16, design §6.2): by time; at one time infusion
 * ends first, then doses the scheduler generated (additional doses, doses
 * shifted by their lag), then the records themselves in data order, so
 * that an observation written before a dose at the same time is read
 * before it, as PREDPP reads a data set. */
#ifndef ZUPMXODE_SCHED_H
#define ZUPMXODE_SCHED_H

#include <stddef.h>

/* The most timeline entries one subject may have (assurance Limits row):
 * ADDL can ask for any number of doses. */
#define ZPX_SCHED_MAX_ENTRIES 1000000

enum zpx_entry_kind {
  ZPX_E_OBS = 0,       /* an observation record */
  ZPX_E_BOLUS = 1,     /* amount into cmt */
  ZPX_E_INF_START = 2, /* rate into cmt from here */
  ZPX_E_INF_END = 3,   /* that rate stops here */
  ZPX_E_RESET = 4,     /* every compartment to 0, every infusion stopped */
  ZPX_E_SS = 5,        /* steady state: see ss, amount, rate, ii */
  ZPX_E_OFF = 6,       /* cmt switched off (emptied, held at 0) */
  ZPX_E_ON = 7,        /* cmt switched back on */
  ZPX_E_MARK = 8       /* an EVID 2 record: a time where covariates change */
};

typedef struct zpx_record {
  double time, amt, rate, dur, ii;
  int evid, cmt, addl, ss; /* cmt is 1-based; negative switches off */
  /* From the model's pre block for this record's compartment: */
  double lag, f, mrate, mdur;
} zpx_record;

typedef struct zpx_entry {
  double time;
  int kind;
  int cmt;       /* 1-based; 0 when not about one compartment */
  double amount; /* bolus and steady-state amounts, after bioavailability */
  double rate;   /* infusion rate */
  double dur;    /* infusion duration (steady state too) */
  double ii;     /* steady-state interval */
  int ss;        /* 1 or 2 on ZPX_E_SS */
  int rec;       /* 0-based record the entry came from */
  int phase;     /* order at one time: 0 infusion end, 1 generated, 2 record */
  long seq;      /* order within a phase */
} zpx_entry;

/* Errors: 0 is success. *where is the 0-based record. */
enum {
  ZPX_SCHED_OK = 0,
  ZPX_SCHED_E_TIME,      /* TIME missing, not finite, or decreasing */
  ZPX_SCHED_E_EVID,      /* EVID not 0 to 4 */
  ZPX_SCHED_E_AMT,       /* AMT negative or not finite on a dose; nonzero on EVID 0, 2 or 3 */
  ZPX_SCHED_E_CMT,       /* a compartment the model does not have */
  ZPX_SCHED_E_RATE,      /* RATE other than >= 0, -1 or -2 */
  ZPX_SCHED_E_MRATE,     /* RATE -1 and the model gives no positive rate() */
  ZPX_SCHED_E_MDUR,      /* RATE -2 and the model gives no positive dur() */
  ZPX_SCHED_E_DUR,       /* DUR negative, or given with a RATE */
  ZPX_SCHED_E_ADDL,      /* ADDL negative, or ADDL without II > 0 */
  ZPX_SCHED_E_II,        /* II negative or not finite */
  ZPX_SCHED_E_SS,        /* SS not 0 to 2, on a record that is no dose, or
                            without II > 0 for a bolus */
  ZPX_SCHED_E_LAGF,      /* the model's lag negative, or bioavailability
                            negative, or either not finite */
  ZPX_SCHED_E_TOO_MANY,  /* more than ZPX_SCHED_MAX_ENTRIES entries */
  ZPX_SCHED_E_NO_MEMORY,
  ZPX_SCHED_E_INTERNAL   /* the counting and building passes disagree: a bug */
};

const char *zpx_sched_error_message(int error);

/* How many entries the records make, or -1 with *error and *where set.
 * Checks every record; allocates nothing. */
long zpx_sched_count(const zpx_record *rec, int n, int ncmt, int *error, int *where);

/* Writes the timeline into out[0..count), sorted, and returns count, or -1
 * with *error and *where. out must hold zpx_sched_count() entries. */
long zpx_sched_build(const zpx_record *rec, int n, int ncmt, zpx_entry *out, int *error,
                     int *where);

#endif
