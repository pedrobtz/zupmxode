# RFC 0047: `zupmxode`, PK/PD simulation without a compiler

Adopted from [RFC 0047](https://github.com/pedrobtz/packages/blob/72859143843cf3d7407bcbcaec2fd41d677978ae/rfcs/not-started/0047-zupmxode-event-driven-ode-engine.md)
on 2026-11-15. This file is the single source of truth for the design; the RFC
is now a pointer here.

- **Status:** Proposed. Merging this file adopts the design as the starting
  point for the package's own `.agents/design.md`; each decision in section
  13 stays the maintainer's until the package repository records it.
  D2, models compiled to native code by MIR through `zumir`, is the
  decision the whole package rests on. It is provisional on the MIR spike
  that is `zumir`'s Stage 0 (RFC 0073 §10), and §11.1 says what happens if
  the spike fails.
- **Date:** 2026-10-28; amended 2026-11-12 for the rename from `zufaode`
  (D10), the per-subject evaluation entry point that RFC 0072
  `zupmxfit` needs (§4.6, D12), and MIR in place of a hand-written
  bytecode machine (D2, D13, D14). [`foce.md`](https://github.com/pedrobtz/packages/blob/main/foce.md) gives the
  order of the packages that lead to an estimator.
- **Packages:** `zupmxode` (not yet created; crandb returns `not_found`
  and there is no `pedrobtz/zupmxode` repository, both checked
  2026-11-12). `Imports: zumir` only (D20). `LinkingTo: zumir` (zuexact left with D21), and `zuspecial`
  once its C API carries Lambert W
  ([zuspecial#61](https://github.com/pedrobtz/zuspecial/issues/61)) and
  it is on CRAN (D11). `Suggests: testthat, rxode2, mrgsolve, zupmxio,
  dazdiff`. One shared object; vendors `liblsoda` (MIT, D5). MIT.
- **Source:** the re-check in
  [`archive/reviews/pharmacometrics-review-2026-10-28.md`](https://github.com/pedrobtz/packages/blob/main/archive/reviews/pharmacometrics-review-2026-10-28.md);
  the rxode2 documentation on installation requirements, parallel solving
  and model syntax; the mrgsolve user guide and the nlmixr2 blog's
  comparison of the two; the EJCP atezolizumab comparison (2025) for the
  workload shape; the NONMEM Users Guide for PREDPP event semantics;
  `sdwfrost/liblsoda` (MIT, last pushed 2018) as the stiff solver;
  Hairer, Nørsett and Wanner for Dormand–Prince 5(4); Moler and Van Loan,
  "Nineteen dubious ways" (2003) for the matrix exponential; the standard
  closed-form solutions for one-, two- and three-compartment models with
  bolus, zero-order and first-order input; the `LinkingTo` headers of
  `zuexact` (API v1), `zurand` (API v3) and `zuspecial` (API v1) as read
  on 2026-10-28 (§3); the MIR repository at `master` commit `dc2c763`
  (2026-10-23) as read on 2026-11-12 for RFC 0073, whose §1 and §4 this
  RFC relies on for what MIR does. Nothing has been built.

## 1. What zupmxode is

`zupmxode` simulates population PK/PD models written in the rxode2 model
language, or the subset of it in §6, over NONMEM-style event records, on
every core, with no C compiler on the machine. The model is parsed in
R and compiled to a small intermediate representation of the package's
own (§4.1). `zumir` (RFC 0073) lowers that representation to native code
with MIR in well under a millisecond, or runs it in MIR's interpreter
where native code cannot be made (§4.7). Closed-form linear models are
solved exactly; linear ODE systems detected at compile time are solved
by a matrix exponential; the rest go to LSODA or Dormand–Prince, which
call the generated right-hand side as an ordinary C function. Subjects
are independent and solved in parallel with OpenMP.

The one-line statement that governs every decision below:

> **Any PK/PD model, every core, no compiler.**

Why that matters, measured against what exists:

- rxode2 and mrgsolve both generate C or C++ per model and compile it
  with the system compiler. A compiler with OpenMP support must be on
  the analysis machine, which validated GxP environments often forbid
  and locked-down laptops lack. FeRx needs a Rust toolchain and a first
  build of one to two hours.
- Each compile costs seconds. Automated model development runs hundreds
  of candidate models; the compile is a fixed tax per candidate. MIR's
  README measures about 250 µs to generate code for a small C program,
  a hundredth of `gcc -O2`'s time.
- A model that needs no compiler can be run inside a Shiny app on a
  server without a toolchain, in a container with only R, or by a
  reviewer with a stock R installation.

What it is not:

- **Not an estimator.** No likelihoods, no FOCEI, no SAEM. nlmixr2 and
  FeRx own estimation today. `zupmxode` is a simulation engine with a C
  API that an estimator can call; RFC 0072 `zupmxfit` is the family's
  gated attempt at one, and lives in its own package (D1, D12).
- **Not as fast as `gcc -O2` everywhere.** MIR's own figures put its
  generated code within about ten per cent of `gcc -O2`; its interpreter,
  the fallback tier, is six to fourteen times slower than that. §10
  measures both, and the README says which tier a platform gets.
- **Not a new model language.** The syntax is rxode2's, so every model in
  `nlmixr2lib` and every `nonmem2rx` output runs unchanged (D3).
- **Not a Julia replacement.** Pumas's JIT compiles the whole program;
  this package compiles only the model's blocks and keeps the solvers in
  C.

## 2. Scope

In v0.1.0:

- **The model language** (§6.1): the rxode2 subset.
- **Events** (§6.2): NONMEM data-set semantics for `ID`, `TIME`, `EVID`,
  `AMT`, `RATE`, `DUR`, `II`, `ADDL`, `SS`, `CMT`, covariates as columns
  with last-observation-carried-forward, and model-side `alag()`, `f()`,
  `rate()`, `dur()`.
- **Two execution tiers** (§4.7): native code generated by MIR, and
  MIR's interpreter on the same lowered code, chosen automatically with
  override.
- **Solvers** (§4.3): closed-form linear, matrix exponential, LSODA,
  Dormand–Prince 5(4), chosen at compile time with override.
- **Population simulation** (§5): etas drawn in R from omega with R's
  RNG, residual error applied in R, so every result is reproducible from
  `set.seed()`.
- **Parallelism**: OpenMP over subjects, thread count from an option.
- **Output**: a long data frame of subject, time, compartments and
  captured expressions.
- **A C API** header for siblings: compile a model, run one subject over
  an event table into caller-owned buffers, and evaluate one subject's
  predictions and residual variances at given parameters and etas
  (§4.6), the call an estimator makes millions of times.

Out of v0.1.0 (§14): forward sensitivities, dual-number code, vector
lanes, stochastic differential equations, delay differential equations,
`linCmt()` with transit compartments beyond what the closed forms
cover, time-varying matrix exponential, GPU, inductive linearisation,
a hand-written interpreter of the package's own representation.

## 3. Position in the `zu*` family

- **`zumir`** (RFC 0073): `LinkingTo`. It vendors MIR once for the
  family and exposes it through a versioned C API table. This package
  lowers its representation to MIR through that table and never vendors
  MIR itself (D14). `zumir` must be on CRAN before this package's
  submission (D11).
- **`zupmxio`** (RFC 0046): `Suggests`, both ways. `zupmxio` writes
  rxode2 text that `zupmxode` runs, and uses `zupmxode` as the engine
  that belongs to neither side of a translation. `zupmxio`'s corpus is
  this package's conformance suite.
- **`zupmxfit`** (RFC 0072): a consumer of `zupmxode_subject_eval()`
  (§4.6), from R in its first version and through `LinkingTo` later.
  The FOCE-I loop stays in that package; nothing about estimation enters
  this one.
- **`zuexact`** (RFC 0023; CRAN 0.1.0): `LinkingTo`. Its C API gives an
  opaque accumulator of `acc_size` bytes with `init`, `add1`, `addv`,
  `round` and `merge`, usable from OpenMP workers each with their own.
  Superposition over hundreds of additional doses and steady-state sums
  of exponentials are long sums where error accumulates; the closed-form
  solver keeps one accumulator per thread and rounds once (§4.3).
  Checked on 2026-10-28: nothing is missing.
- **`zuspecial`** (RFC 0027; not yet on CRAN): `LinkingTo` once two things
  land. Its C API already has `gammainc` and `gammaincc`, which the math
  table uses for Weibull absorption. It lacks Lambert W, which
  `sc_lambertw()` offers only at the R level; zuspecial#61 asks for a
  real-line `lambertw(x, k, code)` in the table. With it, one-compartment
  Michaelis–Menten elimination gets the Schnell–Mendoza closed form
  (§4.3). Until both land, that model classifies as general.
- **`zurand`** (not yet on CRAN): the v0.2.0 route for draws inside the
  parallel region (D7). Its C API v3 has `fold_int(key, i)` to derive a
  per-subject key and `fill_normal_at(key, start, n, ...)` for positioned
  draws, so a million subjects' etas are a pure function of one key,
  whatever the thread count. Checked on 2026-10-28: nothing is missing.
  v0.1.0 does not use it; draws stay in R with R's RNG.
- **`zuxad`** (RFC 0017): no role. Derivatives come from dual-number
  code the compiler emits (§14), which needs no AD library and none of
  `zuxad`'s AGPL terms.
- **`dazdiff`** (not yet on CRAN): `Suggests`. The three-engine
  conformance tests (§7, §9) compare prediction tables with
  `dazdiff::expect_table_equal()`. RFC 0046 Q5 covers its CRAN timing;
  the same answer applies here.
- **`zufuzz`**: development only. Fuzzes the parser, the lowering to MIR
  and the event scheduler (§9).
- **`zupeg`** (RFC 0039): a possible parser for the model language later;
  v0.1.0 parses R syntax with R's own parser because the language is R
  syntax (D3).
- **`nanowasm`**: no role in v0.1.0. A WebAssembly export of the
  representation is a possible later target, not an execution tier.

## 4. Architecture

### 4.1 Parser and compiler (R)

The model is an R expression or function body. R's parser reads it; a
walker in R checks it against the subset (§6.1), resolves names into
parameters, states, covariates, locals and captures, and classifies the
system:

- **closed-form** when the model is a `linCmt()` call with a recognised
  parameterisation and no ODE;
- **linear** when every `d/dt()` right-hand side is linear in the states
  with coefficients that depend only on parameters and covariates that
  are constant within an interval between events;
- **general** otherwise.

The compiler then emits the package's own intermediate representation
(the IR): straight-line register code over doubles, with typed slots for
states, parameters, covariates, locals and temporaries; instructions for
load, store, arithmetic, comparison, branch, and a fixed table of math
calls (`exp`, `log`, `sqrt`, `pow`, `abs`, `min`, `max`, `floor`,
`ceiling`, `sin`, `cos`, `tan`, `gamma`, `lgamma`, `logit`, `expit`,
`probit`, `probitInv`). Before emission it folds constants, merges
common subexpressions (`cl / v` once), and removes dead code.

Five code blocks come out:

| Block | Runs | Holds |
| --- | --- | --- |
| `pre` | at each event | individual parameters, lag, bioavailability, rate, duration |
| `ivl` | once per interval between events | expressions of `rhs` that depend on no state |
| `rhs` | at every solver call | the derivatives |
| `post` | at each output time | the captured expressions |
| `err` | at each observation | the residual variance (§4.6) |

`ivl` is the largest single optimisation for PK models: most of a
right-hand side's arithmetic is parameter algebra that changes only when
an event or a covariate does. The IR is an integer vector plus a
constant pool, cached on the model object and serialisable. It does not
depend on MIR (D13).

The format is specified where it is verified, in `src/core/ir.h`
(2026-11-15, Stage 2): a header with the counts, five block descriptors,
fixed-width instructions of an opcode and four operands, forward jumps
only, every register 0.0 on entry. The verifier accepts an IR only when
every access stays inside the arrays its header declares, and caps each
count: 10^6 instructions and 65,536 registers per block, 65,536 states,
parameters, work slots and outputs, 10^6 constants. A saved IR is
verified again before it is lowered.

**The math table's values** (decided 2026-11-15). At poles, infinities,
NaN and zero each function gives what R gives, since those are the values
a model written for rxode2 sees; `^` follows R's rules (a zero base has
no sign, a negative base to an infinite power is NaN). For finite
arguments each function is held to the true value, computed with Rmpfr
at 200 bits, within a bound per function, in units in the
last place; for the functions whose result passes through zero (`sin`,
`cos`, `tan`, `log`, `lgamma`, `logit`, `probit`) in units of
eps × max(|f(x)|, |x f′(x)|), the result's movement under one rounding
of x, since no library rounds correctly at a zero (half a unit in the last place
for `sqrt`, a few for the rest), not to R's result: R's `gamma()` and
`lgamma()` are further off than the C library's (by 228 units at
`gamma(100.5)`). The sign of a zero result is C's, as rxode2's compiled
code gives it. Comparisons follow C on NaN. Where a platform's C library is wrong,
the table does not take it: `sqrt` is compiled into the package (mingw-w64's
library rounds `sqrt(1 - 2^-53)` to 1), and `lgamma` is computed directly
near zero and from 1e10 (mingw-w64's overflows at both). `sin`, `cos`
and `tan` are held to the truth for |x| ≤ 1e5; beyond that they are the
platform library's, unchecked: mingw-w64's reduces 1e300 wrongly, and no
model reaches that range.

### 4.2 Event scheduler (C)

Per subject, the scheduler turns the data rows into a sorted timeline of
dose starts, infusion ends, lag-shifted doses, additional doses expanded
from `II` and `ADDL`, steady-state blocks, observation times and
covariate changes, following PREDPP's rules for `EVID` 0 to 4, `RATE`
including `-1` and `-2`, `SS` 1 and 2, and `CMT` with negative values for
turning compartments off. Steady state for closed-form and linear
systems is exact; for general systems it is reached by repeated dosing
until the pre-dose state changes by less than a tolerance, with a cap
and a warning. Lag and bioavailability come from `pre`, so they may
depend on individual parameters and covariates.

### 4.3 Solvers (C)

- **Closed-form**: one-, two- and three-compartment models with bolus,
  zero-order and first-order input, the state carried exactly between
  events (D21), exact steady state. Parameterisations: clearance and volume, micro-constants,
  and the `linCmt()` set. Oracle: rxode2's `linCmt()` and mrgsolve's
  `PKMODEL`. A compartment switched off is not supported here; the
  general solvers take such a model. One-compartment
  Michaelis–Menten elimination has the Schnell–Mendoza closed form
  through the principal branch of Lambert W; it ships when `zuspecial`
  exposes `lambertw` in C (zuspecial#61) and is on CRAN, and until then
  the classifier sends that model to the general solver. Oracle for it:
  LSODA at `rtol = 1e-10`.
- **Matrix exponential**: Padé with scaling and squaring on the system
  matrix, recomputed when a parameter or covariate changes, with the
  particular solution for constant infusion. For systems up to 20
  states; above that the model is solved by LSODA. Built in Stage 7
  (2026-11-16): Padé 13 (Higham 2005) of the augmented matrix
  [[A, c + r], [0, 0]], which carries the model's constant input c and
  the infusions r exactly; A and c come from the rhs block at y = 0 and
  at each unit vector, exact for a model affine in its states. Steady
  state is the fixed point of one interval, solved directly; SS = 2 adds
  the dose's own steady state without c, which the state already carries.
  A switched-off compartment takes a subject to LSODA.
- **LSODA**: `liblsoda` vendored (D5), one instance per thread, with
  the Jacobian by finite differences in v0.1.0. It calls `rhs` through a
  function pointer, as it would call any C function.
- **Dormand–Prince 5(4)**: written for the package, dense output for
  observation times, for non-stiff models where it beats LSODA.

How the general solvers run (2026-11-16, Stage 6): each interval between
timeline entries is integrated afresh, so a dose, an infusion's start or
end, a reset or a parameter change is never stepped across; a run of
observations with nothing between them is one Dormand–Prince integration
read by its dense output. A switched-off compartment is emptied and its
derivative held at 0 until a dose turns it on. Steady state is repeated
dosing from empty until the state before a dose moves by less than
`ss_tol` (1e-8 relative, with `atol`), at most `ss_cap` (1000) intervals,
with a warning past it. A failure, or a state that stops being finite,
leaves the subject NA from the time the solver stopped, with a warning;
LSODA can step over a finite-time singularity (dx/dt = x²) and report
success, which no general test separates from a large solution.

The default is the most specific solver the classification allows;
`pmx_ode_method()` overrides it within what the class allows, and
`pmx_ode_info(model)` says what was chosen and why.

### 4.4 Parallel runtime (C)

Subjects are independent. The runtime gives each thread one arena
(state, locals, solver workspace), splits subjects across threads, and
writes into a preallocated output matrix sized from the event table, so
no allocation happens inside the parallel region and no R object is
touched. The runtime checks for user interrupts between subjects on the
main thread. Without OpenMP the same code runs on one thread.

How it is built (Stage 8, 2026-11-16):

- One subject's solve is `src/core/solve.c`, which builds without R's
  headers, so the C tests compile it and nothing in it can reach the R
  API. `src/zpx_solve.c` does everything R-facing before and after the
  region: copies the records, counts each subject's observations (its
  rows of the output) and timeline entries (the arena it needs), sizes
  one arena per thread for the largest subject, and turns per-subject
  status structs into conditions afterwards, in subject order, so the
  first error reported does not depend on the thread count.
- Subjects go to threads dynamically, one at a time, in chunks of eight
  per thread; `R_CheckUserInterrupt()` runs between chunks, outside the
  region. Each subject writes only its own rows, so the output is the
  same bits at any thread count.
- On the interpreter tier each thread past the first gets its own MIR
  context, lowered from the IR before the region (zumir D10), owned by an
  external pointer so an interrupt cannot strand them. Native code is
  shared: it is re-entrant.
- `lgamma()` from the C library writes the sign of gamma to a global
  (`signgam`) on glibc and others; the math table's `lgamma` is the
  package's own (log |tgamma| below 15, Stirling's series from 15, the
  reflection formula below -15), checked against the Rmpfr reference.
- `pmx_ode_threads()` sets the count, 1 by default: CRAN allows two
  cores in checks, and a package that takes every core by default
  surprises a user already running in parallel.

### 4.5 R layer

`pmx_ode_model()` parses, compiles and lowers; `pmx_ode_sim()`
validates the event table, draws etas and samples residual error in R,
calls C once with every subject's parameters, and returns a data frame.
`pmx_ode_solve()` is the deterministic core: parameters per subject in,
predictions out.

How `pmx_ode_sim()` reads its arguments (decided 2026-11-16, Stage 9):

- **Etas** are the parameters `omega` names, or, with an unnamed
  `omega`, the model's parameters whose names start with `eta`, in the
  order the model uses them (nlmixr2's convention). A named vector is a
  diagonal `omega`. They reach the solve as covariate columns, one value
  per subject, so the C side is called once as for any solve.
- **Draws**, in this order: the etas as `MASS::mvrnorm(n, 0, omega)` makes
  them from the same seed (an eigendecomposition, so a positive
  semi-definite `omega` works), one row per subject; then, for each
  `sigma` entry in turn, one standard normal per observation record.
  The count depends on the table alone, never on the solve, so a failed
  subject leaves the other draws where they were.
- **Error models** take rxode2's spellings, as the language does:
  `add(sd)`, `prop(sd)`, `lnorm(sd)` for the exponential, and
  `add(a) + prop(b)` for the combined model, whose variances add
  (nlmixr2's `combined2`, its default). Each takes a standard deviation.
- **Output:** the solve's columns unchanged, then each subject's etas,
  then each `sigma` entry with residual error as `<name>_obs`, so the
  prediction and the observation sit side by side.
- **Subjects:** a table of one subject is copied `n` times, numbered 1
  to `n`; a table of several is simulated once per subject.

### 4.6 Per-subject evaluation (C API)

An estimator does not want a data frame back. It wants, for one subject
at one parameter vector, the predictions and the residual variances at
that subject's observation rows, written into buffers it owns, as fast
as the engine can produce them. The header carries one entry point for
that:

```c
int zupmxode_subject_eval(const zupmxode_model *m,
                          const zupmxode_subject *s,
                          const double *theta, const double *eta,
                          double *pred, double *var,
                          const zupmxode_eval_opts *opts,
                          zupmxode_arena *arena);
```

- `m` is a compiled model and `s` a subject whose timeline the scheduler
  has already built, so neither is rebuilt per call.
- `pred` and `var` have one entry per observation row of `s`. `var` is
  the residual variance under the error model and Sigma carried in
  `opts`, one of the four forms `pmx_ode_sim()` applies (Q2): additive,
  proportional, combined, exponential, evaluated by the `err` block. The
  caller forms the likelihood without re-implementing the error models.
- `opts` carries the solver choice and tolerances. An estimator taking
  finite differences through LSODA asks for tighter tolerances than a
  simulation does (RFC 0072 §5), and must be able to say so per call.
- `arena` is caller-owned workspace from `zupmxode_arena_new()`, one per
  thread, so the call allocates nothing and touches no R object.
- The call is re-entrant and keeps no global state. Two threads may
  evaluate different subjects of the same model at once.
- It returns a status code (§8) and never calls `error()`.

The struct reserves room for sensitivity outputs (§14) without an ABI
break: `opts` has a version field, and a later version adds
`dpred_deta` and `dvar_deta` buffers behind it. v0.1.0 fills neither.

The same call is exported to R as `pmx_ode_subject_eval(model, subject,
theta, eta, sigma, tol)`, on a subject prepared once by
`pmx_ode_subjects(model, data)`. An estimator prototyped in R pays one
R call per subject per evaluation, not per solver step, which is the
trade RFC 0072 makes in its first version.

This costs little now and a rewrite later: the solvers already take
their workspace from a per-thread arena (§4.4), so the entry point is a
thin layer over the code `pmx_ode_solve()` runs.

### 4.7 Lowering and execution tiers (C, through `zumir`)

`pmx_ode_model()` passes the IR to C, which lowers each block to one MIR
function with a fixed C signature, for example

```c
void rhs(double t, const double *y, double *dydt,
         const double *par, double *work);
```

and asks `zumir` for code:

- **Native tier**: MIR's generator turns the module into machine code.
  The blocks are then ordinary function pointers. The solvers, the
  scheduler and every OpenMP thread call them directly; there is no
  boundary to cross per call.
- **Interpreter tier**: MIR's interpreter runs the same lowered module.
  It is used where `zumir` reports that native code cannot be made: an
  architecture MIR does not target, or a system that refuses memory that
  is both writable and executable. The blocks are reached through a thin
  C trampoline with the same signature, so the solvers do not know which
  tier is running.

Lowering and generation happen once per model, on the main thread,
before any parallel region, in a MIR context owned by the model object
and released by its finalizer. MIR contexts are not shared between
threads; the generated code is, because it is plain re-entrant code
over the arena it is given. `pmx_ode_info()` reports the tier, the time
lowering took and the size of the generated code. `pmx_ode_tier(m,
"interp")` forces the interpreter, which the test suite uses to compare
the two (§9).

## 5. Public R API (v0.1.0 shape)

```r
m <- pmx_ode_model({
  ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
  d/dt(depot)   <- -ka * depot
  d/dt(central) <-  ka * depot - cl / v * central
  cp <- central / v
})
pmx_ode_info(m)    # states, parameters, solver, tier, lowering time

ev <- pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
                      time = seq(0, 168, by = 1))
out <- pmx_ode_solve(m, ev, params = c(tka = 0.5, tcl = 1, tv = 3,
                                        eta.ka = 0, eta.cl = 0,
                                        eta.v = 0))

set.seed(1)
sim <- pmx_ode_sim(m, ev, theta = c(tka = 0.5, tcl = 1, tv = 3),
                    omega = diag(c(0.1, 0.1, 0.1)),
                    sigma = list(cp = prop(0.1)), n = 1000)
pmx_ode_threads(8)
pmx_ode_solve(m, nonmem_data)   # a NONMEM-style data frame as events
pmx_ode_method(m, "lsoda")      # override the solver
pmx_ode_tier(m, "interp")       # override the tier
```

## 6. Semantics

### 6.1 The language

The subset of rxode2 syntax: assignments with `<-` or `=`, `d/dt(x)`,
`x(0)`, `alag(x)`, `f(x)`, `rate(x)`, `dur(x)`, `if`/`else`, the math
table in §4.1, comparison and logical operators, `linCmt()`, and
`cmt()`. Anything else is a compile error naming the construct and the
rxode2 feature it belongs to. Models that compile under `zupmxode`
compile unchanged under rxode2; the reverse is not promised.

Details decided when the compiler was built (2026-11-15, Stage 3):

- **Spellings.** rxode2's `lag(x)` and `F(x)` are accepted for `alag(x)`
  and `f(x)`, and its bounded `logit(x, low, high)` and
  `expit(x, low, high)` as the one-argument functions rescaled; `**` is
  `^`, as R parses it; `min()` and `max()` take two or more arguments;
  `t` and `time` are the time, `pi`, `Inf`, `NaN`, `TRUE` and `FALSE`
  constants. `ifelse()` is not accepted: `if`/`else` is the subset's form.
- **Names.** Compartments are named by `cmt()`, then by `d/dt()` in order
  of appearance, with `linCmt()`'s first. Every name assigned is a local
  and an output; every name read and never assigned is a parameter, in
  order of first use. A local read on a path that has not assigned it is
  a compile error, rather than rxode2's value from an earlier call.
- **`linCmt()`** becomes rxode2's compartments (`depot`, `central`,
  `peripheral1`, `peripheral2`) and their equations, from the
  parameterisations §4.3 names; the model is closed-form only when it has
  no other derivative, and the general solvers solve the same equations.
- **`pre`** holds lags, bioavailabilities, rates, durations and initial
  amounts, which may depend on parameters and covariates but not on a
  compartment or on time.
- **Hoisting** takes a top-level assignment into `ivl` when it is the
  only assignment to its name and reads no compartment and no time,
  directly or through a local.
- **Folding** evaluates a constant subexpression with the math table
  itself, so a folded value is the double the generated code computes.

### 6.2 Events

PREDPP semantics as documented in the NONMEM Users Guide, with rxode2's
documented behaviour as tie-breaker where PREDPP leaves a choice. Where
`zupmxode`, rxode2 and mrgsolve disagree on a corner (a lagged dose that
straddles an infusion end, `SS=2` with a changed interval), the test
suite records all three and the documentation says which one
`zupmxode` follows and why.

How the scheduler reads the records (decided 2026-11-16, Stage 4;
`src/core/sched.h`). The three-engine record of Stage 5 says where
rxode2 or mrgsolve read a case otherwise.

- **Order at one time:** infusion ends, then doses the scheduler made
  (additional doses, doses moved by their lag), then the records in data
  order. An observation written before a dose at the same time sees the
  amount before it.
- **Bioavailability and infusions:** a given `RATE` keeps its rate and F
  changes the duration (PREDPP); a given `DUR`, or the model's `dur()`
  under `RATE = -2`, keeps its duration and F changes the rate.
  `RATE = -1` takes the model's `rate()`. `DUR` with a `RATE` is an error.
- **Additional doses** are at `TIME + k * II`, computed so, not summed,
  and take the lag, bioavailability, rate and duration of the record
  that made them.
- **Steady state** is an entry before the dose it belongs to: with
  `SS = 1` the compartments are reset and set to the steady state just
  before a dose of that interval, with `SS = 2` that steady state is
  added to what is there, and the dose then follows. A steady-state
  infusion with `II = 0` runs on with no end.
- **`EVID 2`** marks a time; a negative `CMT` switches that compartment
  off, a positive one on. `EVID 3` resets every compartment; `EVID 4`
  resets and doses at once.
- **`MDV`:** with no `EVID` column, a record with an `AMT` above 0 is a
  dose, one with `MDV` 1 marks a time (`EVID 2`), and the rest are
  observations (decided 2026-11-16, Stage 9, as NONMEM reads them).
  Given `EVID`, `MDV` is not read. Other columns, `DV` among them, are
  ignored unless named like a parameter.
- **Covariates** are the table's other columns named like parameters,
  each record carrying its value, so a value holds from its record until
  the next (last observation carried forward); the model's lag and
  bioavailability are taken at the dose's record.
- **Limits:** one subject's timeline holds at most 10^6 entries.

### 6.3 Numerics

Default tolerances match rxode2 (`rtol = 1e-6`, `atol = 1e-8`) so
comparisons are like for like. Closed-form and matrix-exponential
results agree with the general solver on the same model to within those
tolerances; this is a test, not an assumption.

The two tiers are expected to give identical doubles, since the IR has
no operation whose rounding is left to the code generator. That is
tested on every model (§9), and a difference is a bug unless it is
traced to a documented contraction such as a fused multiply-add on one
architecture, in which case the tolerance and the cause are recorded.

## 7. The oracle

rxode2 and mrgsolve, in `Suggests`, run in CI when installed. Every
model in the conformance suite (§9) is solved by all three on the same
events and the results compared at the §6.3 tolerances. The suite is
`zupmxio`'s corpus and the `nlmixr2lib` model library plus mrgsolve's
`modlib()`. Where the two oracles disagree with each other, the case is
recorded and `zupmxode` follows PREDPP. Each tier is compared against
the oracles, and the tiers against each other.

How the suite is held, decided 2026-11-15 when it was built (roadmap
Stage 1):

- **Recorded, then compared.** `tools/make-corpus.R` solves the corpus
  with both oracles at rtol 1e-10, atol 1e-12, and writes their answers
  to `tests/testthat/fixtures/oracle/` with the engine versions. Tests
  compare with those fixtures, so the `nosuggests` leg still checks
  answers; rxode2 is in `Suggests` only so a test can show the fixtures
  are still what rxode2 gives.
- **Paired models.** mrgsolve's `modlib()` and the §6.2 event features
  enter as models written in both languages (`tools/corpus/`): pk1, pk2
  and pk3 as ODEs and closed forms, Michaelis–Menten, an indirect
  response, an effect compartment, lag and bioavailability, modelled rate
  and duration, a time-varying covariate. Running `modlib()`'s own text
  needs a translator, which is `zupmxio`'s; until then these ports are
  the mrgsolve side.
- **nlmixr2lib** models are rxode2-only cases: residual-error lines
  removed, initial thetas, etas at 0, reference covariates. Those needing
  an input the generator has no value for are listed as skipped in the
  manifest, not dropped silently.
- **Settings made explicit**, each of which caused a disagreement before
  it was set: mrgsolve reads covariates last-observation-carried-forward
  (`nocb = FALSE`); both engines iterate steady state to the solve's own
  tolerance; doses sort before observations at the same time. Remaining
  disagreements are listed in `disagreements.csv`, each with a PREDPP
  note in `disagreements.md`; there are none on the 2026-11-15 corpus.
- **Tables** are compared by `tests/testthat/helper-tables.R`: a value
  matches when within atol or within rtol of the larger magnitude, the
  rule `dazdiff` uses, through `dazdiff` when installed (D16).

## 8. Errors

- Compile errors name the construct, its position in the model and the
  nearest supported alternative.
- Lowering errors, which are bugs in this package or in MIR, name the
  block and the IR instruction, and the model falls back to the
  interpreter tier with a warning if only generation failed.
- Event-table errors name the row and the rule (`AMT` without `EVID`,
  `ADDL` without `II`, unknown compartment).
- Solver failures return `NA` for the subject with a warning naming the
  subject and the time, never an error that loses the other subjects.
- The C API returns status codes; it never calls `error()`. MIR errors
  reach this package as `zumir` status codes, never as a process exit
  (RFC 0073 D4).

## 9. Testing

- Closed-form against general solver on every linear model in the suite,
  and the Michaelis–Menten closed form against LSODA at tight tolerance.
- Three-engine comparisons and the closed-form checks use
  `dazdiff::expect_table_equal()` on the prediction tables.
- Native tier against interpreter tier on every model in the suite:
  identical output, or the recorded exception of §6.3.
- Superposition over 10³ additional doses: the `zuexact` sum against a
  long-double reference, and against the naive sum to show the gap.
- Matrix exponential against the general solver on linear ODE models.
- Every event-semantics rule in §6.2 as one test with a hand-computed
  expectation and the three-engine record.
- Thread count 1 against thread count 8: bitwise identical output, on
  each tier.
- The `ivl` hoisting against a build with hoisting switched off:
  identical output.
- Fuzz targets over the parser, the lowering and the scheduler, the last
  two built standalone without R.
- Allocation test: the parallel region allocates nothing on the R heap.
- `pmx_ode_subject_eval()` against `pmx_ode_solve()` on every model in
  the suite: identical predictions. The same call from eight threads on
  different subjects of one model: identical to the serial results.

## 10. Performance targets

Measured by `tools/run-benchmarks` against rxode2 and mrgsolve on the
same machine, recorded in `.agents/benchmarks.md`, with the oracles' baseline recorded
before any solver is built.

| Workload | Target |
| --- | --- |
| compile and lower a two-compartment model | under 20 ms; rxode2, mrgsolve: seconds |
| closed-form, 10⁶ subjects, 8 threads | matches rxode2 `linCmt()` run time |
| 2-cpt ODE, 10³ subjects, 10³ times, 1 thread, native | within 1.5× of rxode2 |
| same, interpreter tier | recorded |
| same, 8 threads against rxode2 8 threads, native | within 1.5× |
| EJCP seven-cycle workload, synthetic | above 5,000 obs/s (NONMEM's figure) |
| 200 candidate models, compiled and run once | faster than rxode2 end to end |

The second row excludes rxode2's compile; the first and last rows charge
it. The last row is the workload automated model development pays. The
interpreter row has no target because it exists for platforms that have
no alternative; the README states the figure.

## 11. Roadmap

The stages, their order and exit criteria are in
[roadmap.md](roadmap.md). Until `zumir` is on CRAN this package builds
against its GitHub version through `Remotes:`.

### 11.1 If the MIR spike fails

`zumir`'s Stage 0 (RFC 0073 §10) is the first analysis, and it looks at
MIR alone. Its outcome sets this package's back-end before the lowering
is built:

| Spike outcome | Back-end for v0.1.0 |
| --- | --- |
| native code builds and runs on every r-actions leg | both tiers as §4.7 says |
| native code fails on Windows only | Windows gets the interpreter tier while the patches go upstream; the README says so |
| MIR's interpreter is unusable as a fallback, or native code fails more widely | a second analysis: a small interpreter of the IR written for this package, and LuaJIT through `luajr` as its comparison |

Because the IR is the package's own (D13), the third row replaces the
lowering and nothing else.

**Outcome, 2026-11-13: the first row.** MIR builds and runs on every
r-actions leg
([mirspike RESULTS.md](https://github.com/pedrobtz/mirspike/blob/main/RESULTS.md),
`zumir` Stage 0). Two `zumir` decisions taken since bear on §4.7: its
D10 (the interpreter tier needs one MIR context per thread, compiled
before the parallel region; native code is still compiled once and
shared) and its Q5 (where W+X memory is refused, `zumir_tier()` is
`"none"` and neither tier runs). Q8 asks what this package does then.

## 12. Acceptance criteria for v0.1.0

- Every behaviour in §6 and §8 has a passing test.
- Every model in the conformance suite matches rxode2 and mrgsolve at
  §6.3 tolerances, or the disagreement is recorded with the PREDPP
  reading that `zupmxode` follows.
- The two tiers agree on every model, as §6.3 defines.
- Output is identical across thread counts.
- The parallel region allocates nothing on the R heap.
- Every §10 row is measured and recorded, met or not; the two rows the
  package exists for (compile time, candidate loop) are met.
- `R CMD check --as-cran` is clean on all r-actions legs, including the
  legs without OpenMP; no `Imports` but `zumir` (D20). The README names the tier each leg
  runs.

## 13. Decisions

Marked **(M)** when the maintainer decided it and **(C)** when this RFC
decides it and the maintainer may overturn it.

- **D1 (C).** The package is a simulation engine, not an estimator.
  Estimation, if the family builds it, is RFC 0072's package, on top of
  this one.
- **D2 (M), 2026-11-12, provisional on `zumir`'s spike.** Models compile
  to native code with MIR, through `zumir`, and fall back to MIR's
  interpreter. This replaces the original D2, a bytecode virtual machine
  written for the package. The trade the original D2 accepted, two to
  three times slower right-hand sides for no toolchain, mostly goes
  away: MIR needs no system compiler and generates code close to
  `gcc -O2` in microseconds, and the generated functions are called by
  the C solvers with no boundary to cross, which LuaJIT could not offer.
  The costs are a dependency of about 23,000 lines for MIR's core and
  generator, plus its interpreter, with one principal maintainer;
  memory that is both writable and executable on the native tier; and
  Windows under Rtools, which MIR's own CI does not test.
  §11.1 says what replaces it if the spike shows those costs are too
  high.
- **D3 (C).** The language is a subset of rxode2's and is parsed by R's
  parser. Existing models run unchanged; no second syntax to learn.
- **D4 (C).** Event semantics follow PREDPP, with rxode2 as tie-breaker.
- **D5 (C).** LSODA is `liblsoda` vendored (MIT, about 2,000 lines),
  made re-entrant per thread. Dormand–Prince is written for the package.
  CVODE is not vendored: SUNDIALS is large and BSD-3, and its gain over
  LSODA on models of this size is not shown.
- **D6 (C).** Solver selection is automatic from the classification,
  with override. Tier selection is automatic from what `zumir` reports,
  with override.
- **D7 (C).** Random draws happen in R with R's RNG in v0.1.0; C is
  deterministic. A simulation is reproducible from `set.seed()` on every
  platform and thread count. v0.2.0 may add draws inside the parallel
  region through `zurand`, keyed per subject with `fold_int()` and
  positioned with `fill_normal_at()`, reproducible from one key at any
  thread count. The R route stays the default so `set.seed()` keeps
  working; the `zurand` route is for simulations too large to draw in R.
- **D8 (C).** OpenMP over subjects only; no parallelism within a
  subject.
- **D9 (C).** A C API header ships at v0.1.0 so `zupmxio` and estimators
  can drive the engine without R objects.
- **D10 (M), 2026-11-12.** The name is `zupmxode`: `zupmx`, the
  pharmacometrics sub-family mark set by `zupmxio`'s D24, then the job.
  It replaces `zufaode` (`fa` meant pharmacology, as `pmx` now does) and
  the original D10's `zu` + `pk` + role. A general name, `zuode`, was
  weighed and set aside: RFC 0050 `zussa`, RFC 0053 `zufmi` and the SBML
  and PBPK sketch in `base-layer-candidates.md` will use the engine too,
  but PK/PD is its first use and its design (PREDPP events, closed-form
  compartment models) is shaped by it, and D24 is kept without an
  exception. Those packages depend on a pharmacometrics-named engine,
  which is a cost of naming only. R functions are written `zupmxode_*`
  in this RFC; Q6 asks whether they take D24's `pmx_` prefix instead.
- **D11 (M).** Siblings that already offer what the engine needs are
  linked, not re-implemented: `zumir` for code generation and `zuexact`
  for exact sums at v0.1.0, `zuspecial` for Lambert W and incomplete
  gamma once zuspecial#61 and its CRAN release land, `zurand` for in-C
  draws at v0.2.0. A sibling not on CRAN cannot sit in `LinkingTo` of a
  CRAN submission, so each is gated on its release. The engine degrades
  without the optional ones: plain sums, Michaelis–Menten as a general
  model, draws in R. `zumir` is not optional; its release gates this
  package's submission.
- **D12 (C), 2026-11-12.** `zupmxode_subject_eval()` (§4.6) ships at
  v0.1.0, re-entrant, allocation-free and with an extensible options
  struct, even though no estimator exists yet. RFC 0072 depends on it,
  and adding it after the solvers settle would mean reworking how they
  take workspace.
- **D13 (C), 2026-11-12.** The IR is the package's own, and MIR is a
  target of it, not its format. The compiler's optimisations, the
  serialised model, the fuzz targets and the tests are written against
  the IR. MIR's own IR can change without touching any of them, and if
  MIR has to be replaced (§11.1), only the lowering is rewritten.
- **D14 (C), 2026-11-12.** MIR is reached through `zumir`'s C API, never
  vendored here. One copy of MIR for the family means one set of
  platform patches and one place where MIR's process-exiting error
  handler is replaced.

- **D15 (M), 2026-11-15.** Where `zumir_tier()` is `"none"` (W+X memory
  refused, `zumir` Q5), v0.1.0 reports the system as unsupported with a
  classed error naming `zumir_tier()`. An IR interpreter for that case is
  §14 work; D13 keeps the IR ready for it. Answers Q8.
- **D16 (M), 2026-11-15.** Tests compare prediction tables with a helper
  in `tests/testthat/`, through `dazdiff::expect_table_equal()` when it
  is installed; `dazdiff` leaves `Suggests` before submission unless it
  is on CRAN by then, so its release does not gate this package. Answers
  Q9.

- **D17 (M), 2026-11-15.** The R functions take the sub-family prefix
  of `zupmxio`'s D24: `pmx_ode_model()`, `pmx_ode_solve()`,
  `pmx_ode_sim()` and so on, beside RFC 0072's `pmx_fit()`. The C API
  and its types keep `zupmxode_`. Answers Q6.
- **D18 (M), 2026-11-15.** A saved model serialises only the IR; code is
  regenerated on load, which takes under a millisecond, so nothing
  specific to a MIR version is ever stored. Answers Q7.
- **D19 (M), 2026-11-15.** `pmx_ode_sim()` ships the four residual error
  models (additive, proportional, combined, exponential); any other is
  a function the user applies to the predictions, as the documentation
  shows. Answers Q2.

- **D20 (C), 2026-11-15.** `zumir` is in `Imports` as well as
  `LinkingTo`, the only package in `Imports`. Its C API table is fetched
  with `R_GetCCallable()`, which needs `zumir`'s namespace loaded, and
  `Imports` is what loads it; `zumir`'s own stand-in consumer does the
  same. This adds no dependency: §3 and D11 already make `zumir`
  required. Found in Stage 2, when the first lowering ran.

- **D21 (M), 2026-11-16.** The closed forms carry the compartments'
  state exactly from one event to the next, through the system's
  eigenvalues (closed-form roots for one to three compartments), instead
  of summing every dose's response. Superposition is valid only while a
  subject's parameters stay constant, which a time-varying covariate, a
  reset or a switched-off compartment breaks. With no long sums left,
  `zuexact` leaves `LinkingTo`. Steady state is the fixed point of one
  dosing interval. Where two eigenvalues nearly meet (ka at an
  elimination rate), a matrix exponential takes over. Tested against
  200-bit superposition where superposition holds.

## 14. Later versions, in the order they are expected

- **Dual-number code.** The compiler emits each block a second time over
  dual numbers, a value and *n* derivatives, and `zumir` generates it
  like any other block. That gives forward sensitivities of states with
  respect to parameters and etas, without an AD library and through
  `if`/`else` as written. They fill the `dpred_deta` and `dvar_deta`
  slots §4.6 reserves. RFC 0072's Gate C decides whether they are
  needed: if finite differences through LSODA give stable FOCE-I
  estimates on ODE models, they wait; if not, they are the next version
  of this package.
- **Analytic Jacobian** for LSODA from the same dual-number code.
- **Vector lanes.** Blocks lowered to operate on *K* parameter sets at
  once: closed-form models across subjects, and the `2·nη + 1`
  finite-difference perturbations `zupmxfit` evaluates per subject, in
  lockstep on one solver step sequence, which also removes most of the
  noise adaptive steps put into those differences.
- **Draws in C through `zurand`** (D7), and the **Michaelis–Menten
  closed form** if zuspecial#61 lands after v0.1.0 (D11).
- **Adjoint sensitivities for large models** through `r2sundials`, the
  CRAN binding of CVODES with forward and adjoint modes, as an optional
  backend behind `method = "cvodes"` rather than a vendored SUNDIALS.
  The adjoint costs about one extra solve regardless of parameter
  count, which pays when parameters far outnumber states, as QSP and
  systems-biology models do and population PK models with five to ten
  individual parameters do not. Forward sensitivities above stay the
  NLME route; D5 stands. `sundialr` is the other CRAN binding and has
  forward mode only. A second conformance corpus for this backend is
  the PEtab benchmark collection, once `zupmxio` reads SBML
  ([zupmxio#17](https://github.com/pedrobtz/zupmxio/issues/17)).
- **Stochastic and delay differential equations.**
- **Transit-compartment closed forms** and time-varying linear systems.
- **Inductive linearisation** if a model class shows the need.
- **A `zupmxio` writer target** so a model held in the representation
  runs without passing through rxode2 text.
- **Streaming output** for simulations larger than memory.

## 15. Open questions

- **Q1.** Should the IR be exposed and documented, so that another
  package (an SBML or FMI engine, RFC 0053) could emit it directly, or
  stay internal until the instruction set settles? This RFC keeps it
  internal for v0.1.0.
- **Q2. Decided by D19, 2026-11-15.** Should `pmx_ode_sim()` support residual error models beyond
  additive, proportional, combined and exponential at v0.1.0, or defer
  to the user applying a function to predictions? This RFC ships the
  four and documents the function route.
- **Q3.** Is the EJCP workload reproducible from the paper's description
  closely enough to be a benchmark row, or should the row be a
  synthetic equivalent? This RFC uses a synthetic equivalent of the
  stated size and says so.
- **Q4.** Is v0.1.0 a CRAN submission or a GitHub release? The board's
  `CRAN intent` for this card decides; this RFC assumes `Planned`.
- **Q5.** Should v0.1.0 wait for zuspecial#61 and `zuspecial` on CRAN so
  the Michaelis–Menten closed form ships at release? This RFC does not
  wait; the general solver is the fallback and §14 carries the form.
- **Q6. Decided by D17, 2026-11-15.** `zupmxio`'s D24 gives the sub-family the `pmx_` function
  prefix, and RFC 0072 follows it with `pmx_fit()`. Should this
  package's R functions be `pmx_ode_model()`, `pmx_ode_solve()` and so
  on, rather than `zupmxode_*`? This RFC keeps `zupmxode_*` so the
  engine's API reads as its own; the C API keeps `zupmxode_` either
  way.
- **Q7. Decided by D18, 2026-11-15.** Should a compiled model's generated code be cached across
  sessions, with MIR's binary format, or regenerated on load? Generation
  is fast enough that this RFC regenerates, and serialises only the IR.
- **Q8. Decided by D15, 2026-11-15.** `zumir`'s Q5 (decided
  2026-11-13) leaves no tier where W+X memory is refused, and says
  `zupmxode` brings its own IR interpreter for that case; §2 puts that interpreter out of v0.1.0 and §4.7 assumed
  MIR's interpreter covered it. Does v0.1.0 ship an IR interpreter, or
  report such a system as unsupported with a classed error naming
  `zumir_tier()`? Proposed: the classed error at v0.1.0, the
  interpreter in §14, since D13 keeps the IR ready for it.
- **Q9. Decided by D16, 2026-11-15.** `dazdiff`, which §9 uses for table comparisons, is not on CRAN
  (checked 2026-11-15), and a CRAN submission cannot suggest a package
  CRAN lacks. Wait for its release, or compare tables with a helper in
  `tests/testthat/` until then? Proposed: the helper, with `dazdiff`
  used when installed, so RFC 0046 Q5's answer does not gate this
  package.
