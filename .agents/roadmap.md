# zupmxode: roadmap

Status: draft 2026-11-15. Implements [design.md](design.md); replaces the
ten-row plan that was its §11, re-cut into one-pull-request stages with
the lowering, the riskiest boundary, before the parser.
Target: v0.1.0, the first CRAN release.
Tracking: [v0.1.0 (#1)](https://github.com/pedrobtz/zupmxode/issues/1), with
one sub-issue per stage;
[CRAN-submission (#18)](https://github.com/pedrobtz/zupmxode/issues/18).

Precondition met: `zumir`'s MIR spike passed on 2026-11-13 with the first
outcome of [§11.1](design.md#111-if-the-mir-spike-fails), so the back-end
is the one §4.7 describes.

## Stages

| Stage | Title | CI | Issue | Pull request | Done |
| --- | --- | --- | --- | --- | --- |
| 0 | [Package identity](#stage-0--package-identity) | quick | [#2](https://github.com/pedrobtz/zupmxode/issues/2) | | 2026-11-15 |
| 1 | [Oracle corpus and baseline](#stage-1--oracle-corpus-and-baseline) | quick | [#3](https://github.com/pedrobtz/zupmxode/issues/3) | | |
| 2 | [IR and lowering through zumir](#stage-2--ir-and-lowering-through-zumir) | full | [#4](https://github.com/pedrobtz/zupmxode/issues/4) | | |
| 3 | [Parser and compiler](#stage-3--parser-and-compiler) | quick | [#5](https://github.com/pedrobtz/zupmxode/issues/5) | | |
| 4 | [Event scheduler](#stage-4--event-scheduler) | quick | [#6](https://github.com/pedrobtz/zupmxode/issues/6) | | |
| 5 | [Closed-form solver](#stage-5--closed-form-solver) | quick | [#7](https://github.com/pedrobtz/zupmxode/issues/7) | | |
| 6 | [LSODA and Dormand–Prince](#stage-6--lsoda-and-dormandprince) | full | [#8](https://github.com/pedrobtz/zupmxode/issues/8) | | |
| 7 | [Matrix exponential and solver selection](#stage-7--matrix-exponential-and-solver-selection) | quick | [#9](https://github.com/pedrobtz/zupmxode/issues/9) | | |
| 8 | [Parallel runtime](#stage-8--parallel-runtime) | full | [#10](https://github.com/pedrobtz/zupmxode/issues/10) | | |
| 9 | [Population simulation](#stage-9--population-simulation) | quick | [#11](https://github.com/pedrobtz/zupmxode/issues/11) | | |
| 10 | [Per-subject evaluation and the C API](#stage-10--per-subject-evaluation-and-the-c-api) | quick | [#12](https://github.com/pedrobtz/zupmxode/issues/12) | | |
| 11 | [Benchmarks and README](#stage-11--benchmarks-and-readme) | quick | [#13](https://github.com/pedrobtz/zupmxode/issues/13) | | |
| 12 | [Assurance](#stage-12--assurance) | quick | [#14](https://github.com/pedrobtz/zupmxode/issues/14) | | |
| 13 | [Full CI](#stage-13--full-ci) | full | [#15](https://github.com/pedrobtz/zupmxode/issues/15) | | |
| 14 | [v0.1.0 readiness](#stage-14--v010-readiness) | full | [#16](https://github.com/pedrobtz/zupmxode/issues/16) | | |
| 15 | [Review](#stage-15--review) | full | [#17](https://github.com/pedrobtz/zupmxode/issues/17) | | |

The oracle every stage compares with, per boundary, is in
[§7](design.md#7-the-oracle): rxode2 and mrgsolve, live when installed,
and their recorded outputs (Stage 1) when not, so the `nosuggests` leg
still checks answers rather than skipping them.

## Stage 0 — Package identity

**Issue.** [#2](https://github.com/pedrobtz/zupmxode/issues/2), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** The package has its real identity, the CI every later stage
runs under, and an assurance file that declares what it will owe.

**Work.**

- `DESCRIPTION`: real Title, Description and `Authors@R` (`aut`, `cre`,
  `cph`); `LinkingTo: zumir, zuexact`; `Remotes: pedrobtz/zumir`; no
  `Imports`.
- `R-CMD-check.yaml` with the Windows R-devel runner row, as `zumir`'s.
- `.agents/assurance.md` from the template with kinds `compiled`,
  `untrusted-input` (a serialised model is bytes the package reads
  back), `spec` (PREDPP, LSODA, Dormand–Prince, closed forms) and
  `lifecycle` (the compiled model and the prepared subjects); every row
  `gap #<issue>` or `n/a` until its stage; the `assurance.yml` caller
  with `ungated-jobs: fail`.
- `tools/ci-run-tests.R`: shuffled suite, seed printed, coverage on.
- Labels `full-ci`, `ci:minimal`, `ci:docs`, `stage`.
- `src/core/` for C that builds without R headers, with a
  `tests/c/` Makefile that compiles it standalone; `src/zupmxode.c`
  stays the only file including R's headers until Stage 2 adds the
  glue.

**Tests.**

- Standalone C build: `make -C tests/c` builds `src/core/` with no R
  headers on the path, in the quick R CMD check job.
- Suite order: two runs with different printed seeds both pass, in CI.

**Exit criteria.**

- `R CMD check` passes on the quick r-actions profile, every check under
  five minutes.
- `assurance.yml` green with every gap tracked by an issue.
- The standalone C build runs in CI and fails when `src/core/` includes
  `R.h` (seen once).

**Design.** [§2](design.md#2-scope), [§3](design.md#3-position-in-the-zu-family).

## Stage 1 — Oracle corpus and baseline

**Issue.** [#3](https://github.com/pedrobtz/zupmxode/issues/3), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** The answers every later stage is checked against exist before
the code that must reproduce them.

**Work.**

- `tools/make-corpus.R`: the conformance models (the `nlmixr2lib`
  library, mrgsolve's `modlib()`, and the §6.2 event cases) with event
  tables, solved by rxode2 and mrgsolve; outputs written to
  `tests/testthat/fixtures/oracle/` with the engine versions, the
  tolerances and the date. The script is the provenance; the fixtures
  are regenerated, never edited.
- `helper-tables.R`: compare prediction tables at the §6.3 tolerances,
  through `dazdiff::expect_table_equal()` when installed (D16).
- Where rxode2 and mrgsolve disagree, the case is listed in
  `fixtures/oracle/disagreements.md` with the PREDPP reading.
- `tools/run-benchmarks` and `.agents/benchmarks.md` with the oracle
  rows of §10 recorded on one machine.

**Tests.**

- The helper: a table differing by more than `rtol` in one cell fails,
  one within it passes, in `test-helper-tables.R`.
- Fixture freshness: when rxode2 is installed, re-solving three corpus
  models reproduces their fixtures, in `test-oracle-fixtures.R`.

**Exit criteria.**

- The corpus fixtures and the disagreement list are committed, generated
  by the script.
- The oracle rows of §10 are recorded in `.agents/benchmarks.md`.
- Both tests exist and were seen to fail on a perturbed fixture.

**Design.** [§7](design.md#7-the-oracle),
[§10](design.md#10-performance-targets).

## Stage 2 — IR and lowering through zumir

**Issue.** [#4](https://github.com/pedrobtz/zupmxode/issues/4), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `full — first compiled code holding an external pointer with a
finalizer (the MIR context), and the first link to zumir's C API:
valgrind, ASan and a dispatched step-20 gctorture before merging.`

**Goal.** A hand-built IR model runs on both tiers and gives the same
doubles on each; the riskiest boundary is crossed before anything is
built on it.

**Work.**

- The IR (§4.1): integer vector plus constant pool, a verifier in
  `src/core/` that rejects malformed IR, and serialisation of the IR
  only (D18).
- An internal R builder for hand-written IR, used by tests until Stage 3.
- Lowering of each block to one MIR function through `zumir`'s API,
  the math table imported with `zumir_import()`; the model object owns
  its context and frees it in a finalizer.
- `pmx_ode_tier()`; on the interpreter tier, one context per thread
  (zumir D10) is left to Stage 8, which owns threads.
- Where `zumir_tier()` is `"none"`, a classed error naming it (D15).
- Fuzz harness over the IR verifier, built from `src/core/` without R,
  and the `fuzz.yml` caller with a `reach-floor`.

**Tests.**

- Tiers: every hand-built model's blocks give identical doubles on
  native and interpreter, in `test-tiers.R`.
- Math table: each function against R's own on a grid that includes
  `NaN`, `±Inf`, `-0` and the domain edges, on both tiers, in
  `test-math-table.R`.
- Serialised IR: `saveRDS()`/`readRDS()` round trip gives identical
  output, and a truncated or altered IR is refused by the verifier with a
  classed error, never lowered, in `test-ir-serialise.R`.
- Model object: `rm()` then `gc()` frees the context (finalizer runs
  once), a restored model whose context is gone regenerates, and a model
  used in a forked child gives the parent's output, in `test-model-object.R`.
- Verifier fuzz harness: no input crashes it, and every input it accepts
  lowers, in `fuzz/fuzz_ir_verify.c`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.
- A full-profile run green on the pull request, and a dispatched step-20
  gctorture run green on its last commit.
- `assurance.md` rows for the serialised IR and the model object filled.

**Design.** [§4.1](design.md#41-parser-and-compiler-r),
[§4.7](design.md#47-lowering-and-execution-tiers-c-through-zumir),
[§6.3](design.md#63-numerics), D13, D14.

## Stage 3 — Parser and compiler

**Issue.** [#5](https://github.com/pedrobtz/zupmxode/issues/5), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** Every corpus model written in the §6.1 subset compiles to IR,
classified, folded and hoisted, and runs on both tiers.

**Work.**

- `pmx_ode_model()`: the walker over R's parse, name resolution, the
  classifier (closed-form, linear, general).
- Constant folding, common subexpressions, dead code, `ivl` hoisting;
  each can be switched off internally.
- Compile errors naming the construct, its position and the nearest
  alternative (§8).

**Tests.**

- Each construct of §6.1 compiles and evaluates to R's own value of the
  same expression, in `test-language.R`.
- Each unsupported construct in the corpus's rxode2-only features errors
  with construct, position and alternative, in `test-compile-errors.R`.
- Classifier: each corpus model's class against a hand label, in
  `test-classify.R`.
- Optimisations: every corpus model gives identical block output with
  folding, CSE and hoisting each switched off, in `test-optimise.R`.
- Subset promise: every model `pmx_ode_model()` accepts also compiles
  under rxode2, when installed, in `test-language.R`.

**Exit criteria.**

- IR builds for every corpus model in the §6.1 subset.
- Every Tests line exists, uses its oracle, and was seen to fail.

**Design.** [§4.1](design.md#41-parser-and-compiler-r),
[§6.1](design.md#61-the-language), [§8](design.md#8-errors), D3.

## Stage 4 — Event scheduler

**Issue.** [#6](https://github.com/pedrobtz/zupmxode/issues/6), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** Event tables become per-subject timelines with PREDPP's
semantics, checked rule by rule before any solver uses them.

**Work.**

- `pmx_ode_events()` and validation of NONMEM-style columns, errors
  naming the row and the rule (§8).
- The scheduler in `src/core/`: dose starts, infusion ends, lags,
  `II`/`ADDL` expansion, `SS` 1 and 2, `EVID` 0–4, `RATE` −1 and −2,
  negative `CMT`, covariate changes with LOCF.
- An internal accessor returning a subject's timeline, for tests.
- Fuzz harness over the scheduler, built without R.

**Tests.**

- Each §6.2 rule: one case with a hand-computed timeline, in
  `test-events-rules.R`.
- Each event-table error: the row and the rule named, in
  `test-events-errors.R`.
- `ADDL` expansion: counts at 0, 1 and 10³ doses give the expected
  number of timeline entries, in `test-events-rules.R`.
- Scheduler fuzz harness: every timeline it produces is sorted and its
  dose count equals the expansion of the input, in
  `fuzz/fuzz_scheduler.c`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.
- The scheduler harness runs in `fuzz.yml` above its reach floor.

**Design.** [§4.2](design.md#42-event-scheduler-c),
[§6.2](design.md#62-events), D4.

## Stage 5 — Closed-form solver

**Issue.** [#7](https://github.com/pedrobtz/zupmxode/issues/7), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** `pmx_ode_solve()` returns predictions for closed-form models,
matching both oracles, with exact sums.

**Work.**

- `pmx_ode_solve()` for one or more subjects, single-threaded.
- One-, two- and three-compartment closed forms with bolus, zero- and
  first-order input, in the three parameterisations; superposition and
  steady state through a `zuexact` accumulator.
- The three-engine record of the §6.2 corner cases starts here.

**Tests.**

- Every closed-form corpus model against rxode2 `linCmt()` and
  mrgsolve `PKMODEL` fixtures, in `test-closed-form.R`.
- Superposition over 10³ doses: the `zuexact` sum against a long-double
  reference, and the naive sum shown further off, in
  `test-superposition.R`.
- Steady state: the exact `SS=1` state equals repeated dosing to
  convergence, in `test-closed-form.R`.
- The §6.2 corner cases recorded against both fixtures, disagreements
  listed, in `test-events-engines.R`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.
- `assurance.md` boundary rows for the closed forms and PREDPP events
  filled.

**Design.** [§4.3](design.md#43-solvers-c), [§4.5](design.md#45-r-layer),
[§6.3](design.md#63-numerics).

## Stage 6 — LSODA and Dormand–Prince

**Issue.** [#8](https://github.com/pedrobtz/zupmxode/issues/8), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `full — vendors liblsoda: sanitizers, valgrind, rchk and every
platform leg on code a re-entrant per-thread solver is built from.`

**Goal.** General models solve and match both oracles at the §6.3
tolerances.

**Work.**

- `liblsoda` vendored under `src/vendor/` with its licence, pinned
  commit and every patch recorded in `CLAUDE.md`; made re-entrant per
  instance.
- Dormand–Prince 5(4) with dense output, written for the package.
- Steady state for general models by repeated dosing, with a cap and a
  warning; solver failure as `NA` for the subject with a warning naming
  subject and time.

**Tests.**

- Every general corpus model against both fixtures at §6.3 tolerances,
  for each solver, in `test-ode.R`.
- Every linear corpus model: closed form against LSODA, in
  `test-ode.R`.
- Dense output against the solution at the same time stepped to
  directly, in `test-ode.R`.
- Steady-state cap: convergence on the last allowed iteration, the
  warning one iteration past it, in `test-steady-state.R`.
- A stiff model forced to fail: `NA` for that subject only, the others
  intact, the warning naming subject and time, in `test-solver-failure.R`.
- Two LSODA instances interleaved give the results each gives alone, in
  `test-ode.R`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.
- A full-profile run green on the pull request.

**Design.** [§4.3](design.md#43-solvers-c), [§6.3](design.md#63-numerics),
[§8](design.md#8-errors), D5.

## Stage 7 — Matrix exponential and solver selection

**Issue.** [#9](https://github.com/pedrobtz/zupmxode/issues/9), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** Linear systems use the matrix exponential, every model gets the
most specific solver, and the user can see and override the choice.

**Work.**

- Padé with scaling and squaring, recomputed on parameter or covariate
  change, with the constant-infusion particular solution; the state
  limit (about 20, fixed here) above which the classifier falls back.
- Automatic selection; `pmx_ode_method()`; `pmx_ode_info()` with
  states, parameters, solver, tier and lowering time.

**Tests.**

- Every linear corpus model: matrix exponential against LSODA at tight
  tolerance and against the closed form where one exists, in
  `test-expm.R`.
- State limit: a linear model at the limit uses the matrix exponential,
  one state above it falls back to general, in `test-expm.R`.
- Selection: each class gets its documented solver, and `method=`
  overrides it, as `pmx_ode_info()` reports, in `test-select.R`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.

**Design.** [§4.3](design.md#43-solvers-c), D6.

## Stage 8 — Parallel runtime

**Issue.** [#10](https://github.com/pedrobtz/zupmxode/issues/10), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `full — C running in OpenMP worker threads, which must never
touch the R API (threadcheck), and the legs built without OpenMP.`

**Goal.** Subjects run across threads with output identical to one
thread, on both tiers, with nothing allocated on the R heap.

**Work.**

- Per-thread arenas, a preallocated output matrix, subjects split across
  threads, interrupts checked between subjects on the main thread.
- On the interpreter tier, one MIR context per thread, compiled before
  the parallel region (zumir D10).
- `pmx_ode_threads()`; the serial build where OpenMP is absent.

**Tests.**

- 1 against 8 threads: bitwise identical output on every corpus model,
  on each tier, in `test-threads.R`.
- Allocation: `Rprofmem()` records no allocation between entering and
  leaving the parallel region, in `test-alloc.R`.
- Interrupt: an interrupt raised during a long run returns control
  without a crash, in `test-threads.R`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.
- A full-profile run green on the pull request, threadcheck included.

**Design.** [§4.4](design.md#44-parallel-runtime-c),
[§4.7](design.md#47-lowering-and-execution-tiers-c-through-zumir), D8.

## Stage 9 — Population simulation

**Issue.** [#11](https://github.com/pedrobtz/zupmxode/issues/11), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** `pmx_ode_sim()` simulates a population reproducibly from
`set.seed()`, and NONMEM-style data frames run as events.

**Work.**

- `pmx_ode_sim()`: etas from omega and residual error drawn in R; the
  four error models (D19), `prop()` and its siblings.
- `pmx_ode_solve()` on a NONMEM-style data frame.

**Tests.**

- Reproducibility: the same seed gives identical results at 1 and 8
  threads, in `test-sim.R`.
- Error models: each of the four against its formula applied by hand to
  `pmx_ode_solve()` output with the same draws, in `test-sim.R`.
- Etas: the sample covariance of 10⁴ subjects' etas within sampling
  error of omega, in `test-sim.R`.
- NONMEM input: a corpus data set solved through the data frame and
  through `pmx_ode_events()` gives identical output, in
  `test-nonmem-data.R`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.

**Design.** [§4.5](design.md#45-r-layer),
[§5](design.md#5-public-r-api-v010-shape), D7.

## Stage 10 — Per-subject evaluation and the C API

**Issue.** [#12](https://github.com/pedrobtz/zupmxode/issues/12), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** An estimator can evaluate one subject's predictions and
residual variances, from R or C, without allocation.

**Work.**

- `pmx_ode_subjects()`, `pmx_ode_subject_eval()` in R and
  `zupmxode_subject_eval()` in C,
  `zupmxode_arena_new()`, the versioned options struct, status codes.
- `inst/include/zupmxode.h`, exported through `R_RegisterCCallable`.

**Tests.**

- Every corpus model: `pmx_ode_subject_eval()` predictions identical to
  `pmx_ode_solve()`'s, and `var` equal to the error model's formula, in
  `test-subject-eval.R`.
- The same call from 8 threads on different subjects of one model gives
  the serial results, through a C test caller, in `test-subject-eval.R`.
- An unknown options version and a too-small buffer return their status
  codes and never raise an R error, in `test-c-api.R`.
- The header compiles standalone as C99, in `tests/c/`.

**Exit criteria.**

- Every Tests line exists, uses its oracle, and was seen to fail.
- `assurance.md` row for the prepared subjects object filled.

**Design.** [§4.6](design.md#46-per-subject-evaluation-c-api), D9, D12.

## Stage 11 — Benchmarks and README

**Issue.** [#13](https://github.com/pedrobtz/zupmxode/issues/13), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** Every §10 row is measured, and the README and Get started page
show the package as it runs.

**Work.**

- Every §10 row recorded in `.agents/benchmarks.md`, met or not.
- README: what the package is, the tier each r-actions leg runs, the
  interpreter figure.
- `vignettes/zupmxode.Rmd` evaluated (`pkg-docs`).

**Tests.**

- The vignette knits with `eval = TRUE` in R CMD check.

**Exit criteria.**

- The compile-time and candidate-loop rows of §10 are met.
- R CMD check passes with the vignette evaluated.

**Design.** [§10](design.md#10-performance-targets),
[§12](design.md#12-acceptance-criteria-for-v010).

## Stage 12 — Assurance

**Issue.** [#14](https://github.com/pedrobtz/zupmxode/issues/14), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `quick`.

**Goal.** No new features. Every required check is wired, every row of
`.agents/assurance.md` filled, and the gap lists read.

**Work.**

- `alloc-failure` and, with the C tests, `mutation` on a schedule;
  `fuzz` with its `reach-floor` for both harnesses.
- Mutation survivors and the fuzz-reach and native-coverage gap lists
  read; the ones missing an oracle filed as issues.

**Exit criteria.**

- `assurance.yml` green.
- Every remaining `gap #<issue>` closed or moved under "After 0.1.0" with
  a reason.

**Design.** [§9](design.md#9-testing), [§12](design.md#12-acceptance-criteria-for-v010).

## Stage 13 — Full CI

**Issue.** [#15](https://github.com/pedrobtz/zupmxode/issues/15), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `full`.

**Goal.** No new features. The full r-actions profile finds nothing on
the default branch's code.

**Work.**

- Run the full profile (`full-ci`) and fix what it finds: valgrind, ASan,
  every platform leg, `nosuggests`, `nold`, cran-special.
- Dispatch gctorture at step 20 on the branch; fix it and any open
  `gctorture` issue.

**Exit criteria.**

- A full-profile run green on the pull request, a dispatched one green on
  the merge commit, and a dispatched step-20 gctorture run green on the
  pull request's last commit.

**Design.** [§12](design.md#12-acceptance-criteria-for-v010).

## Stage 14 — v0.1.0 readiness

**Issue.** [#16](https://github.com/pedrobtz/zupmxode/issues/16), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `full`.

**Goal.** The package is versioned and documented for submission.

**Work.**

- `Version: 0.1.0` and its `NEWS.md` heading.
- `cran-comments.md` listing the r-actions legs as test environments.
- `urlchecker::url_check()` and `spelling::spell_check_package()` clean.
- `zumir` on CRAN and `Remotes:` removed; `dazdiff` out of `Suggests`
  unless it is on CRAN (D16).

**Exit criteria.**

- A full-profile r-actions run green on the merge commit.

**Design.** [§12](design.md#12-acceptance-criteria-for-v010), D11.

## Stage 15 — Review

**Issue.** [#17](https://github.com/pedrobtz/zupmxode/issues/17), sub-issue of
[v0.1.0](https://github.com/pedrobtz/zupmxode/issues/1).

**CI.** `full`.

**Goal.** No new features. The release candidate passes a critical
review and a CRAN review.

**Work.**

1. `pkg-critical-review`: its blocking and required findings become
   sub-issues of this stage's issue, fixed together in one pull request
   that closes each.
2. `pkg-cran-review`, run after those fixes, every FAIL and WARN
   resolved.

**Exit criteria.**

- Every blocking and required review issue closed, the `pkg-cran-review`
  table clean, `NEWS.md` reflecting any user-visible fix, and a
  full-profile run green on the default branch, whose head becomes the
  release commit.

**Design.** [§12](design.md#12-acceptance-criteria-for-v010).

## After 0.1.0

- [§14](design.md#14-later-versions-in-the-order-they-are-expected): all of
  it, as the design orders it; v0.1.0's scope is §2.
- The Michaelis–Menten closed form (§4.3): waits for zuspecial#61 and
  `zuspecial` on CRAN (Q5); until then the model classifies as general.
- Draws in C through `zurand` (D7): v0.2.0.
- An IR interpreter for systems where `zumir_tier()` is `"none"` (D15).

## Blocking open questions

- [Q4](design.md#15-open-questions) blocks Stage 14: CRAN or a GitHub
  release; this roadmap assumes CRAN.
