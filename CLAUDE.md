# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Shared conventions

The portfolio’s shared rules (`workflow`, `r-package`, `ci`,
`native-code`, `correctness`, `assurance`) come from the
[r-skills](https://github.com/pedrobtz/r-skills) plugin’s `rules/`,
loaded at session start because of the marker above. Never copy them
here. Where this file and a shared rule disagree, this file wins for
this package. Cloud sessions do not load the plugin and run without
them.

## What this is

zupmxode simulates population PK/PD models, written in a subset of the
rxode2 model language, over NONMEM-style event records, in parallel over
subjects, with no C compiler on the machine. Models are parsed in R,
compiled to the package’s own IR, and turned into native code (or run in
MIR’s interpreter) through `zumir`. It is a simulation engine with a C
API for estimators; it is not an estimator, not a new model language,
and imports nothing but `zumir` (design D20).

The specification is
[`.agents/design.md`](https://pedrobtz.github.io/zupmxode/.agents/design.md),
adopted from RFC 0047. It outranks this file: `CLAUDE.md` orients, the
design decides. The whole package rests on D2 (MIR through `zumir`),
which is provisional on `zumir`’s Stage 0 spike; §11.1 of the design
says what changes if the spike fails.

## Current state

Stage 0 (package identity and CI) done 2026-11-15; no engine
functionality yet. The roadmap in `.agents/roadmap.md` takes it to
v0.1.0.

## Versions and dependents

The current version’s parent issue is
[\#1](https://github.com/pedrobtz/zupmxode/issues/1) for `v0.1.0`.
Siblings that will depend on it: `zupmxfit` (RFC 0072, through
`zupmxode_subject_eval()`) and `zupmxio` (RFC 0046, `Suggests` both
ways). It will link to `zumir` and `zuexact`; `zumir` is not on CRAN
yet, so it comes through `Remotes:` until it is, and its release gates
this package’s submission.

## Commands

- `make -C tests/c check`: the C core built and tested without R’s
  headers. Part of the default loop once a stage touches `src/core/`.
- `bash tools/ci/platform.sh <os>`: every per-OS check, as CI runs it.
- `Rscript tools/make-corpus.R`: regenerates the oracle fixtures (about
  2 minutes; needs rxode2, mrgsolve and nlmixr2lib). Occasional: when
  the corpus or an oracle changes.
- `tools/run-benchmarks [engine]`: the §10 rows, appended to
  `.agents/benchmarks.csv`; run alone on an idle machine.

## Architecture

See `.agents/design.md` §4 for the whole. Built so far:

- `src/core/ir.[ch]`: the IR format, its verifier, a reference evaluator
  (the tests’ oracle and the fuzz harness’s, never an execution tier)
  and the math table. No R headers.
- `src/zpx_lower.c`: IR to MIR, one module and one function per block,
  through `zumir`’s C API and emitter; `zumir_import()` registers the
  table’s functions the C library lacks.
- `src/zpx_model.c`: the model’s external pointer and finalizer, and the
  `.Call` entries. Failures of the IR or of MIR come back as status
  lists; `R/model.R` raises the classed condition.
- `R/compile.R`: the compiler, R’s parse to IR text for the assembler:
  statements, names, `linCmt()`’s expansion, the classifier, and one
  generator for every block (DCE by statement, folding through the math
  table, CSE by value, `ivl` hoisting; `compile_model(opts =)` switches
  each off for the tests). `R/pmx_ode_model.R` is the export.
- `R/ir.R`: the IR assembler; `R/model.R`: tiers, lowering, regeneration
  after [`readRDS()`](https://rdrr.io/r/base/readRDS.html) (D18).
- `src/zupmxode.c`: registration only.

## Invariants that are easy to break

These are fixed by the design before any code exists; each names the
section that decides it.

- No `Imports` but `zumir` (§2, §12, D20). Another runtime dependency
  breaks the “stock R installation” promise the package exists for.
  `zumir` is there only because `R_GetCCallable()` needs its namespace
  loaded.
- Nothing under `src/core/` includes R’s headers (`tests/c` fails to
  build if it does). The IR verifier in `src/core/ir.c` is the only
  thing between a saved IR and lowered code: every new opcode or operand
  needs its range check there, its refusal test in `test-ir-verify.R`,
  and a seed in `tests/fuzz/corpus/`.
- C code never calls `error()` and never touches an R object inside the
  parallel region or `zupmxode_subject_eval()`; it returns status codes
  and writes into caller-owned buffers (§4.4, §4.6, §8). Breaking it
  makes the engine unsafe for estimators calling it from threads.
- Random draws happen only in R, with R’s RNG (D7). C is deterministic,
  so results are reproducible from
  [`set.seed()`](https://rdrr.io/r/base/Random.html) at any thread
  count.
- MIR is reached only through `zumir`’s C API and never vendored here
  (D14); the IR is the package’s own and nothing outside the lowering
  depends on MIR (D13).
- The two tiers must give identical doubles, and output must be bitwise
  identical across thread counts (§6.3, §9). A difference is a bug
  unless traced to a recorded contraction.

## Testing

Conformance suite: the oracle fixtures in
`tests/testthat/fixtures/oracle/`, written by `tools/make-corpus.R`
(design §7 says how). Regenerate them, never edit them; a new
`disagreements.csv` row needs its note in `disagreements.md`, or
`test-oracle-fixtures.R` fails.

Helpers:

- `helper-tables.R`: `expect_table_equal()` at the §6.3 tolerances,
  through `dazdiff` when installed; `table_differences()` is the local
  route, tested on its own.
- `helper-oracle.R`: reads fixtures and solves with rxode2. The
  generator sources it, so tests and fixtures solve the same way.
- `helper-models.R`: hand-built IR models, each with its answer written
  in R; `available_tiers()` lists the tiers this session runs.
- `helper-compile.R`: `run_blocks()` runs a compiled model’s blocks as
  the solvers will; `r_eval_model()` evaluates the same statements in R,
  the compiler’s oracle.

`fixtures/math-reference.csv` is the math table’s oracle, from
`tools/make-math-reference.R` (needs Rmpfr); `tests/fuzz/corpus/` the
fuzz seeds, from `tools/make-fuzz-corpus.R`.

The suite runs shuffled, top-level code included: build test data inside
`test_that()` or in a helper, never at file level. Runtime: about 20 s
locally, most of it rxode2 compiling three models.

## Continuous integration

Shared r-actions workflows, called from `.github/workflows/`:
`R-CMD-check.yaml` (with the Windows R-devel row), `native-checks.yaml`
(valgrind, LTO, rchk, analyzers, cran-special; sanitizers from Stage 2,
r-actions#68), `gctorture.yaml` (wide step on pull requests, step 20
daily and on dispatch), `coverage.yaml` (`native: true`),
`assurance.yaml` and `platforms.yaml`. `pkgdown.yaml` builds and deploys
the site.

`platforms.yaml` runs `tools/ci/platform.sh`: today the C core built and
tested without R (`make -C tests/c`), Ubuntu only until Stage 2 adds
macOS and Windows. Dispatch step-20 gctorture on the branch
(`gh workflow run gctorture.yml --ref <branch>`) before merging a change
to `PROTECT`s, finalizers or external pointers. The roadmap’s `full`
stages are 2, 6, 8 and the three closing ones.

## Vendored native code

`liblsoda` (`sdwfrost/liblsoda`, MIT), made re-entrant per thread (D5).
Not vendored yet: record the pinned commit and every local patch here
when it is.
