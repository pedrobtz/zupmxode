# CLAUDE.md

<!-- markdownlint-disable-next-line MD013 -->
This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Shared conventions

<!-- r-skills: shared-rules -->

The portfolio's shared rules (`workflow`, `r-package`, `ci`, `native-code`,
`correctness`, `assurance`) come from the
[r-skills](https://github.com/pedrobtz/r-skills) plugin's `rules/`, loaded at
session start because of the marker above. Never copy them here. Where this
file and a shared rule disagree, this file wins for this package. Cloud
sessions do not load the plugin and run without them.

## What this is

zupmxode simulates population PK/PD models, written in a subset of the
rxode2 model language, over NONMEM-style event records, in parallel over
subjects, with no C compiler on the machine. Models are parsed in R,
compiled to the package's own IR, and turned into native code (or run in
MIR's interpreter) through `zumir`. It is a simulation engine with a C API
for estimators; it is not an estimator, not a new model language, and
imports nothing but `zumir` (design D20).

The specification is [`.agents/design.md`](.agents/design.md), adopted from
RFC 0047. It outranks this file: `CLAUDE.md` orients, the design decides.
The whole package rests on D2 (MIR through `zumir`), which is provisional
on `zumir`'s Stage 0 spike; §11.1 of the design says what changes if the
spike fails.

## Current state

Stages 0 to 10 done (2026-11-16): models compile, solve and simulate
populations, and estimators evaluate subjects through the C API. The roadmap in `.agents/roadmap.md` takes it to v0.1.0.

## Versions and dependents

The current version's parent issue is
[#1](https://github.com/pedrobtz/zupmxode/issues/1) for `v0.1.0`.
Siblings that will depend on it: `zupmxfit` (RFC 0072, through
`zupmxode_subject_eval()`) and `zupmxio` (RFC 0046, `Suggests` both
ways). It links to `zumir` only (zuexact
left with design D21); `zumir` is not on CRAN yet, so it comes through `Remotes:`
until it is, and its release gates this package's submission.

## Commands

- `make -C tests/c check`: the C core built and tested without R's
  headers. Part of the default loop once a stage touches `src/core/`.
- `bash tools/ci/platform.sh <os>`: every per-OS check, as CI runs it.
- `Rscript tools/make-corpus.R`: regenerates the oracle fixtures (about
  2 minutes; needs rxode2, mrgsolve and nlmixr2lib). Occasional: when the
  corpus or an oracle changes.
- `tools/run-benchmarks [engine]`: the §10 rows, appended to
  `.agents/benchmarks.csv`; run alone on an idle machine.

## Architecture

See `.agents/design.md` §4 for the whole. Built so far:

- `src/core/ir.[ch]`: the IR format, its verifier, a reference evaluator
  (the tests' oracle and the fuzz harness's, never an execution tier) and
  the math table. No R headers.
- `src/core/sched.[ch]`: the event scheduler, one subject's records to a
  sorted timeline (design §6.2 has the rules); `src/zpx_sched.c` its
  `.Call`; `R/events.R` `pmx_ode_events()`, the table's checks and
  `sched_timeline()`, the tests' view of a timeline.
- `src/core/solve.[ch]`: one subject's timeline walked with the closed
  form (`src/core/linear.[ch]`, D21), the matrix exponential
  (`src/core/expm.[ch]`, probing A and c from rhs), LSODA or
  Dormand-Prince (`src/core/dopri.[ch]`, dense output across runs of
  observations), steady state, and a failure as NA from where it
  happened. No R headers: it runs on worker threads.
- `src/zpx_solve.c`: the parallel runtime (design §4.4), every subject
  split across OpenMP threads with per-thread arenas and, on the
  interpreter tier, per-thread MIR contexts; `R/solve.R` is
  `pmx_ode_solve()`, `pmx_ode_threads()`, `auto_method()` (the
  selection) and `subject_method()` (a switched-off compartment goes to
  LSODA). `R/method.R`: `pmx_ode_method()`, `pmx_ode_info()`.
- `src/zpx_lower.c`: IR to MIR, one module and one function per block,
  through `zumir`'s C API and emitter; `zumir_import()` registers the
  table's functions the C library lacks.
- `src/zpx_model.c`: the model's external pointer and finalizer, and the
  `.Call` entries. Failures of the IR or of MIR come back as status lists;
  `R/model.R` raises the classed condition.
- `R/sim.R`: `pmx_ode_sim()`, the etas and residual draws in R (D7) put
  on the table as covariate columns for one `pmx_ode_solve()`, and the
  error models `add()`, `prop()`, `lnorm()` and `add() + prop()` (D19).
- `src/zpx_eval.[ch]`: the C API (design §4.6, `inst/include/zupmxode.h`):
  prepared subjects, arenas and `zupmxode_subject_eval()`, on any
  thread; `src/zpx_api.c` its R glue, the registered `zupmxode_api`
  table and the tests' C callers; `R/subjects.R` is
  `pmx_ode_subjects()` and `pmx_ode_subject_eval()`.
- `R/compile.R`: the compiler, R's parse to IR text for the assembler:
  statements, names, `linCmt()`'s expansion, the classifier, and one
  generator for every block (DCE by statement, folding through the math
  table, CSE by value, `ivl` hoisting; `compile_model(opts =)` switches
  each off for the tests). `R/pmx_ode_model.R` is the export.
- `R/ir.R`: the IR assembler;
  `R/model.R`: tiers, lowering, regeneration after `readRDS()` (D18).
- `src/zupmxode.c`: registration only.

## Invariants that are easy to break

These are fixed by the design before any code exists; each names the
section that decides it.

- No `Imports` but `zumir` (§2, §12, D20). Another runtime dependency
  breaks the "stock R installation" promise the package exists for.
  `zumir` is there only because `R_GetCCallable()` needs its namespace
  loaded.
- Nothing under `src/core/` includes R's headers (`tests/c` fails to
  build if it does). The IR verifier in `src/core/ir.c` is the only thing
  between a saved IR and lowered code: every new opcode or operand needs
  its range check there, its refusal test in `test-ir-verify.R`, and a
  seed in `tests/fuzz/corpus/`.
- C code never calls `error()` and never touches an R object inside the
  parallel region or `zupmxode_subject_eval()`; it returns status codes
  and writes into caller-owned buffers (§4.4, §4.6, §8). Breaking it
  makes the engine unsafe for estimators calling it from threads.
- Random draws happen only in R, with R's RNG (D7). C is deterministic,
  so results are reproducible from `set.seed()` at any thread count.
- MIR is reached only through `zumir`'s C API and never vendored here
  (D14); the IR is the package's own and nothing outside the lowering
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
- `helper-oracle.R`: reads fixtures and solves with rxode2. The generator
  sources it, so tests and fixtures solve the same way.
- `helper-models.R`: hand-built IR models, each with its answer written
  in R; `available_tiers()` lists the tiers this session runs.
- `helper-compile.R`: `run_blocks()` runs a compiled model's blocks as
  the solvers will; `r_eval_model()` evaluates the same statements in R,
  the compiler's oracle.

`fixtures/math-reference.csv` is the math table's oracle, from
`tools/make-math-reference.R` (needs Rmpfr); `tests/fuzz/corpus/` and
`tests/fuzz/corpus-scheduler/` the fuzz seeds, from
`tools/make-fuzz-corpus.R`.

The suite runs shuffled, top-level code included: build test data inside
`test_that()` or in a helper, never at file level. Runtime: about 20 s
locally, most of it rxode2 compiling three models.

## Continuous integration

Shared r-actions workflows, called from `.github/workflows/`:
`R-CMD-check.yaml` (with the Windows R-devel row), `native-checks.yaml`
(valgrind, LTO, rchk, analyzers, cran-special; sanitizers from Stage 2,
r-actions#68; threadcheck from Stage 8, for C on OpenMP worker threads),
`gctorture.yaml` (wide step on pull requests, step 20 daily and on
dispatch), `coverage.yaml` (`native: true`), `assurance.yaml` and
`platforms.yaml`. `pkgdown.yaml` builds and deploys the site.

`platforms.yaml` runs `tools/ci/platform.sh`: today the C core built and
tested without R (`make -C tests/c`), Ubuntu only until Stage 2 adds
macOS and Windows. Dispatch step-20 gctorture on the branch
(`gh workflow run gctorture.yml --ref <branch>`) before merging a change
to `PROTECT`s, finalizers or external pointers. The roadmap's `full`
stages are 2, 6, 8 and the three closing ones.

## Vendored native code

`liblsoda` (<https://github.com/sdwfrost/liblsoda>, MIT), pinned at
`16e49b428482e2efb552b0b4816a8ff8221e7efe` (upstream's last commit, 2018),
in `src/vendor/liblsoda/`. `tools/vendor/liblsoda.sh` makes the tree from
that commit and is the only place it is changed; `tools/vendor/verify`
(the `vendor` workflow) rebuilds it and compares byte for byte. The
patches, all in that script (8 to 11 found by Stage 6's full run):

1. No console output: `printf`/`fprintf(stderr, …)` become `LSODA_LOG`,
   which does nothing (the messages are in `ctx->error`).
2. Only the sources upstream's Makefile builds (not `printcf.c`, nor the
   `#error` stub `ewset.c`).
3. `_strdup_printf()` allocates the NUL's byte and checks `malloc()`.
4. The BLAS routines' K&R definitions become prototypes (C23).
5. `lsoda()`'s `ihit` starts at 0, and `softfailure()`'s loop variable no
   longer shadows the `i` its message prints (an uninitialised read).
6. `(void)` for empty parameter lists.
7. `solsy()` returns 0 rather than `abort()` for the banded Jacobian
   liblsoda never builds (CRAN: no C that can end R).
8. A stray comma in `lsoda()`'s message for a repeatedly failed error test
   made `%g` read a `char *` and `%s` a double: a crash on that path.
9. `ERROR()`, `hardfailure()` and `softfailure()` are C99 variadic macros,
   not GNU's `, ## __VA_ARGS__` (a significant warning under clang 23).
10. `_strdup_printf()` has printf's format attribute, so the compiler
    checks every message, measures into a one-byte buffer, not NULL, and
    returns NULL for a NULL format.
11. The BLAS routines are `lsoda_daxpy` … `lsoda_idamax`: under the plain
    names the linker bound them to the BLAS R links (OpenBLAS's Fortran
    `idamax_`, by pointer), so an `int` was read as an address. Every
    symbol is hidden too (`$(C_VISIBILITY)`), but `R_init_zupmxode`.

To change the pin or a patch: edit the script or `manifest.tsv`, run
`sh tools/vendor/liblsoda.sh`, regenerate `tools/vendor/checksums.sha256`
(`cd src/vendor && find . -type f | LC_ALL=C sort | xargs shasum -a 256`),
and run `sh tools/vendor/verify`. One instance per thread: each solve
makes its own context. `make -C tests/c` builds it standalone and checks
two contexts in turn against each alone.
