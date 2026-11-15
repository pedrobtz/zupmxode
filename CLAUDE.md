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
for estimators; it is not an estimator, not a new model language, and has
no `Imports`.

The specification is [`.agents/design.md`](.agents/design.md), adopted from
RFC 0047. It outranks this file: `CLAUDE.md` orients, the design decides.
The whole package rests on D2 (MIR through `zumir`), which is provisional
on `zumir`'s Stage 0 spike; §11.1 of the design says what changes if the
spike fails.

## Current state

Scaffolded 2026-11-15 with `create-pkg.sh` (C via `use_c()`); no
functionality yet. The roadmap in `.agents/roadmap.md` takes it to v0.1.0.
`DESCRIPTION`'s `Title`, `Description` and `Authors@R` are still the
scaffold placeholders.

## Versions and dependents

The current version's parent issue is
[#1](https://github.com/pedrobtz/zupmxode/issues/1) for `v0.1.0`.
Siblings that will depend on it: `zupmxfit` (RFC 0072, through `zupmxode_subject_eval()`) and
`zupmxio` (RFC 0046, `Suggests` both ways). It will link to `zumir` and
`zuexact`; `zumir` is not on CRAN yet, so it comes through `Remotes:`
until it is, and its release gates this package's submission.

## Commands

Only the shared commands so far.

## Architecture

To be written as stages land; see `.agents/design.md` §4. The scaffold set
up a C layer (`src/zupmxode.c`, registered through `useDynLib`).

## Invariants that are easy to break

These are fixed by the design before any code exists; each names the
section that decides it.

- No `Imports`, ever (§2, §12). A runtime dependency breaks the
  "stock R installation" promise the package exists for.
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

Conformance suite: none yet. The design's oracle is rxode2 and mrgsolve
on `zupmxio`'s corpus, `nlmixr2lib` and mrgsolve's `modlib()` (§7, §9).
Helpers: none yet. Runtime: none yet.

## Continuous integration

`R-CMD-check.yaml` and `coverage.yaml` call the shared r-actions
workflows (`r-cmd-check.yml@v1`, `coverage.yml@v1`); `pkgdown.yaml`
builds and deploys the site. Stage 9 of the design's roadmap is marked
`full`.

## Vendored native code

`liblsoda` (`sdwfrost/liblsoda`, MIT), made re-entrant per thread (D5).
Not vendored yet: record the pinned commit and every local patch here
when it is.
