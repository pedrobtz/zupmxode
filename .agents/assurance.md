# Assurance

How this package meets the shared `assurance.md` rules. CI checks this
file (r-actions `assurance.yml`), so keep the headings and table
columns as they are.

- **Statuses:** a cell is a fact, `gap #<issue>` or `n/a: <reason>`.
- **Kinds:** delete the ones that do not apply.
- **Tables:** delete a table only when no kind requires it. Delete the
  example rows.

Written 2026-11-15 in Stage 0, before any code beyond the scaffold existed;
Stage 2 filled the rows of the IR, the lowering, the math table and the
compiled model.
Every row comes from [design.md](design.md) and [roadmap.md](roadmap.md),
so almost every cell is a gap naming the stage issue whose pull request
writes its test. Each stage turns its rows into facts.

## Kinds

`compiled` `untrusted-input` `spec` `lifecycle`

<!-- Evidence: compiled, src/ and vendored liblsoda (design §4.3);
untrusted-input, a saved model's IR read back from disk (§4.1, D18) and
event tables from user data (§4.2); spec, PREDPP event semantics, the
closed-form solutions, LSODA and Dormand-Prince (§4.3, §6.2); lifecycle,
the compiled model holding a MIR context and the prepared subjects
(§4.6, §4.7). Not emits: results are R objects, and no bytes are written
for another program to read. -->

## Required checks

`Where` names the workflow file under `.github/workflows/` that runs the
check. Add `(own)` when that workflow is the package's own rather than a
call to r-actions.

| Check         | Status  | Where              |
| ------------- | ------- | ------------------ |
| sanitizers    | done    | native-checks.yaml |
| valgrind      | done    | native-checks.yaml |
| rchk          | done    | native-checks.yaml |
| gctorture     | done    | gctorture.yaml     |
| analyzers     | done    | native-checks.yaml |
| alloc-failure | gap #14 |                    |
| mutation      | gap #14 |                    |
| fuzz          | done    | fuzz.yaml          |
| fuzz-reach    | done    | fuzz.yaml          |

## Boundaries

One row per boundary the package crosses. An oracle knows the answer:

- a conformance suite or test vectors;
- a reference implementation;
- a round trip;
- equivalent inputs compared;
- the bytes the other side received.

"No crash" is not an oracle.

<!-- markdownlint-disable MD013 -->
| Boundary                         | Oracle                                                         | Carried by |
| -------------------------------- | -------------------------------------------------------------- | ---------- |
| serialised IR                    | round trip; altered or truncated IR refused by the verifier    | tests/testthat/test-ir-serialise.R |
| serialised IR                    | an accepted IR touches only its declared arrays (ASan)         | tests/fuzz/fuzz_ir_verify.c |
| lowering to MIR, two tiers       | native and interpreter give identical doubles                  | tests/testthat/test-tiers.R |
| lowering to MIR, two tiers       | the reference evaluator, and R's arithmetic on the same inputs | tests/testthat/test-tiers.R |
| math table                       | Rmpfr at 200 bits for finite arguments                         | tests/testthat/test-math-table.R |
| math table                       | R's own functions at poles, infinities, NaN and zero           | tests/testthat/test-math-table.R |
| rxode2 model language subset     | R's value of each construct; rxode2 compiles what we accept    | tests/testthat/test-language.R |
| PREDPP event semantics           | hand-computed timelines, rule by rule                          | tests/testthat/test-events-rules.R |
| PREDPP event semantics           | the timeline sorted, every dose expanded, every infusion ended | tests/fuzz/fuzz_scheduler.c |
| PREDPP event semantics           | rxode2 and mrgsolve on the same events (three-engine record)   | tests/testthat/test-closed-form.R |
| closed-form solutions            | rxode2 `linCmt()` and mrgsolve `PKMODEL`                       | tests/testthat/test-closed-form.R |
| closed-form solutions            | 200-bit superposition over 1,000 doses                         | tests/testthat/test-closed-form.R |
| LSODA and Dormand–Prince         | rxode2 and mrgsolve; closed form on linear models              | tests/testthat/test-ode.R |
| LSODA, vendored                  | two contexts in turn against each alone, and the exact solution | tests/c/test_lsoda.c |
| matrix exponential               | LSODA at tight tolerance; closed form where one exists         | gap #9     |
| C API `zupmxode_subject_eval()`  | `pmx_ode_solve()` on the same model and subjects               | gap #12    |
<!-- markdownlint-enable MD013 -->

## Objects

One row per object that outlives the call that made it. "Sequences" is
the test of its out-of-order use (abort then read, use after end, use
after free). The guard columns name the test, or say why the guard
cannot apply.

<!-- markdownlint-disable MD013 -->
| Object                          | Sequences | Fork    | Unwind  | Owner   |
| ------------------------------- | --------- | ------- | ------- | ------- |
| compiled model (MIR context)    | tests/testthat/test-model-object.R | tests/testthat/test-model-object.R | tests/testthat/test-model-object.R | n/a: the lowered code embeds its constants and points into no IR or R object |
| prepared subjects               | gap #12   | gap #12 | gap #12 | gap #12 |
| evaluation arena (C API)        | gap #12   | gap #12 | gap #12 | gap #12 |
<!-- markdownlint-enable MD013 -->

## Limits

One row per limit the package enforces. "Straddled by" is the test of
`limit - 1`, `limit` and `limit + 1` on every path that reaches it.

| Limit                         | Unit       | Straddled by |
| ----------------------------- | ---------- | ------------ |
| matrix-exponential state cap  | states     | gap #9       |
| general steady-state cap      | iterations | tests/testthat/test-steady-state.R |
| IR instructions per block     | instructions | tests/testthat/test-ir-verify.R |
| IR registers per block        | registers  | tests/testthat/test-ir-verify.R |
| IR states, parameters, work and outputs | slots | tests/testthat/test-ir-verify.R |
| IR constants                  | constants  | tests/testthat/test-ir-verify.R |
| timeline entries per subject  | entries    | tests/testthat/test-events-rules.R |
