# Benchmarks

The rows of [design §10](design.md#10-performance-targets), written by
`tools/run-benchmarks` from `benchmarks.csv`, which keeps every run.
This table shows the latest run of each row and engine. Regenerated,
never edited.

| Row | Target | Engine | Value | Unit | Detail | Date | Machine |
| --- | --- | --- | --- | --- | --- | --- | --- |
| compile_2cpt | under 20 ms | rxode2 | 1.365e+00 | s |  | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| closed_form_1e6 | matches rxode2 `linCmt()` | rxode2 | 3.983e+01 | s | 1e6 subjects x 1 obs, 6 threads | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_1thread | within 1.5x of rxode2 | rxode2 | 1.607e+00 | s | 1e3 subjects x 1e3 obs | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_threads | within 1.5x of rxode2 | rxode2 | 1.549e+00 | s | 1e3 subjects x 1e3 obs, 6 threads | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ejcp_obs_per_s | above 5,000 obs/s | rxode2 | 2.527e+05 | obs/s | 1e3 subjects x 70 obs, 1 thread, 0.28 s | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| candidates_200 | faster than rxode2 end to end | rxode2 | 2.054e+02 | s | compile and solve 10 subjects | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| compile_2cpt | under 20 ms | mrgsolve | 2.035e+00 | s |  | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| closed_form_1e6 | matches rxode2 `linCmt()` | mrgsolve | 1.251e+01 | s | 1e6 subjects x 1 obs, 1 thread (mrgsolve has no threads) | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_1thread | within 1.5x of rxode2 | mrgsolve | 6.830e-01 | s | 1e3 subjects x 1e3 obs | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ejcp_obs_per_s | above 5,000 obs/s | mrgsolve | 4.605e+05 | obs/s | 1e3 subjects x 70 obs, 1 thread, 0.15 s | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| candidates_200 | faster than rxode2 end to end | mrgsolve | 3.410e+02 | s | compile and solve 10 subjects | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
