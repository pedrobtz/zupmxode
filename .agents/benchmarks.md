# Benchmarks

The rows of [design §10](design.md#10-performance-targets), written by
`tools/run-benchmarks` from `benchmarks.csv`, which keeps every run.
This table shows the latest run of each row and engine. Regenerated,
never edited.

| Row | Target | Engine | Value | Unit | Detail | Date | Machine |
| --- | --- | --- | --- | --- | --- | --- | --- |
| closed_form_1e6 | matches rxode2 `linCmt()` | rxode2 | 3.983e+01 | s | 1e6 subjects x 1 obs, 6 threads | 2026-11-15 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| compile_2cpt | under 20 ms | zupmxode | 1.100e-02 | s | parse, compile and lower to native code; median of 11 | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| closed_form_1e6 | matches rxode2 `linCmt()` | zupmxode | 3.712e+00 | s | 1e6 subjects x 1 obs, 1 thread (this build of R has no OpenMP) | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_1thread | within 1.5x of rxode2 | zupmxode | 2.862e+00 | s | 1e3 subjects x 1e3 obs, expm | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_interp | recorded | zupmxode | 3.903e+00 | s | 1e3 subjects x 1e3 obs, expm | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_threads | within 1.5x of rxode2 | zupmxode | 2.717e+00 | s | 1e3 subjects x 1e3 obs, expm 1 thread (this build of R has no OpenMP) | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ejcp_obs_per_s | above 5,000 obs/s | zupmxode | 2.979e+05 | obs/s | 1e3 subjects x 70 obs, 1 thread, expm, 0.23 s | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_1thread | within 1.5x of rxode2 | zupmxode/lsoda | 4.966e+00 | s | 1e3 subjects x 1e3 obs, lsoda | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_interp | recorded | zupmxode/lsoda | 8.109e+00 | s | 1e3 subjects x 1e3 obs, lsoda | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_threads | within 1.5x of rxode2 | zupmxode/lsoda | 4.707e+00 | s | 1e3 subjects x 1e3 obs, lsoda 1 thread (this build of R has no OpenMP) | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ejcp_obs_per_s | above 5,000 obs/s | zupmxode/lsoda | 1.284e+05 | obs/s | 1e3 subjects x 70 obs, 1 thread, lsoda, 0.55 s | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| candidates_200 | faster than rxode2 end to end | zupmxode | 2.590e+00 | s | compile and solve 10 subjects | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| compile_2cpt | under 20 ms | mrgsolve | 1.666e+00 | s | median of 11 | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| closed_form_1e6 | matches rxode2 `linCmt()` | mrgsolve | 1.110e+01 | s | 1e6 subjects x 1 obs, 1 thread (mrgsolve has no threads) | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_1thread | within 1.5x of rxode2 | mrgsolve | 7.650e-01 | s | 1e3 subjects x 1e3 obs | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ejcp_obs_per_s | above 5,000 obs/s | mrgsolve | 3.977e+05 | obs/s | 1e3 subjects x 70 obs, 1 thread, 0.18 s | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| candidates_200 | faster than rxode2 end to end | mrgsolve | 3.211e+02 | s | compile and solve 10 subjects | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| compile_2cpt | under 20 ms | rxode2 | 9.460e-01 | s | median of 11 | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_1thread | within 1.5x of rxode2 | rxode2 | 1.521e+00 | s | 1e3 subjects x 1e3 obs | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ode_2cpt_threads | within 1.5x of rxode2 | rxode2 | 1.607e+00 | s | 1e3 subjects x 1e3 obs, 6 threads | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| ejcp_obs_per_s | above 5,000 obs/s | rxode2 | 2.893e+05 | obs/s | 1e3 subjects x 70 obs, 1 thread, 0.24 s | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
| candidates_200 | faster than rxode2 end to end | rxode2 | 1.979e+02 | s | compile and solve 10 subjects | 2026-11-16 | Darwin x86_64 Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz 6 cores |
