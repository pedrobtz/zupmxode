# Oracle disagreements

Written by hand; `tools/make-corpus.R` never touches this file. Every case
in `disagreements.csv` (rxode2 and mrgsolve further apart than the §6.3
tolerances) has a section here, headed by its case name, giving the PREDPP
reading zupmxode follows and why (design §6.2, §7).
`test-oracle-fixtures.R` fails on a case without one.

None as of 2026-11-15: on every paired case of the corpus, rxode2 5.1.8 and
mrgsolve 2.0.1 agree within rtol 1e-6 and atol 1e-8, solved at rtol 1e-10.
That agreement rests on four settings `tools/make-corpus.R` makes explicit,
each of which produced a disagreement before it was set:

- mrgsolve's covariates are read last-observation-carried-forward
  (`nocb = FALSE`); its default carries them backward, against PREDPP.
- Both engines' steady-state iterations run at the solve's own
  tolerances; rxode2's default stops at a relative change of 1e-6, which
  left `pk2_oral__ss1_bolus` 3e-6 away from mrgsolve.
- rxode2 returns `EVID 2` records as output rows; they are dropped so both
  engines report the observation records only.
- Doses sort before observations at the same time, so an observation at a
  dose time is read after the dose in both.

## zupmxode against the oracles

Cases where zupmxode's answer is not an oracle's, each with the reading
zupmxode follows. `test-ode.R` asserts every one.

### lib_indirect_1cpt_inhi_kin_r0rmaxcrmax

### lib_indirect_1cpt_inhi_kout_r0rmaxcrmax

Recorded 2026-11-16 (roadmap Stage 6). With the library's initial
parameters, `imax <- (r0 - rmax) / r0` makes `(1 - imax) * r0` equal to
`rmax`, so `IC50 <- crmax * (rmax - (1 - imax) * r0) / (r0 - rmax)` is
exactly 0, as R computes it too (-0). The effect's derivative then holds
`Cc / (Cc + IC50)`, which is 0 / 0 at t = 0: NaN. zupmxode evaluates the
model as written, so its solver fails at the start and the subject is NA
with a warning (design §8). rxode2's C code, compiled with reordering,
found a tiny nonzero IC50 and went on; its numbers are an accident of
rounding, not the model's answer. mrgsolve has no port of these models.
