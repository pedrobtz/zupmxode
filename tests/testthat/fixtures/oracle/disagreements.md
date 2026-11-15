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
