#!/usr/bin/env Rscript
# The workload of the alloc-failure job (.github/workflows/alloc-failure.yaml):
# each run fails one allocation, and must end in an R error naming it, never
# a crash. `baseline` stops before the workload, so the sweep starts at the
# package's own allocations: the solve's records, LSODA's workspace,
# prepared subjects, arenas, and the interpreter tier's model copy.
#
#   Rscript tools/alloc-exercise.R baseline
#   Rscript tools/alloc-exercise.R
#
# The run checks its results, so a failure absorbed into a wrong answer
# exits non-zero rather than passing.
mode <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(zupmxode))
ev <- pmx_ode_events(amt = 100, cmt = 1L, time = c(1, 4, 12), id = 1)
ev2 <- ev
ev2$ID <- 2
ev <- rbind(ev, ev2)
ev$wt <- c(rep(60, 4), rep(80, 4))
p <- c(ka = 0.8, vmax = 10, km = 2)
# The model is compiled under the floor too: R's parser and regular
# expressions, which compiling runs, crash on a failed allocation (an R
# finding, not the package's), and would take the sweep's whole window.
m <- pmx_ode_model("d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - vmax * central / (km + central) * wt / 70")
if (identical(mode, "baseline")) quit(save = "no")

out <- pmx_ode_solve(m, ev, p)
stopifnot(nrow(out) == 6L, all(is.finite(out$central)))
subs <- pmx_ode_subjects(m, ev)
r <- pmx_ode_subject_eval(m, subs[[2]], p, sigma = list(central = add(1)))
stopifnot(identical(r$pred, out$central[4:6]))
mi <- pmx_ode_tier(m, "interp")
stopifnot(identical(pmx_ode_solve(mi, ev, p)$central, out$central))
cat("alloc-exercise: ok\n")
