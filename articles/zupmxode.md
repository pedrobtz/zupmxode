# Get started with zupmxode

``` r

library(zupmxode)
```

We take one model from text to a simulated population: compile it, give
it a week of doses, solve it for one subject, simulate a hundred, and
evaluate one subject the way an estimator does. By the end you will know
the package’s four objects (a model, an event table, a solution and a
prepared subject) and the functions that make them. How fast this is,
and where it runs, is left to the [Performance and
platforms](https://pedrobtz.github.io/zupmxode/articles/performance.html)
article.

## Compile a model

A **model** is written in rxode2’s syntax and read by R’s own parser, so
an rxode2 model in the supported subset runs unchanged. This one has a
depot, a central compartment, and three parameters with a random effect
each. Compiling it takes milliseconds and no C compiler.

``` r

m <- pmx_ode_model({
  ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
  d/dt(depot)   <- -ka * depot
  d/dt(central) <-  ka * depot - cl / v * central
  cp <- central / v
})
m
#> <pmx_ode_model> linear, 2 compartments on the native tier
#>   compartments: depot, central
#>   parameters:   tka, eta.ka, tcl, eta.cl, tv, eta.v
#>   outputs:      ka, cl, v, cp
```

[`pmx_ode_info()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_info.md)
says what the compiler found: the compartments, the parameters, which
solver the model takes and why, and its tier.

``` r

pmx_ode_info(m)
#> <pmx_ode_info>
#>   class:        linear
#>   solver:       expm (linear in its compartments, with coefficients constant between events (2 of at most 20 compartments))
#>   tier:         native, lowered in 1.0 ms
#>   compartments: depot, central
#>   parameters:   tka, eta.ka, tcl, eta.cl, tv, eta.v
#>   outputs:      ka, cl, v, cp
#>   instructions: pre 12, ivl 15, rhs 13, post 9, err 0
```

## Describe the doses

An **event table** holds NONMEM-style records: doses, and the times to
report. Here it is a week of daily 100 mg doses into the depot, observed
every hour. A data frame with NONMEM’s columns (`ID`, `TIME`, `AMT`,
`CMT`, `EVID`, `MDV` and the rest) works the same way.

``` r

ev <- pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
                      time = seq(0, 168, by = 1))
head(ev, 3)
#>   ID TIME EVID AMT   CMT RATE DUR II ADDL SS
#> 1  1    0    1 100 depot    0   0 24    6  0
#> 2  1    0    0   0 depot    0   0  0    0  0
#> 3  1    1    0   0 depot    0   0  0    0  0
```

## Solve one subject

[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md)
is the deterministic core: every parameter in, the compartments and
outputs at each observation out, with no random draws. Look for `cp`
rising after each dose and settling into a daily cycle.

``` r

p <- c(tka = 0.5, tcl = 1, tv = 3, eta.ka = 0, eta.cl = 0, eta.v = 0)
out <- pmx_ode_solve(m, ev, params = p)
out[c(1, 2, 25, 26, 169), ]
#>     id time        depot   central       ka       cl        v        cp
#> 1    1    0 1.000000e+02  0.000000 1.648721 2.718282 20.08554 0.0000000
#> 2    1    1 1.922956e+01 74.203752 1.648721 2.718282 20.08554 3.6943873
#> 25   1   24 1.000000e+02  4.232420 1.648721 2.718282 20.08554 0.2107198
#> 26   1   25 1.922956e+01 77.900445 1.648721 2.718282 20.08554 3.8784348
#> 169  1  168 6.535319e-16  4.403496 1.648721 2.718282 20.08554 0.2192371
```

The mistake new users make most is to leave a parameter out. The error
names it rather than solving with a default:

``` r

pmx_ode_solve(m, ev, params = c(tka = 0.5, tcl = 1, tv = 3))
#> Error:
#> ! No value for the model's parameters `eta.ka`, `eta.cl`, `eta.v`.
```

## Simulate a population, reproducibly

[`pmx_ode_sim()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_sim.md)
draws each subject’s etas from `omega` and adds residual error to the
outputs `sigma` names. Every draw is made in R with R’s random number
generator, so the same seed gives the same result at any number of
threads.

``` r

old <- pmx_ode_threads(2)
set.seed(1)
sim <- pmx_ode_sim(m, ev, theta = c(tka = 0.5, tcl = 1, tv = 3),
                    omega = diag(c(0.1, 0.1, 0.1)),
                    sigma = list(cp = prop(0.1)), n = 100)
pmx_ode_threads(old)
head(sim[c("id", "time", "cp", "cp_obs", "eta.ka", "eta.cl", "eta.v")])
#>   id time       cp   cp_obs    eta.ka     eta.cl      eta.v
#> 1  1    0 0.000000 0.000000 0.1294642 -0.1961772 -0.1981021
#> 2  1    1 4.710940 4.217564 0.1294642 -0.1961772 -0.1981021
#> 3  1    2 4.834858 5.787972 0.1294642 -0.1961772 -0.1981021
#> 4  1    3 4.332211 4.166014 0.1294642 -0.1961772 -0.1981021
#> 5  1    4 3.799775 4.428313 0.1294642 -0.1961772 -0.1981021
#> 6  1    5 3.320535 3.822670 0.1294642 -0.1961772 -0.1981021
```

The etas are the parameters `omega` names, or, as here, the model’s
parameters named `eta...` in order. Residual error comes from
[`add()`](https://pedrobtz.github.io/zupmxode/reference/add.md),
[`prop()`](https://pedrobtz.github.io/zupmxode/reference/add.md), their
sum `add() + prop()`, or
[`lnorm()`](https://pedrobtz.github.io/zupmxode/reference/add.md); the
observation with error sits next to the prediction as `cp_obs`.

## Evaluate one subject for an estimator

An estimator tries thousands of parameter vectors on the same data. It
prepares the subjects once with
[`pmx_ode_subjects()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_subjects.md),
which splits the model’s parameters into theta and eta, and then asks
[`pmx_ode_subject_eval()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_subject_eval.md)
for one **prepared subject’s** predictions and residual variances at
each vector.

``` r

subs <- pmx_ode_subjects(m, ev)
attr(subs, "eta")
#> [1] "eta.ka" "eta.cl" "eta.v"
r <- pmx_ode_subject_eval(m, subs[[1]], theta = c(tka = 0.5, tcl = 1, tv = 3),
                          eta = c(eta.ka = 0, eta.cl = 0.2, eta.v = 0),
                          sigma = list(cp = prop(0.1)))
head(cbind(time = r$time, pred = r$pred, var = r$var), 3)
#>      time     pred       var
#> [1,]    0 0.000000 0.0000000
#> [2,]    1 3.626342 0.1315036
#> [3,]    2 3.771163 0.1422167
```

The predictions are the doubles
[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md)
gives for the same parameters. The same call exists in C, for estimators
written there, through the header `inst/include/zupmxode.h`.

## Choose the solver and the tier

The solver is the most specific one the model allows (closed form,
matrix exponential, or a general ODE solver), and the tier is native
code unless the platform cannot run it. Both can be overridden; each
override returns a new model and leaves `m` as it was.

``` r

pmx_ode_method(pmx_ode_method(m, "lsoda"))
#> [1] "lsoda"
pmx_ode_tier(pmx_ode_tier(m, "interp"))
#> [1] "interp"
```

## Where next

- [Performance and
  platforms](https://pedrobtz.github.io/zupmxode/articles/performance.html):
  how these calls compare with rxode2 and mrgsolve, and which tier runs
  where.
- [`?pmx_ode_model`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_model.md):
  the subset of rxode2’s language a model may use.
- [`?pmx_ode_events`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_events.md):
  every column of an event table and how it is read.
- The
  [reference](https://pedrobtz.github.io/zupmxode/reference/index.html)
  lists every function.
