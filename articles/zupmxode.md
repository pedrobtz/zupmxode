# Get started with zupmxode

``` r

library(zupmxode)
```

rxode2 and mrgsolve compile every model with the system C compiler,
which locked-down and validated machines often lack. zupmxode runs the
same models with nothing but R. This page compiles a one-compartment
model with first-order absorption, and gives it a week of daily 100 mg
doses:

``` r

ev <- pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
                      time = seq(0, 168, by = 1))
```

## Compile a model

The model is written in rxode2’s syntax and read by R’s own parser, so
existing rxode2 models run unchanged. It is compiled to native code in
milliseconds, without a compiler.

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

## Solve one subject

[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md)
is the deterministic core: parameters in, predictions out, with no
random draws.

``` r

out <- pmx_ode_solve(m, ev, params = c(tka = 0.5, tcl = 1, tv = 3,
                                        eta.ka = 0, eta.cl = 0,
                                        eta.v = 0))
head(out)
#>   id time        depot  central       ka       cl        v       cp
#> 1  1    0 100.00000000  0.00000 1.648721 2.718282 20.08554 0.000000
#> 2  1    1  19.22956455 74.20375 1.648721 2.718282 20.08554 3.694387
#> 3  1    2   3.69776153 79.08032 1.648721 2.718282 20.08554 3.937177
#> 4  1    3   0.71106344 71.81445 1.648721 2.718282 20.08554 3.575431
#> 5  1    4   0.13673440 63.25203 1.648721 2.718282 20.08554 3.149133
#> 6  1    5   0.02629343 55.34724 1.648721 2.718282 20.08554 2.755577
```

## Simulate a population, reproducibly

Between-subject variability and residual error are drawn in R with R’s
random number generator, and the C engine stays deterministic. The same
seed gives the same result at any thread count.

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

## Choose the solver and the tier

The solver is the most specific one the model allows (closed form,
matrix exponential, or a general ODE solver), and the tier is native
code unless the platform cannot run it.
[`pmx_ode_info()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_info.md)
says which was chosen, and both can be overridden.

``` r

pmx_ode_method(pmx_ode_method(m, "lsoda"))
#> [1] "lsoda"
pmx_ode_tier(pmx_ode_tier(m, "interp"))
#> [1] "interp"
```

## Where next

Every function is listed in the
[reference](https://pedrobtz.github.io/zupmxode/reference/index.html).
