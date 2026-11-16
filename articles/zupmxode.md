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
```

## Solve one subject

[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md)
is the deterministic core: parameters in, predictions out, with no
random draws.

``` r

out <- pmx_ode_solve(m, ev, params = c(tka = 0.5, tcl = 1, tv = 3,
                                        eta.ka = 0, eta.cl = 0,
                                        eta.v = 0))
```

## Simulate a population, reproducibly

Between-subject variability and residual error are drawn in R with R’s
random number generator, and the C engine stays deterministic. The same
seed gives the same result at any thread count.

``` r

pmx_ode_threads(8)
set.seed(1)
sim <- pmx_ode_sim(m, ev, theta = c(tka = 0.5, tcl = 1, tv = 3),
                    omega = diag(c(0.1, 0.1, 0.1)),
                    sigma = list(cp = prop(0.1)), n = 1000)
```

## Choose the solver and the tier

The solver is the most specific one the model allows (closed form,
matrix exponential, or a general ODE solver), and the tier is native
code unless the platform cannot run it. `pmx_ode_info()` says which was
chosen, and both can be overridden.

``` r

pmx_ode_method(m, "lsoda")
pmx_ode_tier(m, "interp")
```

## Where next

Every function is listed in the
[reference](https://pedrobtz.github.io/zupmxode/reference/index.html).
