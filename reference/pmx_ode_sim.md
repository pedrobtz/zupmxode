# Simulate a population

Draws each subject's etas from a multivariate normal with covariance
`omega`, solves every subject at `theta` and its etas with
[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md),
and adds residual error to the outputs `sigma` names. Every draw is made
in R with R's random number generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) makes a simulation
reproducible, at any thread count.

## Usage

``` r
pmx_ode_sim(
  model,
  events,
  theta = double(),
  omega = NULL,
  sigma = NULL,
  n = NULL
)
```

## Arguments

- model:

  A model from
  [`pmx_ode_model()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_model.md).

- events:

  An event table: from
  [`pmx_ode_events()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_events.md),
  or a NONMEM-style data frame with columns `ID`, `TIME`, `EVID`, `AMT`,
  `CMT`, `RATE`, `DUR`, `II`, `ADDL` and `SS` (any case; all but `TIME`
  optional). With no `EVID`, a record with an `AMT` above 0 is a dose,
  one with `MDV` 1 marks a time (`EVID` 2), and the rest are
  observations; given `EVID`, `MDV` is not read. `CMT` is a
  compartment's name or number. Columns named like a model parameter are
  covariates, each record carrying its value; other columns, such as
  `DV`, are ignored.

- theta:

  The model's fixed parameters, a named numeric vector.

- omega:

  The covariance of the etas: a symmetric positive semi-definite matrix,
  or a named vector of variances for a diagonal one. Its names, or with
  none the model's parameters whose names start with `eta`, in the order
  the model uses them, say which parameters are etas. `NULL` draws none.

- sigma:

  Residual error: a named list of error models, from
  [`add()`](https://pedrobtz.github.io/zupmxode/reference/add.md),
  [`prop()`](https://pedrobtz.github.io/zupmxode/reference/add.md) or
  [`lnorm()`](https://pedrobtz.github.io/zupmxode/reference/add.md),
  each name a compartment or an output of the model. `NULL` adds none.

- n:

  The number of subjects. A table of one subject is copied `n` times,
  the copies numbered 1 to `n`; a table of several subjects is simulated
  once per subject, and `n`, if given, must be their number.

## Value

A data frame with
[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md)'s
columns, then one column per eta holding the subject's draw, then, for
each entry of `sigma`, the output with residual error, named after it
with `_obs` appended.

## Details

The draws are made in a fixed order: first the etas, as
`MASS::mvrnorm(n, 0, omega)` would draw them, one row per subject; then,
for each entry of `sigma` in turn, one standard normal per observation
record. A subject the solver fails for still takes its draws, so the
others' are unchanged.

Four residual error models are built in:
[`add()`](https://pedrobtz.github.io/zupmxode/reference/add.md),
[`prop()`](https://pedrobtz.github.io/zupmxode/reference/add.md), their
sum `add(a) + prop(b)` (the combined model) and
[`lnorm()`](https://pedrobtz.github.io/zupmxode/reference/add.md). Any
other is a function applied to the predictions this returns.

## Errors

`pmx_ode_error_argument` for an argument that is not what is described
here; the errors of
[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md).

## See also

[`add()`](https://pedrobtz.github.io/zupmxode/reference/add.md) for the
error models.

## Examples

``` r
m <- pmx_ode_model({
  ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
  cp <- linCmt()
})
ev <- pmx_ode_events(amt = 100, cmt = "depot", time = c(1, 2, 4, 8, 12, 24))
set.seed(1)
sim <- pmx_ode_sim(m, ev, theta = c(tka = 0.5, tcl = 1, tv = 3),
                   omega = diag(c(0.1, 0.1)), sigma = list(cp = prop(0.1)),
                   n = 10)
head(sim)
#>   id time        depot  central       ka       cl        v        cp     eta.cl
#> 1  1    1 1.922956e+01 75.74679 1.648721 1.685282 16.47589 4.5974332 -0.4780672
#> 2  1    2 3.697762e+00 82.94768 1.648721 1.685282 16.47589 5.0344898 -0.4780672
#> 3  1    4 1.367344e-01 70.66900 1.648721 1.685282 16.47589 4.2892382 -0.4780672
#> 4  1    8 1.869630e-04 47.03596 1.648721 1.685282 16.47589 2.8548363 -0.4780672
#> 5  1   12 2.556427e-07 31.24207 1.648721 1.685282 16.47589 1.8962299 -0.4780672
#> 6  1   24 6.535319e-16  9.15511 1.648721 1.685282 16.47589 0.5556672 -0.4780672
#>        eta.v    cp_obs
#> 1 -0.1981021 5.0199269
#> 2 -0.1981021 5.4282556
#> 3 -0.1981021 4.3212209
#> 4 -0.1981021 2.2869090
#> 5 -0.1981021 2.0137631
#> 6 -0.1981021 0.5525483
```
