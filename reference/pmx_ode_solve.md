# Solve a model for one or more subjects

Runs `model` over the event table `events` with the parameters `params`,
and returns the compartments and the model's outputs at every
observation record. Nothing is drawn at random: this is the
deterministic core that the population simulation builds on.

## Usage

``` r
pmx_ode_solve(model, events, params = double())
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

- params:

  The model's parameters, a named numeric vector. A covariate column of
  `events` takes the place of its parameter.

## Value

A data frame with `id`, `time`, one column per compartment and one per
model output, one row per observation record.

## Errors

`pmx_ode_error_events` for a record that breaks a rule of the event
table, naming the row; `pmx_ode_error_argument` for a parameter with no
value; `pmx_ode_error_solver` for a closed-form model whose rate
constants are not positive.

## Warnings

`pmx_ode_warning_solver` when the solver fails for a subject: that
subject's outputs are `NA` from where it failed, and the others are
solved (design §8). `pmx_ode_warning_steady_state` when a steady state
of a model with no closed form does not converge.

## Examples

``` r
m <- pmx_ode_model({
  cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
  cp <- linCmt()
})
ev <- pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
                     time = seq(0, 168, by = 12))
head(pmx_ode_solve(m, ev, params = c(tcl = 0, tv = 3, tka = -0.7)))
#>   id time       depot  central cl        v        ka       cp
#> 1  1    0 100.0000000  0.00000  1 20.08554 0.4965853 0.000000
#> 2  1   12   0.2582432 60.86565  1 20.08554 0.4965853 3.030322
#> 3  1   24 100.0006669 33.64642  1 20.08554 0.4965853 1.675157
#> 4  1   36   0.2582449 79.37885  1 20.08554 0.4965853 3.952040
#> 5  1   48 100.0006669 43.83267  1 20.08554 0.4965853 2.182300
#> 6  1   60   0.2582449 84.98348  1 20.08554 0.4965853 4.231079
```
