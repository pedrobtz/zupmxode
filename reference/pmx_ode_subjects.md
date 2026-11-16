# Prepare subjects for repeated evaluation

Reads an event table once and keeps each subject's records ready for
[`pmx_ode_subject_eval()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_subject_eval.md),
which an estimator calls at every parameter vector it tries. The same
subjects serve the C API in `inst/include/zupmxode.h`.

## Usage

``` r
pmx_ode_subjects(model, data, eta = NULL)
```

## Arguments

- model:

  A model from
  [`pmx_ode_model()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_model.md).

- data:

  An event table, as for
  [`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md).

- eta:

  The names of the parameters that are etas. By default the model's
  parameters whose names start with `eta` that are not covariates.

## Value

A list of class `pmx_ode_subjects`, one prepared subject per ID in the
order they first appear, each a list with `id`, `time` (its observation
times) and the subject itself. The attributes `theta` and `eta` name the
parameters
[`pmx_ode_subject_eval()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_subject_eval.md)
expects. A prepared subject lives for the session: one saved and loaded
again must be prepared again.

## Details

The model's parameters split three ways. Those named in `eta` are the
subject's random effects; columns of `data` named like a parameter are
covariates; the rest are theta, in the order the model uses them. The
split is stored on the result as the attributes `theta` and `eta`.

## Errors

`pmx_ode_error_argument` for an `eta` that is not a parameter of the
model or is a covariate; `pmx_ode_error_events` as for
[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md).

## See also

[`pmx_ode_subject_eval()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_subject_eval.md).

## Examples

``` r
m <- pmx_ode_model({
  cl <- exp(tcl + eta.cl); v <- exp(tv); ka <- exp(tka)
  cp <- linCmt()
})
ev <- pmx_ode_events(amt = 100, cmt = "depot", time = c(1, 4, 12, 24))
subs <- pmx_ode_subjects(m, ev)
attr(subs, "theta")
#> [1] "tcl" "tv"  "tka"
pmx_ode_subject_eval(m, subs[[1]], theta = c(tcl = 0, tv = 3, tka = -0.5),
                     eta = c(eta.cl = 0.1), sigma = list(cp = prop(0.1)))
#> $time
#> [1]  1  4 12 24
#> 
#> $pred
#> [1] 2.196872 3.909813 2.825413 1.461862
#> 
#> $var
#> [1] 0.04826248 0.15286634 0.07982960 0.02137041
#> 
```
