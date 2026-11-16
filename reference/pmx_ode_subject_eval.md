# Evaluate one prepared subject

One subject's predictions and residual variances at its observation
records, at `theta` and `eta`: what an estimator needs to form the
likelihood. The predictions are the doubles
[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md)
gives for the same subject and parameters. This is
`zupmxode_subject_eval()` of the C API (`inst/include/zupmxode.h`), one
R call per subject.

## Usage

``` r
pmx_ode_subject_eval(
  model,
  subject,
  theta = double(),
  eta = double(),
  sigma,
  tol = NULL
)
```

## Arguments

- model:

  A model from
  [`pmx_ode_model()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_model.md).

- subject:

  One element of
  [`pmx_ode_subjects()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_subjects.md)'s
  result.

- theta, eta:

  Named numeric vectors with the parameters the subject was prepared
  with (`attr(subjects, "theta")` and `attr(subjects, "eta")`), in any
  order.

- sigma:

  The output and its error model, a list of one element, as in
  [`pmx_ode_sim()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_sim.md):
  `list(cp = prop(0.1))`.

- tol:

  The general solvers' tolerances, `c(rtol, atol)`; `NULL` for the
  model's.

## Value

A list with `time`, `pred` and `var`, one value each per observation
record.

## Details

The variance follows the error model, `f` the prediction: `add(a)` gives
`a^2`, `prop(b)` gives `(b * f)^2`, `add(a) + prop(b)` their sum, and
`lnorm(s)` gives `s^2`, the variance on the log scale.

## Errors

`pmx_ode_error_argument` for arguments that do not match the subject and
model, or a subject saved and loaded again; `pmx_ode_error_events` and
`pmx_ode_error_solver` as for
[`pmx_ode_solve()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_solve.md).

## Warnings

`pmx_ode_warning_solver` when the solver fails: `pred` and `var` are
`NA` from there. `pmx_ode_warning_steady_state` when a steady state does
not converge.

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
