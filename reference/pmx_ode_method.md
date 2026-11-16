# The solver of a model

Every model has a default solver, the most specific its class allows
(design §4.3, D6): the closed form for a `linCmt()` model, the matrix
exponential for a linear one of up to 20 compartments, LSODA otherwise.
`pmx_ode_method(model)` says which a model uses;
`pmx_ode_method(model, method)` returns a copy that uses `method`.
[`pmx_ode_info()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_info.md)
says why.

## Usage

``` r
pmx_ode_method(model, method)
```

## Arguments

- model:

  A model.

- method:

  `"auto"` (the default choice), `"closed"` (closed-form models only),
  `"expm"` (closed-form or linear models of up to 20 compartments),
  `"lsoda"` or `"dopri5"` (Dormand-Prince 5(4)).

## Value

The method in use, or, given `method`, the model with it.

## Errors

`pmx_ode_error_argument` for a method the model's class cannot take.

## Examples

``` r
m <- pmx_ode_model({
  d/dt(central) <- -k * central
})
pmx_ode_method(m)
#> [1] "expm"
pmx_ode_method(pmx_ode_method(m, "lsoda"))
#> [1] "lsoda"
```
