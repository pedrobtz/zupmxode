# Compile a model

Reads a model written in the subset of rxode2's model language that
zupmxode runs (design §6.1), and compiles it to native code through MIR,
with no C compiler. The model is R syntax, read by R's parser:
assignments with `<-` or `=`, `d/dt(x)` for a derivative, `x(0)` for an
initial amount, `alag(x)`, `f(x)`, `rate(x)` and `dur(x)` for dosing
into `x`, `if`/`else`, comparison and logical operators, `linCmt()`, and
the functions `exp`, `log`, `sqrt`, `abs`, `floor`, `ceiling`, `sin`,
`cos`, `tan`, `gamma`, `lgamma`, `logit`, `expit`, `probit`,
`probitInv`, `min` and `max`. Every name the model reads but never
assigns is a parameter, given at solve time.

## Usage

``` r
pmx_ode_model(model, tier = "auto")
```

## Arguments

- model:

  The model: a braced expression, a function whose body is the model, or
  its text.

- tier:

  The tier to run on: `"auto"`, `"native"` or `"interp"`; see
  [`pmx_ode_tier()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_tier.md).

## Value

A model of class `pmx_ode_model`.

## Errors

`pmx_ode_error_compile` for a construct outside the subset, naming it,
its position, and the nearest supported alternative.

## Examples

``` r
m <- pmx_ode_model({
  ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
  d/dt(depot) <- -ka * depot
  d/dt(central) <- ka * depot - cl / v * central
  cp <- central / v
})
m
#> <pmx_ode_model> linear, 2 compartments on the native tier
#>   compartments: depot, central
#>   parameters:   tka, tcl, tv
#>   outputs:      ka, cl, v, cp
```
