# What a model is and how it is solved

What a model is and how it is solved

## Usage

``` r
pmx_ode_info(model)
```

## Arguments

- model:

  A model.

## Value

A list of class `pmx_ode_info`: the compartments, parameters and
outputs; the class the compiler found; the solver and why it was chosen;
the execution tier; the seconds lowering took; and the IR's instructions
per block.

## Examples

``` r
m <- pmx_ode_model({
  cl <- exp(tcl); v <- exp(tv)
  cp <- linCmt()
})
pmx_ode_info(m)
#> <pmx_ode_info>
#>   class:        closed-form
#>   solver:       closed (a linCmt() model with no other derivative)
#>   tier:         native, lowered in 0.0 ms
#>   compartments: central
#>   parameters:   tcl, tv
#>   outputs:      cl, v, cp
#>   instructions: pre 7, ivl 6, rhs 7, post 7, err 0
```
