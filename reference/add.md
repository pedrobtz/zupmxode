# Residual error models

The four error models
[`pmx_ode_sim()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_sim.md)
applies to a prediction `f`, each with `e` a standard normal draw:

## Usage

``` r
add(sd)

prop(sd)

lnorm(sd)
```

## Arguments

- sd:

  The standard deviation: a single non-negative finite number.

## Value

An error model, for the `sigma` argument of
[`pmx_ode_sim()`](https://pedrobtz.github.io/zupmxode/reference/pmx_ode_sim.md).

## Details

- `add(a)`: additive, `f + a * e`;

- `prop(b)`: proportional, `f + f * b * e`;

- `add(a) + prop(b)`: combined, `f + sqrt(a^2 + b^2 * f^2) * e`, the
  variances added;

- `lnorm(s)`: exponential, `f * exp(s * e)`.

The spellings are rxode2's.

## Examples

``` r
prop(0.1)
#> proportional error, sd 0.1
add(0.5) + prop(0.1)
#> combined error, additive sd 0.5, proportional sd 0.1
```
