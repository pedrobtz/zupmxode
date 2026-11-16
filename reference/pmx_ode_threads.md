# The number of threads a solve uses

Subjects are independent and are solved in parallel with OpenMP (design
§4.4, D8); the result is the same, bit for bit, at any thread count.
Without OpenMP (as on some macOS builds of R) everything runs on one
thread whatever is set here.

## Usage

``` r
pmx_ode_threads(n)
```

## Arguments

- n:

  The number of threads, a positive whole number; leave it out to read
  the current setting.

## Value

The current setting, or, given `n`, the one it replaces, invisibly.

## Examples

``` r
old <- pmx_ode_threads(2)
pmx_ode_threads()
#> [1] 2
pmx_ode_threads(old)
```
