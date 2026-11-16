# zupmxode

zupmxode simulates population pharmacokinetic and pharmacodynamic
models, written in [rxode2](https://nlmixr2.github.io/rxode2/)’s model
language, with no C compiler: for locked-down and validated machines,
where rxode2 and mrgsolve, which need one, cannot run.

- A model becomes native code in about 10 ms through
  [zumir](https://github.com/pedrobtz/zumir)’s
  [MIR](https://github.com/vnmakarov/mir); rxode2 and mrgsolve take one
  to two seconds. MIR’s interpreter, the fallback, gives the same
  doubles 1.4 to 1.6 times slower; every R CMD check leg runs both.
- Doses are NONMEM-style records, read as PREDPP reads them.
- Closed forms, the matrix exponential, LSODA or Dormand-Prince, with
  subjects on parallel threads and the same result at any thread count.

## Installation

zupmxode is not on CRAN yet. Install the development version with:

``` r

# install.packages("pak")
pak::pak("pedrobtz/zupmxode")
```

## Usage

Five subjects of a one-compartment model, each with its own clearance,
after 100 mg into the depot at time 0; one row per observation:

``` r

library(zupmxode)

m <- pmx_ode_model({
  ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
  cp <- linCmt()
})
ev <- data.frame(TIME = c(0, 1, 4, 12), AMT = c(100, 0, 0, 0), CMT = "depot")

set.seed(1)
sim <- pmx_ode_sim(m, ev, theta = c(tka = 0.5, tcl = 1, tv = 3),
                   omega = c(eta.cl = 0.1), n = 5)
head(sim[c("id", "time", "cp", "eta.cl")], 8)
#>   id time        cp      eta.cl
#> 1  1    1 3.7507426 -0.19810209
#> 2  1    4 3.4167467 -0.19810209
#> 3  1   12 1.4087652 -0.19810209
#> 4  2    1 3.6758616  0.05807312
#> 5  2    4 3.0649859  0.05807312
#> 6  2   12 0.9753628  0.05807312
#> 7  3    1 3.7673999 -0.26424897
#> 8  3    4 3.4992927 -0.26424897
```

zupmxode also adds residual error, solves without random draws, and lets
an estimator evaluate one subject at a time from R or C; [Get
started](https://pedrobtz.github.io/zupmxode/articles/zupmxode.html)
shows how.

## Learn more

- [Get
  started](https://pedrobtz.github.io/zupmxode/articles/zupmxode.html):
  compile a model, solve one subject, simulate a population, choose the
  solver and the tier.
- [Performance and
  platforms](https://pedrobtz.github.io/zupmxode/articles/performance.html):
  the benchmarks against rxode2 and mrgsolve, and the tier each platform
  runs.
- [Reference](https://pedrobtz.github.io/zupmxode/reference/index.html):
  every function, and the C API in `inst/include/zupmxode.h`.
