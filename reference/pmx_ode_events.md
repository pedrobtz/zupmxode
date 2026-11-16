# Build an event table

Doses and observation times as NONMEM-style records: one dose record per
element of `amt` (recycled with `dose_time`, `cmt`, `rate`, `dur`, `ii`,
`addl` and `ss`), and one observation record per element of `time`.

## Usage

``` r
pmx_ode_events(
  amt = double(),
  cmt = 1L,
  time = double(),
  dose_time = 0,
  ii = 0,
  addl = 0L,
  rate = 0,
  dur = 0,
  ss = 0L,
  id = 1L
)
```

## Arguments

- amt:

  Dose amounts.

- cmt:

  The compartment each dose goes into, by name or number.

- time:

  Observation times.

- dose_time:

  When each dose is given.

- ii, addl:

  Interval and number of additional doses after each one.

- rate, dur:

  Infusion rate or duration; 0 for a bolus. `rate = -1` takes the rate
  from the model's `rate()`, `rate = -2` the duration from its `dur()`.

- ss:

  Steady state: 1 resets and doses to steady state, 2 adds the steady
  state to what is there.

- id:

  The subject.

## Value

A data frame of class `pmx_ode_events` with columns ID, TIME, EVID, AMT,
CMT, RATE, DUR, II, ADDL and SS, sorted by time with doses before
observations at the same time.

## Examples

``` r
pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
               time = seq(0, 168, by = 24))
#>   ID TIME EVID AMT   CMT RATE DUR II ADDL SS
#> 1  1    0    1 100 depot    0   0 24    6  0
#> 2  1    0    0   0 depot    0   0  0    0  0
#> 3  1   24    0   0 depot    0   0  0    0  0
#> 4  1   48    0   0 depot    0   0  0    0  0
#> 5  1   72    0   0 depot    0   0  0    0  0
#> 6  1   96    0   0 depot    0   0  0    0  0
#> 7  1  120    0   0 depot    0   0  0    0  0
#> 8  1  144    0   0 depot    0   0  0    0  0
#> 9  1  168    0   0 depot    0   0  0    0  0
```
