# Profile log likelihood for piecewise exponential regression

Obtains the profile log likelihood value for piecewise exponential
regression.

## Usage

``` r
pllik_pwexp(beta, time, event, J, tcut, x)
```

## Arguments

- beta:

  The regression coefficients with respect to the covariates.

- time:

  The survival time.

- event:

  The event indicator.

- J:

  The number of time intervals.

- tcut:

  A vector that specifies the endpoints of time intervals for the
  baseline piecewise exponential survival distribution. Must start with
  0, e.g., c(0, 60) breaks the time axis into 2 event intervals:
  \[0, 60) and \[60, Inf). By default, it is set to 0.

- x:

  The covariates matrix (including the intercept).

## Value

The profile log likelihood value for piecewise exponential regression.
