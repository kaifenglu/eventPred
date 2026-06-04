# Piecewise exponential regression

Obtains the maximum likelihood estimates for piecewise exponential
regression.

## Usage

``` r
pwexpreg(time, event, J, tcut, q = 0, x = 1)
```

## Arguments

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

- q:

  The number of columns of the covariates matrix (exluding the
  intercept).

- x:

  The covariates matrix (including the intercept).

## Value

The maximum likelihood estimates and the associated covariance matrix,
AIC and BIC.
