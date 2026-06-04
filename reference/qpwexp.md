# Quantile function for piecewise exponential regression

Obtains the quantile function value for piecewise exponential
regression.

## Usage

``` r
qpwexp(p, theta, J, tcut, q = 0, x = 1, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- p:

  The vector of probabilities.

- theta:

  The parameter vector consisting of gamma for log piecewise hazards and
  beta for regression coefficients.

- J:

  The number of time intervals.

- tcut:

  A vector that specifies the endpoints of time intervals for the
  baseline piecewise exponential survival distribution. Must start with
  0, e.g., c(0, 60) breaks the time axis into 2 event intervals:
  \[0, 60) and \[60, Inf). By default, it is set to 0.

- q:

  The number of elements in the vector of covariates (excluding the
  intercept).

- x:

  The vector of covariates (including the intercept).

- lower.tail:

  logical; if TRUE (default), probabilities are the distribution
  function, otherwise, the survival function.

- log.p:

  logical; if TRUE, probabilities p are given as log(p).

## Value

The quantiles t such that P(T \<= t \| X = x) = p.
