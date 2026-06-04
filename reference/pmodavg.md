# Distribution function for model averaging of Weibull and log-normal

Obtains the distribution function value for model-averaging of Weibull
and log-normal regression.

## Usage

``` r
pmodavg(t, theta, w1, q = 0, x = 1, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- t:

  The vector of time points.

- theta:

  The parameter vector consisting of the accelerate failure time (AFT)
  regression coefficients and the logrithm of the AFT regression scale
  parameter for the Weibull and log-normal distributions.

- w1:

  The weight for the Weibull component distribution.

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

The probabilities p = P(T \<= t \| X = x).
