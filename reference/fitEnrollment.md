# Fit enrollment model

Fits a specified enrollment model to the enrollment data.

## Usage

``` r
fitEnrollment(
  df,
  enroll_model = "b-spline",
  nknots = 0,
  accrualTime = 0,
  showplot = TRUE,
  generate_plot = TRUE,
  interactive_plot = TRUE
)
```

## Arguments

- df:

  The subject-level enrollment data, including `trialsdt`, `randdt` and
  `cutoffdt`.

- enroll_model:

  The enrollment model which can be specified as "Poisson",
  "Time-decay", "B-spline", or "Piecewise Poisson". By default, it is
  set to "B-spline".

- nknots:

  The number of inner knots for the B-spline enrollment model. By
  default, it is set to 0.

- accrualTime:

  The accrual time intervals for the piecewise Poisson model. Must start
  with 0, e.g., c(0, 30) breaks the time axis into 2 accrual intervals:
  \[0, 30) and \[30, Inf). By default, it is set to 0.

- showplot:

  A Boolean variable to control whether or not to show the fitted
  enrollment curve. By default, it is set to `TRUE`.

- generate_plot:

  Whether to generate plots.

- interactive_plot:

  Whether to produce interactive plots using plotly or static plots
  using ggplot2.

## Value

A list of results from the model fit including key information such as
the enrollment model, `model`, the estimated model parameters, `theta`,
the covariance matrix, `vtheta`, the Akaike Information Criterion,
`aic`, and the Bayesian Information Criterion, `bic`, as well as the
design matrix `x` for the B-spline enrollment model, and `accrualTime`
for the piecewise Poisson enrollment model.

The fitted enrollment curve is also returned.

## Details

For the time-decay model, the mean function is \$\$\mu(t) =
(\mu/\delta)(t - (1/\delta)(1 - \exp(-\delta t)))\$\$ and the rate
function is \$\$\lambda(t) = (\mu/\delta)(1 - \exp(-\delta t)).\$\$ For
the B-spline model, the daily enrollment rate is \\\lambda(t) =
\exp(B(t)' \theta)\\, where \\B(t)\\ represents the B-spline basis
functions.

## References

Xiaoxi Zhang and Qi Long. Stochastic modeling and prediction for accrual
in clinical trials. Stat in Med. 2010; 29:649-658.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

enroll_fit <- fitEnrollment(
  df = interimData1, enroll_model = "b-spline",
  nknots = 1)
```
