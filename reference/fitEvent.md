# Fit time-to-event model

Fits a specified time-to-event model to the event data.

## Usage

``` r
fitEvent(
  df,
  event_model = "model averaging",
  piecewiseSurvivalTime = 0,
  k = 0,
  scale = "hazard",
  m = 5,
  showplot = TRUE,
  by_treatment = FALSE,
  covariates = NULL,
  generate_plot = TRUE,
  interactive_plot = TRUE
)
```

## Arguments

- df:

  The subject-level event data, including `time` and `event`. The data
  should also include `treatment` coded as 1, 2, and so on, and
  `treatment_description` for fitting the event model by treatment.

- event_model:

  The event model used to analyze the event data which can be set to one
  of the following options: "exponential", "Weibull", "log-logistic",
  "log-normal", "piecewise exponential", "model averaging", "spline", or
  "cox". The model averaging uses the `exp(-bic/2)` weighting and
  combines Weibull and log-normal models. The spline model of Royston
  and Parmar (2002) assumes that a transformation of the survival
  function is modeled as a natural cubic spline function of log time. By
  default, it is set to "model averaging".

- piecewiseSurvivalTime:

  A vector that specifies the time intervals for the piecewise
  exponential survival distribution. Must start with 0, e.g., c(0, 60)
  breaks the time axis into 2 event intervals: \[0, 60) and \[60, Inf).
  By default, it is set to 0.

- k:

  The number of inner knots of the spline. The default `k=0` gives a
  Weibull, log-logistic or log-normal model, if `scale` is "hazard",
  "odds", or "normal", respectively. The knots are chosen as
  equally-spaced quantiles of the log uncensored survival times. The
  boundary knots are chosen as the minimum and maximum log uncensored
  survival times.

- scale:

  The scale of the spline. The default is "hazard", in which case the
  log cumulative hazard is modeled as a spline function. If
  `scale = "odds"`, the log cumulative odds is modeled as a spline
  function. If `scale = "normal"`, `-qnorm(S(t))` is modeled as a spline
  function.

- m:

  The number of event time intervals to extrapolate the hazard function
  beyond the last observed event time when `event_model = "cox"`.

- showplot:

  A Boolean variable to control whether or not to show the fitted
  time-to-event survival curve. By default, it is set to `TRUE`.

- by_treatment:

  A Boolean variable to control whether or not to fit the time-to-event
  data by treatment group. By default, it is set to `FALSE`.

- covariates:

  The names of baseline covariates from the input data frame to include
  in the event model, e.g., c("age", "sex"). Factor variables need to be
  declared in the input data frame.

- generate_plot:

  Whether to generate plots.

- interactive_plot:

  Whether to produce interactive plots using plotly or static plots
  using ggplot2.

## Value

A list of results from the model fit including key information such as
the event model, `model`, the estimated model parameters, `theta`, the
covariance matrix, `vtheta`, as well as the Akaike Information
Criterion, `aic`, and Bayesian Information Criterion, `bic`.

If the piecewise exponential model is used, the location of knots used
in the model, `piecewiseSurvivalTime`, will be included in the list of
results.

If the model averaging option is chosen, the weight assigned to the
Weibull component is indicated by the `w1` variable.

If the spline option is chosen, the `knots` and `scale` will be included
in the list of results.

If the cox option is chosen, the list of results will include `model`,
`theta`, `vtheta`, `aic`, `bic`, and `piecewiseSurvivalTime`. Here
\$\$\theta = (\log(\lambda_1), \ldots, \log(\lambda_M), \beta^T)^T,\$\$
\\M\\ denotes the number of distinct observed event times, \\t_1 \<
\cdots \< t_M\\, \\\lambda_j\\ denotes the estimated baseline hazard
rate in the \\j\\th event time interval, \\(t\_{j-1}, t_j\]\\, and
\\\beta\\ represents the regression coefficients (log hazard ratios)
from the Cox model. For a fair comparison, the estimation of baseline
hazards is incorporated into the `aic` and `bic` values. In addition,
\\\mbox{piecewiseSurvivalTime} = (0, t_1, \ldots, t_M)\\. To extend the
survival curve beyond the last observed event time, a weighted average
of the hazard rates from the final `m` event time intervals is used. The
weights are proportional to the lengths of those intervals, i.e.,
\$\$\lambda\_{M+1} = \sum\_{j=M-m+1}^{M} w_j \lambda_j,\$\$ where \\w_j
= (t_j - t\_{j-1})/(t_M - t\_{M-m})\\ for \\j=M-m+1,\ldots,M\\.

When fitting the event model by treatment, the outcome is presented as a
list of lists, where each list element corresponds to a specific
treatment group.

The fitted time-to-event survival curve is also returned.

## References

Patrick Royston and Mahesh K. B. Parmar. Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data, with application to prognostic modelling and estimation of
treatment effects. Stat in Med. 2002; 21:2175-2197.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

event_fit <- fitEvent(
  df = interimData2,
  event_model = "piecewise exponential",
  piecewiseSurvivalTime = c(0, 180))
```
