# Enrollment and event prediction

Performs enrollment and event prediction by utilizing observed data and
specified enrollment and event models.

## Usage

``` r
getPrediction(
  df = NULL,
  to_predict = "enrollment and event",
  target_n = NA,
  target_d = NA,
  enroll_model = "b-spline",
  nknots = 0,
  lags = 30,
  accrualTime = 0,
  enroll_prior = NULL,
  event_model = "model averaging",
  piecewiseSurvivalTime = 0,
  k = 0,
  scale = "hazard",
  m = 5,
  event_prior = NULL,
  dropout_model = "exponential",
  piecewiseDropoutTime = 0,
  k_dropout = 0,
  scale_dropout = "hazard",
  m_dropout = 5,
  dropout_prior = NULL,
  fixedFollowup = FALSE,
  followupTime = 365,
  pilevel = 0.9,
  nyears = 4,
  target_t = NA,
  nreps = 500,
  showEnrollment = TRUE,
  showEvent = TRUE,
  showDropout = FALSE,
  showOngoing = FALSE,
  showsummary = TRUE,
  showplot = TRUE,
  by_treatment = FALSE,
  ngroups = 1,
  alloc = NULL,
  treatment_label = NULL,
  covariates_event = NULL,
  event_prior_with_covariates = NULL,
  covariates_dropout = NULL,
  dropout_prior_with_covariates = NULL,
  fix_parameter = FALSE,
  generate_plot = TRUE,
  interactive_plot = TRUE,
  nthreads = 0
)
```

## Arguments

- df:

  The subject-level enrollment and event data, including `trialsdt`,
  `usubjid`, `randdt`, and `cutoffdt` for enrollment prediction, and,
  additionally, `time`, `event`, and `dropout` for event prediction. The
  data should also include `treatment` coded as 1, 2, and so on, and
  `treatment_description` for enrollment and event prediction by
  treatment. By default, it is set to `NULL` for enrollment and event
  prediction at the design stage.

- to_predict:

  Specifies what to predict: "enrollment only", "event only", or
  "enrollment and event". By default, it is set to "enrollment and
  event".

- target_n:

  The target number of subjects to enroll in the study.

- target_d:

  The target number of events to reach in the study.

- enroll_model:

  The enrollment model which can be specified as "Poisson",
  "Time-decay", "B-spline", or "Piecewise Poisson". By default, it is
  set to "B-spline".

- nknots:

  The number of inner knots for the B-spline enrollment model. By
  default, it is set to 0.

- lags:

  The day lags to compute the average enrollment rate to carry forward
  for the B-spline enrollment model. By default, it is set to 30.

- accrualTime:

  The accrual time intervals for the piecewise Poisson model. Must start
  with 0, e.g., c(0, 30) breaks the time axis into 2 accrual intervals:
  \[0, 30) and \[30, Inf). By default, it is set to 0.

- enroll_prior:

  The prior of enrollment model parameters.

- event_model:

  The event model used to analyze the event data which can be set to one
  of the following options: "exponential", "Weibull", "log-logistic",
  "log-normal", "piecewise exponential", "model averaging", "spline", or
  "cox". The model averaging uses the `exp(-bic/2)` weighting and
  combines Weibull and log-normal models. By default, it is set to
  "model averaging".

- piecewiseSurvivalTime:

  A vector that specifies the time intervals for the piecewise
  exponential survival distribution. Must start with 0, e.g., c(0, 60)
  breaks the time axis into 2 event intervals: \[0, 60) and \[60, Inf).
  By default, it is set to 0.

- k:

  The number of inner knots of the spline event model of Royston and
  Parmar (2002). The default `k=0` gives a Weibull, log-logistic or
  log-normal model, if `scale` is "hazard", "odds", or "normal",
  respectively. The knots are chosen as equally-spaced quantiles of the
  log uncensored survival times. The boundary knots are chosen as the
  minimum and maximum log uncensored survival times.

- scale:

  If "hazard", the log cumulative hazard is modeled as a spline
  function. If "odds", the log cumulative odds is modeled as a spline
  function. If "normal", -qnorm(S(t)) is modeled as a spline function.

- m:

  The number of event time intervals to extrapolate the hazard function
  beyond the last observed event time.

- event_prior:

  The prior of event model parameters.

- dropout_model:

  The dropout model used to analyze the dropout data which can be set to
  one of the following options: "none", "exponential", "Weibull",
  "log-logistic", "log-normal", "piecewise exponential", "model
  averaging", "spline", or "cox". The model averaging uses the
  `exp(-bic/2)` weighting and combines Weibull and log-normal models. By
  default, it is set to "exponential".

- piecewiseDropoutTime:

  A vector that specifies the time intervals for the piecewise
  exponential dropout distribution. Must start with 0, e.g., c(0, 60)
  breaks the time axis into 2 event intervals: \[0, 60) and \[60, Inf).
  By default, it is set to 0.

- k_dropout:

  The number of inner knots of the spline dropout model of Royston and
  Parmar (2002). The default `k_dropout=0` gives a Weibull, log-logistic
  or log-normal model, if `scale_dropout` is "hazard", "odds", or
  "normal", respectively. The knots are chosen as equally-spaced
  quantiles of the log uncensored survival times. The boundary knots are
  chosen as the minimum and maximum log uncensored survival times.

- scale_dropout:

  If "hazard", the log cumulative hazard for dropout is modeled as a
  spline function. If "odds", the log cumulative odds is modeled as a
  spline function. If "normal", -qnorm(S(t)) is modeled as a spline
  function.

- m_dropout:

  The number of dropout time intervals to extrapolate the hazard
  function beyond the last observed dropout time.

- dropout_prior:

  The prior of dropout model parameters.

- fixedFollowup:

  A Boolean variable indicating whether a fixed follow-up design is
  used. By default, it is set to `FALSE` for a variable follow-up
  design.

- followupTime:

  The follow-up time for a fixed follow-up design, in days. By default,
  it is set to 365.

- pilevel:

  The prediction interval level. By default, it is set to 0.90.

- nyears:

  The number of years after the data cut for prediction. By default, it
  is set to 4.

- target_t:

  The target number of days after the data cutoff used to predict both
  the number of events and the probability of achieving the target event
  count.

- nreps:

  The number of replications for simulation. By default, it is set to
  500.

- showEnrollment:

  A Boolean variable to control whether or not to show the number of
  enrolled subjects. By default, it is set to `TRUE`.

- showEvent:

  A Boolean variable to control whether or not to show the number of
  events. By default, it is set to `TRUE`.

- showDropout:

  A Boolean variable to control whether or not to show the number of
  dropouts. By default, it is set to `FALSE`.

- showOngoing:

  A Boolean variable to control whether or not to show the number of
  ongoing subjects. By default, it is set to `FALSE`.

- showsummary:

  A Boolean variable to control whether or not to show the prediction
  summary. By default, it is set to `TRUE`.

- showplot:

  A Boolean variable to control whether or not to show the plots. By
  default, it is set to `TRUE`.

- by_treatment:

  A Boolean variable to control whether or not to predict by treatment
  group. By default, it is set to `FALSE`.

- ngroups:

  The number of treatment groups for enrollment prediction at the design
  stage. By default, it is set to 1. It is replaced with the actual
  number of treatment groups in the observed data if `df` is not `NULL`.

- alloc:

  The treatment allocation in a randomization block. By default, it is
  set to `NULL`, which yields equal allocation among the treatment
  groups.

- treatment_label:

  The treatment labels for treatments in a randomization block for
  design stage prediction. It is replaced with the treatment_description
  in the observed data if `df` is not `NULL`.

- covariates_event:

  The names of baseline covariates from the input data frame to include
  in the event model, e.g., c("age", "sex"). Factor variables need to be
  declared in the input data frame.

- event_prior_with_covariates:

  The prior of event model parameters in the presence of covariates.

- covariates_dropout:

  The names of baseline covariates from the input data frame to include
  in the dropout model, e.g., c("age", "sex"). Factor variables need to
  be declared in the input data frame.

- dropout_prior_with_covariates:

  The prior of dropout model parameters in the presence of covariates.

- fix_parameter:

  Whether to fix parameters at the maximum likelihood estimates when
  generating new data for prediction. Defaults to `FALSE`, in which
  case, parameters will be drawn from their approximate posterior
  distribution.

- generate_plot:

  Whether to generate plots.

- interactive_plot:

  Whether to produce interactive plots using plotly or static plots
  using ggplot2.

- nthreads:

  Integer number of threads to use for \`data.table' (0 means the
  default data.table behavior).

## Value

A list containing model-fit objects and prediction objects.

The model-fit objects summarize either:

- the fitted models based on the observed data, or

- the posterior distribution of the model parameters when prior
  information is supplied.

The prediction objects may include:

- simulated enrollment data for future subjects, and

- simulated event data for both ongoing subjects and future subjects.

At the design stage, all predictions are based solely on prior
information. In that case, the output includes `enroll_prior`,
`event_prior`, and `dropout_prior`.

At the analysis stage, predictions are based on:

- the observed-data likelihood when no prior is provided, or

- the posterior distribution when prior information is provided.

When prior information is incorporated, the parameter vector `theta` in
`enroll_post`, `event_post`, `event_post_with_covariates`,
`dropout_post`, and `dropout_post_with_covariates` represents a weighted
average of the prior mean and the maximum likelihood estimate. The
corresponding variance-covariance matrix `vtheta` is the inverse of the
total information matrix, where the total information is the sum of:

- the information from the prior distribution, and

- the information from the observed-data likelihood.

In addition to the model-fit objects, the output also includes the
analysis stage at which prediction is performed, the prediction target,
and the enrollment and event prediction results when applicable.

## Details

For the time-decay model, the mean function is \\\mu(t) =
(\mu/\delta)(t - (1/\delta)(1 - \exp(-\delta t)))\\ and the rate
function is \\\lambda(t) = (\mu/\delta)(1 - \exp(-\delta t))\\. For the
B-spline model, the daily enrollment rate is approximated as
\\\lambda(t) = \exp(B(t)' \theta)\\, where `B(t)` represents the
B-spline basis functions.

The `enroll_prior` variable should be a list that includes `model` to
specify the enrollment model (poisson, time-decay, or piecewise
poisson), `theta` and `vtheta` to indicate the parameter values and the
covariance matrix. One can use a very small value of `vtheta` to fix the
parameter values. For the piecewise Poisson enrollment model, the list
should also include `accrualTime`. It should be noted that the B-spline
model is not appropriate for use as prior.

For event prediction by treatment with prior information, the
`event_prior` (`dropout_prior`) variable should be a list with one
element per treatment. For each treatment, the element should include
`model` to specify the event (dropout) model (exponential, weibull,
log-logistic, log-normal, or piecewise exponential), and `theta` and
`vtheta` to indicate the parameter values and the covariance matrix. For
the piecewise exponential event (dropout) model, the list should also
include `piecewiseSurvivalTime` (`piecewiseDropoutTime`) to indicate the
location of knots. It should be noted that the model averaging, spline,
and cox options are not appropriate for use as prior.

If the event prediction is not by treatment while the prior information
is given by treatment, then each element of `event_prior`
(`dropout_prior`) should also include `w` to specify the weight of the
treatment in a randomization block. If the prediction is not by
treatment and the prior is given for the overall study, then
`event_prior` (`dropout_prior`) is a flat list with `model`, `theta`,
and `vtheta`. For the piecewise exponential event (dropout) model, it
should also include `piecewiseSurvivalTime` (`piecewiseDropoutTime`) to
indicate the location of knots.

For analysis-stage enrollment and event prediction, the `enroll_prior`,
`event_prior`, and `dropout_prior` are either set to `NULL` to use the
observed data only, or specify the prior distribution of model
parameters to be combined with observed data likelihood for enhanced
modeling flexibility.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Event prediction after enrollment completion
set.seed(3000)

pred <- getPrediction(
  df = interimData2,
  to_predict = "event only",
  target_d = 200,
  event_model = "weibull",
  dropout_model = "exponential",
  pilevel = 0.90,
  nreps = 100,
  nthreads = 1)
#> Time from cutoff until 200 events: 130 days 
#>  Median prediction date: 2021-02-27 
#>  Prediction interval: 2021-01-08, 2021-05-18 
```
