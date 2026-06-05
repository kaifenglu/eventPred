# Predict event

Utilizes pre-fitted time-to-event and time-to-dropout models to generate
event and dropout times for ongoing subjects and new subjects. It also
provides a prediction interval for the expected time to reach the target
number of events.

## Usage

``` r
predictEvent(
  df = NULL,
  target_d = NA,
  newSubjects = NULL,
  event_fit = NULL,
  m = 5,
  dropout_fit = NULL,
  m_dropout = 5,
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
  covariates_event = NULL,
  event_fit_with_covariates = NULL,
  covariates_dropout = NULL,
  dropout_fit_with_covariates = NULL,
  fix_parameter = FALSE,
  generate_plot = TRUE,
  interactive_plot = TRUE,
  nthreads = 0
)
```

## Arguments

- df:

  The subject-level enrollment and event data, including `trialsdt`,
  `usubjid`, `randdt`, `cutoffdt`, `time`, `event`, and `dropout`. The
  data should also include `treatment` coded as 1, 2, and so on, and
  `treatment_description` for by-treatment prediction. By default, it is
  set to `NULL` for event prediction at the design stage.

- target_d:

  The target number of events to reach in the study.

- newSubjects:

  The enrollment data for new subjects including `draw` and
  `arrivalTime`. The data should also include `treatment` for prediction
  by treatment. By default, it is set to `NULL`, indicating the
  completion of subject enrollment.

- event_fit:

  The pre-fitted event model used to generate predictions.

- m:

  The number of event time intervals to extrapolate the hazard function
  beyond the last observed event time.

- dropout_fit:

  The pre-fitted dropout model used to generate predictions. By default,
  it is set to `NULL`, indicating no dropout.

- m_dropout:

  The number of dropout time intervals to extrapolate the hazard
  function beyond the last observed dropout time.

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

  The number of replications for simulation. By default, it is set
  to 500. If `newSubjects` is not `NULL`, the number of draws in
  `newSubjects` should be `nreps`.

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

  A Boolean variable to control whether or not to show the prediction
  plot. By default, it is set to `TRUE`.

- by_treatment:

  A Boolean variable to control whether or not to predict event by
  treatment group. By default, it is set to `FALSE`.

- covariates_event:

  The names of baseline covariates from the input data frame to include
  in the event model, e.g., c("age", "sex"). Factor variables need to be
  declared in the input data frame.

- event_fit_with_covariates:

  The pre-fitted event model with covariates used to generate event
  predictions for ongoing subjects.

- covariates_dropout:

  The names of baseline covariates from the input data frame to include
  in the dropout model, e.g., c("age", "sex"). Factor variables need to
  be declared in the input data frame.

- dropout_fit_with_covariates:

  The pre-fitted dropout model with covariates used to generate dropout
  predictions for ongoing subjects.

- fix_parameter:

  Whether to fix parameters at the maximum likelihood estimates when
  generating new data for prediction. Defaults to FALSE, in which case,
  parameters will be drawn from their approximate posterior
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

A list of prediction results which includes important information such
as the median, lower and upper percentiles for the estimated day and
date to reach the target number of events, as well as simulated event
data for both ongoing and new subjects. The data for the prediction plot
is also included within this list. If target_t is specified, it
additionally provides the median, lower, and upper percentiles of the
event count at target_t, as well as the predictive probability of
achieving the target number of events by target_t.

## Details

To ensure successful event prediction at the design stage, it is
important to provide the `newSubjects` data set.

To specify the event (dropout) model used during the design-stage event
prediction, the `event_fit` (`dropout_fit`) should be a list with one
element per treatment. For each treatment, the element should include
`model` to specify the event model (exponential, weibull, log-logistic,
log-normal, or piecewise exponential), and `theta` and `vtheta` to
indicate the parameter values and the covariance matrix. For the
piecewise exponential event (dropout) model, the list should also
include `piecewiseSurvivalTime` (`piecewiseDropoutTime`) to indicate the
location of knots. It should be noted that the model averaging and
spline options are not appropriate for use during the design stage.

Following the commencement of the trial, we obtain the event model fit
and the dropout model fit based on the observed data, denoted as
`event_fit` and `dropout_fit`, respectively. These fitted models are
subsequently utilized to generate event and dropout times for both
ongoing and new subjects in the trial.

## References

Emilia Bagiella and Daniel F. Heitjan. Predicting analysis times in
randomized clinical trials. Stat in Med. 2001; 20:2055-2063.

Gui-shuang Ying and Daniel F. Heitjan. Weibull prediction of event times
in clinical trials. Pharm Stat. 2008; 7:107-120.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

# Event prediction after enrollment completion
set.seed(2000)

event_fits <- fitEvent(
  df = interimData2,
  event_model = "piecewise exponential",
  piecewiseSurvivalTime = c(0, 140, 352),
  nthreads = 1)

dropout_fits <- fitDropout(
  df = interimData2,
  dropout_model = "exponential",
  nthreads = 1)

event_pred <- predictEvent(
  df = interimData2,
  target_d = 200,
  event_fit = event_fits$fit,
  dropout_fit = dropout_fits$fit,
  pilevel = 0.90,
  nreps = 100,
  nthreads = 1)
#> Time from cutoff until 200 events: 129 days 
#>  Median prediction date: 2021-02-26 
#>  Prediction interval: 2021-01-09, 2021-05-09 
```
