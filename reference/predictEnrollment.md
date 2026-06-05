# Predict enrollment

Utilizes a pre-fitted enrollment model to generate enrollment times for
new subjects and provide a prediction interval for the expected time to
reach the enrollment target.

## Usage

``` r
predictEnrollment(
  df = NULL,
  target_n = NA,
  enroll_fit = NULL,
  lags = 30,
  pilevel = 0.9,
  nyears = 4,
  nreps = 500,
  showsummary = TRUE,
  showplot = TRUE,
  by_treatment = FALSE,
  ngroups = 1,
  alloc = NULL,
  treatment_label = NULL,
  fix_parameter = FALSE,
  generate_plot = TRUE,
  interactive_plot = TRUE,
  nthreads = 0
)
```

## Arguments

- df:

  The subject-level enrollment data, including `trialsdt`, `randdt` and
  `cutoffdt`. The data should also include `treatment` coded as 1, 2,
  and so on, and `treatment_description` for prediction by treatment
  group. By default, it is set to `NULL` for enrollment prediction at
  the design stage.

- target_n:

  The target number of subjects to enroll in the study.

- enroll_fit:

  The pre-fitted enrollment model used to generate predictions.

- lags:

  The day lags to compute the average enrollment rate to carry forward
  for the B-spline enrollment model. By default, it is set to 30.

- pilevel:

  The prediction interval level. By default, it is set to 0.90.

- nyears:

  The number of years after the data cut for prediction. By default, it
  is set to 4.

- nreps:

  The number of replications for simulation. By default, it is set to
  500.

- showsummary:

  A Boolean variable to control whether or not to show the prediction
  summary. By default, it is set to `TRUE`.

- showplot:

  A Boolean variable to control whether or not to show the prediction
  plot. By default, it is set to `TRUE`.

- by_treatment:

  A Boolean variable to control whether or not to predict enrollment by
  treatment group. By default, it is set to `FALSE`.

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

- fix_parameter:

  Whether to fix parameters at the maximum likelihood estimates when
  generating new data for prediction. Defaults to FALSE, in which case,
  parameters will be drawn from their approximate posterior
  distributions.

- generate_plot:

  Whether to generate plots.

- interactive_plot:

  Whether to produce interactive plots using plotly or static plots
  using ggplot2.

- nthreads:

  Integer number of threads to use for \`data.table' (0 means the
  default data.table behavior).

## Value

A list of prediction results, which includes important information such
as the median, lower and upper percentiles for the estimated time to
reach the target number of subjects, as well as simulated enrollment
data for new subjects. The data for the prediction plot is also included
within the list.

## Details

The `enroll_fit` variable can be used for enrollment prediction at the
design stage. A piecewise Poisson model can be parameterized through the
time intervals, `accrualTime`, which is treated as fixed, and the
enrollment rates in the intervals, `accrualIntensity`, the log of which
is used as the model parameter. For the homogeneous Poisson, time-decay,
and piecewise Poisson models, `enroll_fit` is used to specify the prior
distribution of model parameters, with a very small variance being used
to fix the parameter values. It should be noted that the B-spline model
is not appropriate for use during the design stage.

During the enrollment stage, `enroll_fit` is the enrollment model fit
based on the observed data. The fitted enrollment model is used to
generate enrollment times for new subjects.

## References

Xiaoxi Zhang and Qi Long. Stochastic modeling and prediction for accrual
in clinical trials. Stat in Med. 2010; 29:649-658.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Enrollment prediction at the design stage
set.seed(1000)

enroll_pred <- predictEnrollment(
  target_n = 300,
  enroll_fit = list(
    model = "piecewise poisson",
    theta = log(26/9*seq(1, 9)/30.4375),
    vtheta = diag(9)*1e-8,
    accrualTime = seq(0, 8)*30.4375),
  pilevel = 0.90,
  nreps = 100,
  nthreads = 1)
#> Time from trial start until 300 subjects 
#>  Median prediction day: 469 
#>  Prediction interval: 443, 509 
```
