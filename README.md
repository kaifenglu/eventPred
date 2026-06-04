# eventPred

`eventPred` predicts enrollment and event timing in clinical trials.
It supports both:

- Design stage prediction using model assumptions and optional priors.
- Analysis stage prediction using observed blinded or unblinded data.

The package provides enrollment modeling, time-to-event modeling,
time-to-dropout modeling, simulation-based prediction intervals,
and an interactive Shiny app.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("kaifenglu/eventPred")
```

## Key Features

- Enrollment models: Poisson, time-decay, B-spline, piecewise Poisson.
- Event/dropout models: exponential, Weibull, log-logistic, log-normal,
  piecewise exponential, model averaging, spline, and Cox.
- By-treatment prediction and optional baseline covariates.
- Prediction intervals based on simulation (`nreps`).
- Support for fixed follow-up and variable follow-up designs.
- Built-in example datasets: `interimData1`, `interimData2`, `finalData`.

## Typical Workflow

1. Summarize observed data with `summarizeObserved()`.
2. Fit enrollment, event, and dropout models with `fitEnrollment()`,
   `fitEvent()`, and `fitDropout()`.
3. Generate predictions with:
   - `predictEnrollment()` for enrollment only
   - `predictEvent()` for event timing only
   - `getPrediction()` for end-to-end enrollment and event prediction

## Minimal Example

```r
library(eventPred)

# Event prediction after enrollment completion
set.seed(3000)

pred <- getPrediction(
  df = interimData2,
  to_predict = "event only",
  target_d = 200,
  event_model = "weibull",
  dropout_model = "exponential",
  pilevel = 0.90,
  nreps = 100
)

pred$event_pred$event_pred_summary
```

## Model-Fit + Prediction Example

```r
library(eventPred)

set.seed(2000)

event_fits <- fitEvent(
  df = interimData2,
  event_model = "piecewise exponential",
  piecewiseSurvivalTime = c(0, 140, 352)
)

dropout_fits <- fitDropout(
  df = interimData2,
  dropout_model = "exponential"
)

event_pred <- predictEvent(
  df = interimData2,
  target_d = 200,
  event_fit = event_fits$fit,
  dropout_fit = dropout_fits$fit,
  pilevel = 0.90,
  nreps = 100
)

event_pred$event_pred_summary
```

## Design Stage Example

```r
library(eventPred)

set.seed(1000)

enroll_pred <- predictEnrollment(
  target_n = 300,
  enroll_fit = list(
    model = "piecewise poisson",
    theta = log(26 / 9 * seq(1, 9) / 30.4375),
    vtheta = diag(9) * 1e-8,
    accrualTime = seq(0, 8) * 30.4375
  ),
  pilevel = 0.90,
  nreps = 100
)

enroll_pred$enroll_pred_summary
```

## Run the Shiny App

```r
library(eventPred)
runShinyApp_eventPred()
```

## Time Unit

The package uses days as the primary time unit.
To convert rates per month to rates per day, divide by 30.4375.

## Documentation

- Website: <https://kaifenglu.github.io/eventPred/>
- Function reference: <https://kaifenglu.github.io/eventPred/reference/>
- Issues: <https://github.com/kaifenglu/eventPred/issues>

## Citation

If you use `eventPred` in analysis or reporting, please cite relevant
methodology references included in the package documentation.
