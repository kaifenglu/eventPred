# Summarize observed data

Provides an overview of the observed data, including the trial start
date, data cutoff date, enrollment duration, number of subjects
enrolled, number of events and dropouts, number of subjects at risk,
cumulative enrollment and event data, daily enrollment rates, and
Kaplan-Meier plots for time to event and time to dropout.

## Usage

``` r
summarizeObserved(
  df,
  to_predict = "event only",
  showplot = TRUE,
  by_treatment = FALSE,
  generate_plot = TRUE,
  interactive_plot = TRUE,
  nthreads = 0
)
```

## Arguments

- df:

  The subject-level data, including `trialsdt`, `usubjid`, `randdt`, and
  `cutoffdt` for enrollment prediction, as well as `time`, `event` and
  `dropout` for event prediction, and `treatment` coded as 1, 2, and so
  on, and `treatment_description` for prediction by treatment group.

- to_predict:

  Specifies what to predict: "enrollment only", "event only", or
  "enrollment and event". By default, it is set to "event only".

- showplot:

  A Boolean variable to control whether or not to show the observed data
  plots. By default, it is set to `TRUE`.

- by_treatment:

  A Boolean variable to control whether or not to summarize observed
  data by treatment group. By default, it is set to `FALSE`.

- generate_plot:

  Whether to generate plots.

- interactive_plot:

  Whether to produce interactive plots using plotly or static plots
  using ggplot2.

- nthreads:

  Integer number of threads to use for \`data.table' (0 means the
  default data.table behavior).

## Value

A list that includes a range of summary statistics, data sets, and plots
depending on the value of `to_predict`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

observed1 <- summarizeObserved(
  df = interimData1,
  to_predict = "enrollment and event",
  nthreads = 1)

observed2 <- summarizeObserved(
  df = interimData2,
  to_predict = "event only",
  nthreads = 1)
```
