# Event Prediction Incorporating Prior Information

``` r

library(eventPred)
```

We analyzed interim enrollment and event data from the interimData1
dataset within the eventPred package. As of the cutoff date of March 1,
2019, 224 patients had been enrolled (with a target of 300), and 42
patients had experienced the event of interest (with a target of 200
events).

A time-decay model was used to predict future enrollment patterns. A
Weibull distribution was employed to model the time to event. To
accommodate the limited number of dropouts (n = 1) while allowing for
potential non-constant dropout rates, we modeled time to dropout using a
flexible piecewise exponential prior distribution with a changepoint at
day 140.

``` r

set.seed(2000)

dropout_prior <- list(
  model = "piecewise exponential",
  theta = c(-8, -9),
  vtheta = diag(2)*1e-8,
  piecewiseDropoutTime = c(0, 140))

pred <- getPrediction(
  df = interimData1,
  to_predict = "enrollment and event",
  target_n = 300,
  target_d = 200,
  enroll_model = "time-decay",
  event_model = "weibull",
  dropout_model = "piecewise exponential",
  piecewiseDropoutTime = 0,
  dropout_prior = dropout_prior,
  pilevel = 0.90, 
  nreps = 500)
#> Time from cutoff until 300 subjects: 72 days 
#>  Median prediction date: 2019-05-11 
#>  Prediction interval: 2019-04-26, 2019-06-24 
#> Time from cutoff until 200 events: 933 days 
#>  Median prediction date: 2021-09-18 
#>  Prediction interval: 2020-07-29, 2023-10-24
```

The median predicted date to reach 200 events is September 18, 2021. The
90% prediction interval for event completion is wide, ranging from July
29, 2020 to October 24, 2023, indicating uncertainty in the prediction.
Compared with the prediction without prior information, the point
estimate moved upwards and the prediction interval widened due to a
higher dropout rate in the prior distribution.
