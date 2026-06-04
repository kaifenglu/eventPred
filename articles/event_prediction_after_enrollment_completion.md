# Event Prediction after Enrollment Completion

``` r

library(eventPred)
```

We analyzed interim event data from the interimData2 dataset within the
eventPred package. As of the cutoff date of October 1, 2020, all 300
patients had been enrolled, and 183 patients had experienced the event
of interest (with a target of 200 events).

A Weibull distribution was employed to model the time to event. An
exponential distribution was used to model the time to dropout.

``` r

set.seed(2000)

pred <- getPrediction(
  df = interimData2,
  to_predict = "event only",
  target_d = 200,
  event_model = "weibull",
  dropout_model = "exponential",
  pilevel = 0.90, 
  nyears = 1,
  nreps = 500)
#> Time from cutoff until 200 events: 139 days 
#>  Median prediction date: 2021-03-08 
#>  Prediction interval: 2021-01-12, 2021-05-20
```

The median predicted date to reach 200 events is March 8, 2021. The 90%
prediction interval for event completion is much narrower, ranging from
January 12, 2021 to May 20, 2021, indicating reduced uncertainty in the
prediction attributed to the greatly increased number of events.

Based on the finalData dataset within the eventPred package, the actual
date to reach the 200 events was February 14, 2021. The increased number
of events greatly improved the prediction accuracy and precision.
