# Event Prediction before Enrollment Completion

``` r

library(eventPred)
```

We analyzed interim enrollment and event data from the interimData1
dataset within the eventPred package. As of the cutoff date of March 1,
2019, 224 patients had been enrolled (with a target of 300), and 42
patients had experienced the event of interest (with a target of 200
events).

A time-decay model was used to predict future enrollment patterns. A
Weibull distribution was employed to model the time to event. Dropout
modeling was not performed due to the limited number of dropouts (only
one).

``` r

set.seed(2000)

pred <- getPrediction(
  df = interimData1,
  to_predict = "enrollment and event",
  target_n = 300,
  target_d = 200,
  enroll_model = "time-decay",
  event_model = "weibull",
  dropout_model = "none",
  pilevel = 0.90, 
  nreps = 500)
#> Time from cutoff until 300 subjects: 72 days 
#>  Median prediction date: 2019-05-11 
#>  Prediction interval: 2019-04-26, 2019-06-24 
#> Time from cutoff until 200 events: 835 days 
#>  Median prediction date: 2021-06-12 
#>  Prediction interval: 2020-06-17, 2022-12-18
```

The median predicted date to reach 200 events is June 12, 2021. The 90%
prediction interval for event completion is wide, ranging from June 17,
2020 to December 18, 2022, indicating uncertainty in the prediction. The
wide prediction interval is primarily attributed to the limited number
of events observed thus far (42 out of 200 target events). More events
are needed to refine the prediction and narrow the prediction interval.
