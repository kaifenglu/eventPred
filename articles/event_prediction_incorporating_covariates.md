# Event Prediction Incorporating Covariates

``` r

library(eventPred)
```

We analyzed interim event data from the interimData2 dataset within the
eventPred package. As of the cutoff date of October 1, 2020, all 300
patients had been enrolled, and 183 patients had experienced the event
of interest (with a target of 200 events).

A Weibull distribution was employed to model the time to event
incorporating a factor covariate with three levels. An exponential
distribution was used to model the time to dropout.

``` r

set.seed(2000)

df = interimData2
df$agegrp = as.factor(sample(c(1,2,3), nrow(df), replace = TRUE, 
                             prob = c(0.3, 0.6, 0.1)))

pred <- getPrediction(
  df = df,
  to_predict = "event only",
  target_d = 200,
  event_model = "weibull",
  covariates_event = "agegrp",
  dropout_model = "exponential",
  pilevel = 0.90, 
  nyears = 1,
  nreps = 500)
#> Time from cutoff until 200 events: 132 days 
#>  Median prediction date: 2021-03-01 
#>  Prediction interval: 2021-01-13, 2021-05-13
```

The median predicted date to reach 200 events is March 1, 2021. The 90%
prediction interval for event completion ranges from January 13, 2021 to
May 13, 2021. The negligible impact of baseline covariate inclusion on
the results aligns with its random generation and lack of association
with the time to event distribution.
