# eventPred 0.2.0

- allow the use of treatment labels for by-treatment prediction
- include usubjid in subject-level data sets
- use the quantile method for predicted date if all simulated data sets attain the target number of events
- add log-logistic event model and log-logistic dropout model
- change parameterization of Weibull distribution to be consistent with log-logistic and log-normal distributions in the AFT family 
- add AIC to enrollment, event and dropout model fits
- check the required number of events/dropouts for event/dropout model fits 
- add "spline" as a dropout_model option
- update enroll_fit, event_fit, and dropout_fit for prior incorporation


# eventPred 0.1.5

- update design stage prediction with one treatment arm
- allow ongoing subjects with last known date before data cutoff
- update the calculation of ongoing subjects to accommodate ongoing subjects with last known date before data cutoff
- update time for new subjects to start with day 1 and update totalTime calculation for newEvents to remove double count of day 1
- update predictEnrollment to remove calculation of d0, c0, and r0
- add names to event_pred_day
- add nyears and nreps to prediction results

# eventPred 0.1.4

- add validity checks for input dataset variables
- update totalTime calculation for observed data
- use method="Nelder-Mead" as the default optimization algorithm for flexsurvspline
- add by-treatment prediction

# eventPred 0.1.3

- update the description of internal datasets
- update summarizeObserved to remove adt from adsl
- add Royston and Parmar (2002) spline event model

# eventPred 0.1.2

- add mean and variance to prediction output
- update the BIC weight for model averaging
- add more details for model fit parameters 
- add day 1 to enrollment plot
- allow prior piecewise Poisson enrollment and piecewise exponential event or dropout models to have additional cut points beyond the observed data range
- update internal data sets

# eventPred 0.1.1

- add stage and to_predict information in getPrediction output
- add the cutoff time point to the number of ongoing subjects
- change the default model for dropout to exponential
- require trialsdt in input data set

# eventPred 0.1.0

- added the piecewise Poisson model to fitEnrollment and predictEnrollment at the analysis stage
- added number of dropouts
- added number of subjects at risk
- added a data set when the enrollment has completed
- corrected the x-axis title for predictEnrollment and predictEvent
- updated alogrithm to allow one piece piecewise Poisson enrollment model and one piece piecewise exponential time-to-event model
- modified the weight calculation for model averaging to avoid underflow
- used weighted BIC for model averaging
- removed the dropout_model parameter for summarizeObserved
- changed the default number of knots of the b-spline enrollment model to zero
- replaced first and last with slice of dplyr in summarizeObserved
- improved the initial value for the time-decay enrollment model parameters
- added showplot to fitEnrollment, fitEvent and fitDropout
- sped up the calculations of quantiles
- added target_n to predictEnrollment output and target_d to predictEvent output
- removed the cutoff date from ongoing_pred_df before data cutoff
- restricted enrollment model fitting to the last randomization date
- added piecewise exponential dropout model
- use delta method to obtain the variance of model parameters for pooled population
- replace randomization probabilities with treatment allocation within a randomization block
- allow number of subjects to differ among simulated data sets
- remove custom date axis

# eventPred 0.0.1

- Initial release

