# eventPred 0.2.5

* rename the components of fitEnrollment output to fit and fit_plot
* restructure the outputs of fitEvent and fitDropout into a list for by-treatment analysis of model fitting and visualization, where each element in the list corresponds to a specific treatment group and has a dedicated sub-list containing two components with one for fit and the other for fit_plot
* update predictEvent.R, getPrediction.R and app.R accordingly to accommodate the new structure of fitEnrollment, fitEvent and fitDropout outputs
* update the output of event_prediction_after_enrollment_completion vignette
* add the condition of (!is.null(event_fit())) for event_fit_ic and (!is.null(dropout_fit())) for dropout_fit_ic in the shiny app
* minor change to the ui layout of the shiny app
* ensure that the randomization date for new patients is after the cutoff date and the event date for ongoing subjects is after the cutoff date


# eventPred 0.2.4

* fitEnrollment.R

    - replace round with formatC to retain the zeros after the decimal point

* fitEvent.R

    - parameterize the exponential distribution in terms of log(rate)
    
    - update the requirement for fitting a piecewise exponential model
    
    - update the call to the pwexpreg function

    - ensure the sub plots align on the x axis
    
    - export the sub plots as a list instead of a plotly subplot object
    
    - replace round with formatC to retain the zeros after the decimal point
    
* fitDropout.R

    - parameterize the exponential distribution in terms of log(rate)

    - update the requirement for fitting a piecewise exponential model
    
    - update the call to the pwexpreg function
    
    - ensure the sub plots align on the x axis
    
    - export the sub plots as a list instead of a plotly subplot object

    - replace round with formatC to retain the zeros after the decimal point
    
* predictEnrollment.R

    - add the 'name' parameter to the Plotly traces to ensure proper legends 
    
    - export the sub plots as a list instead of a plotly subplot object
    
* predictEvent.R

    - parameterize the exponential distribution in terms of log(rate)

    - ensure simulated time >= 1
    
    - add the 'name' parameter to the Plotly traces to ensure proper legends 
    
    - export the sub plots as a list instead of a plotly subplot object

* getPrediction.R

    - check the input data to ensure all required columns are present
    
    - check the input data to ensure none of the required columns have missing values
    
    - add treatment_description to the input data when treatment is present but treatment_description is missing
    
    - parameterize the exponential distribution in terms of log(rate)
    
    - obtain event_fit (event_fit_with_covariates) without regard of the existence of event_prior (event_prior_with_covariates)

    - obtain dropout_fit (dropout_fit_with covariates) without regard of the existence of dropout_prior (dropout_prior_with_covariates)
    
    
* utilities.R
    
    - update the pwexpreg function so that its parameters are consistent with other piecewise exponential functions 
    
    - use the Brent method to fit the piecewise exponential regression model with only one interval and no covariates

* launchShinyApp.R

    - newly added to launch the Shiny app for event prediction

* vignettes

    - add event_prediction_at_the_design_stage.Rmd
    
    - add event_prediction_before_enrollment_completion.Rmd
    
    - add event_prediction_after_enrollment_completion.Rmd
    
    - add event_prediction_incorporating_prior_information.Rmd
    
    - add event_prediction_incorporating_covariates.Rmd
    

# eventPred 0.2.3

* predictEvent.R
    
    - check to make sure that dropout_fit is not null before simulating dropout times for new and ongoing patients

# eventPred 0.2.2

* eventPred-package.R

    - remove import tmvtnsim rtnorm
    
    - add import purrr list_c map map_dbl
    
    - add import stats as.formula model.matrix qlnorm rlogis
  
* utilities.R
    - add pmodavg for the distribution of model averaging of Weibull and log-normal
    
    - add ppwexp and qpwexp functions for the piecewise exponential distribution
    
    - add llik_pwexp for the log-likelihood of piecewise exponential regression
    
    - add pwexpreg for the regression analysis of piecewise exponential distribution
  
* fitEnrollment.R
    - use the hessian option in optim to remove the optimHess call in fitEnrollment
  
* fitEvent.R
    - add covariates to the fitEvent function to fit regression models

* fitDropout.R
    - add covariates to the fitDropout function to fit regression models
  
* predictEvent.R
    - add covariates_event, event_fit_with_covariates, covariates_dropout, dropout_fit_with_covariates
    
    - fit the event model with covariates if event_fit_with_covariates is not NULL, and fit the event model without covariates otherwise
    
    - fit the dropout model with covariates if dropout_fit_with_covariates is not NULL, and fit the dropout model without covariates otherwise
    
    - generate the event time for new patients separately from the event time for ongoing patients
    
    - generate the dropout time for new patients separately from the dropout time for ongoing patients
    
    - apply ceiling to the derived time after comparison of generated survivalTime and dropoutTime

* getPrediction.R  

    - add covariates_event, event_prior_with_covariates, covariates_dropout, dropout_prior_with_covariates
    
    - add penalized log-likelihood (posterior) function with covariates for exponential, Weibull, log-logistic, log-normal, and piecewise exponential distributions
    
    - simplify the algorithm for combining prior distributions across treatments
    
    - fit event/dropout models with or without covarites depending on the study stage and the presence/absence of covariates_event and covariates_dropout
    
    - add subject_data to the output
  

# eventPred 0.2.1

* remove the factor attribute of the treatment_description variable

* add pilevel in the output data set for prediction interval level

* replace treatment_label with treatment_description in observed data for enrollment prediction

* update the upper bound of the cutoff reference line in prediction plot

* retain the plots of enroll_fit, event_fit, and dropout_fit in getPrediction output

* add usubjid and treatment_description to the internal data sets 

* round the simulated arrivalTime and time so that the time can be interpreted in days

# eventPred 0.2.0

* allow the use of treatment labels for by-treatment prediction

* include usubjid in subject-level data sets

* use the quantile method for predicted date if all simulated data sets attain the target number of events

* add log-logistic event model and log-logistic dropout model

* change parameterization of Weibull distribution to be consistent with log-logistic and log-normal 
distributions in the AFT family 

* add AIC to enrollment, event and dropout model fits

* check the required number of events/dropouts for event/dropout model fits 

* add "model averaging" and "spline" as additional dropout_model options

* update enroll_fit, event_fit, and dropout_fit for prior incorporation


# eventPred 0.1.5

* update design stage prediction with one treatment arm

* allow ongoing subjects with last known date before data cutoff

* update the calculation of ongoing subjects to accommodate ongoing subjects with last known date before data cutoff

* update time for new subjects to start with day 1 and update totalTime calculation for newEvents to remove double count of day 1

* update predictEnrollment to remove calculation of d0, c0, and r0

* add names to event_pred_day

* add nyears and nreps to prediction results

# eventPred 0.1.4

* add validity checks for input dataset variables

* update totalTime calculation for observed data

* use method="Nelder-Mead" as the default optimization algorithm for flexsurvspline

* add by-treatment prediction

# eventPred 0.1.3

* update the description of internal datasets

* update summarizeObserved to remove adt from adsl

* add Royston and Parmar (2002) spline event model

# eventPred 0.1.2

* add mean and variance to prediction output

* update the BIC weight for model averaging

* add more details for model fit parameters 

* add day 1 to enrollment plot

* allow prior piecewise Poisson enrollment and piecewise exponential event or dropout models to have additional cut points beyond the observed data range

* update internal data sets

# eventPred 0.1.1

* add stage and to_predict information in getPrediction output

* add the cutoff time point to the number of ongoing subjects

* change the default model for dropout to exponential

* require trialsdt in input data set

# eventPred 0.1.0

* added the piecewise Poisson model to fitEnrollment and predictEnrollment at the analysis stage

* added number of dropouts

* added number of subjects at risk

* added a data set when the enrollment has completed

* corrected the x-axis title for predictEnrollment and predictEvent

* updated alogrithm to allow one piece piecewise Poisson enrollment model and one piece piecewise exponential time-to-event model

* modified the weight calculation for model averaging to avoid underflow

* used weighted BIC for model averaging

* removed the dropout_model parameter for summarizeObserved

* changed the default number of knots of the b-spline enrollment model to zero

* replaced first and last with slice of dplyr in summarizeObserved

* improved the initial value for the time-decay enrollment model parameters

* added showplot to fitEnrollment, fitEvent and fitDropout

* sped up the calculations of quantiles

* added target_n to predictEnrollment output and target_d to predictEvent output

* removed the cutoff date from ongoing_pred_df before data cutoff

* restricted enrollment model fitting to the last randomization date

* added piecewise exponential dropout model

* use delta method to obtain the variance of model parameters for pooled population

* replace randomization probabilities with treatment allocation within a randomization block

* allow number of subjects to differ among simulated data sets

* remove custom date axis

# eventPred 0.0.1

* Initial release

