# eventPred 0.0.2

- added the piecewise Poisson model to fitEnrollment and predictEnrollment at the analysis stage
- added number of dropouts
- added number of subjects at risk
- added a data set when the enrollment has completed
- modified the weight calculation for model averaging to avoid underflow
- used weighted BIC for model averaging
- added showplot parameter to summarizeObserved, fitEnrollment, fitEvent, and  fitDropout
- change the default number of knots of the b-spline enrollment model to zero
- replaced first and last with slice of dplyr in summarizeObserved
- improved the initial value for the time-decay enrollment model parameters
- sped up the calculations of quantiles

# eventPred 0.0.1

- Initial release

