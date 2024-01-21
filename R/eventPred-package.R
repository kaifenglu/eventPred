#' @name eventPred-package
#' @aliases eventPred-package
#' @keywords internal
#' "_PACKAGE"
#'
#' @title Event Prediction
#'
#' @description Predicts enrollment and events at the design stage
#' using assumed enrollment and treatment-specific time-to-event models,
#' or at the analysis stage using blinded or unblinded data and
#' specified enrollment and time-to-event models through simulations.
#'
#' @details Accurately predicting the date at which a target number
#' of subjects or events will be achieved is critical for the planning,
#' monitoring, and execution of clinical trials. The \code{eventPred}
#' package provides enrollment and event prediction capabilities
#' using assumed enrollment and treatment-specific time-to-event models
#' at the design stage, using blinded or unblinded data and
#' specified enrollment and time-to-event models at the analysis stage.
#'
#' At the design stage, enrollment is often specified using a
#' piecewise Poisson process with a constant enrollment rate
#' during each specified time interval. At the analysis stage,
#' before enrollment completion, the \code{eventPred} package
#' considers several models, including the homogeneous Poisson
#' model, the time-decay model with an enrollment
#' rate function \code{lambda(t) = mu/delta*(1 - exp(-delta*t))},
#' the B-spline model with the daily enrollment rate
#' \code{lambda(t) = exp(B(t)*theta)}, and the piecewise Poisson model.
#' If prior information exists on the model parameters, it can
#' be combined with the likelihood to yield the posterior distribution.
#'
#' The \code{eventPred} package also offers several time-to-event models,
#' including exponential, Weibull, log-logistic, log-normal, piecewise
#' exponential, model averaging of Weibull and log-normal, and spline.
#' For time to dropout, the same set of model options are considered.
#' If enrollment is complete, ongoing subjects who have not had the event
#' of interest or dropped out of the study before the data cut contribute
#' additional events in the future. Their event times are generated
#' from the conditional distribution given that they have survived
#' at the data cut. For new subjects that need to be enrolled,
#' their enrollment time and event time can be generated from the
#' specified enrollment and time-to-event models with parameters
#' drawn from the posterior distribution. Time-to-dropout can be
#' generated in a similar fashion.
#'
#' The \code{eventPred} package displays the Akaike Information
#' Criterion (AIC), the Bayesian Information
#' Criterion (BIC) and a fitted curve overlaid with observed data
#' to help users select the most appropriate model for enrollment
#' and event prediction. Prediction intervals in the prediction plot
#' can be used to measure prediction uncertainty, and the simulated
#' enrollment and event data can be used for further data exploration.
#'
#' The most useful function in the \code{eventPred} package is
#' \code{getPrediction}, which combines model fitting, data simulation,
#' and a summary of simulation results. Other functions perform
#' individual tasks and can be used to select an appropriate
#' prediction model.
#'
#' The \code{eventPred} package implements a model
#' parameterization that enhances the asymptotic normality of
#' parameter estimates. Specifically, the package utilizes the
#' following parameterization to achieve this goal:
#' \itemize{
#'   \item Enrollment models
#'   \itemize{
#'     \item Poisson: \code{theta = log(rate)}
#'     \item Time-decay: \code{theta = c(log(mu), log(delta))}
#'     \item B-spline: no reparametrization is needed. The
#'     knots as considered fixed.
#'     \item Piecewise Poisson: \code{theta = log(rates)}.
#'     The left endpoints of time intervals, denoted as
#'     \code{accrualTime}, are considered fixed.
#'   }
#'
#'   \item Event or dropout models
#'   \itemize{
#'     \item Exponential: \code{theta = log(rate)}
#'     \item Weibull: \code{theta = c(log(scale), -log(shape))}
#'     \item Log-logistic: \code{theta = c(log(scale), -log(shape))}
#'     \item Log-normal: \code{theta = c(meanlog, log(sdlog))}
#'     \item Piecewise exponential: \code{theta = log(rates)}.
#'     The left endpoints of time intervals, denoted as
#'     \code{piecewiseSurvivalTime} for event model and
#'     \code{piecewiseDropoutTime} for dropout model, are
#'     considered fixed.
#'     \item Model averaging: \code{theta = c(log(weibull$scale),
#'     -log(weibull$shape), lnorm$meanlog, log(lnorm$sdlog))}.
#'     The covariance matrix for \code{theta} is structured
#'     as a block diagonal matrix, with the upper-left block
#'     corresponding to the Weibull component and the
#'     lower-right block corresponding to the log-normal
#'     component. In other words, the covariance matrix is
#'     partitioned into two distinct blocks, with no
#'     off-diagonal elements connecting the two components.
#'     The weight assigned to the Weibull component, denoted as
#'     \code{w1}, is considered fixed.
#'     \item Spline: \code{theta} corresponds to the coefficients
#'     of basis vectors. The \code{knots} and \code{scale}
#'     are considered fixed. The \code{scale} can be hazard,
#'     odds, or normal, corresponding to extensions of Weibull,
#'     log-logistic, and log-normal distributions, respectively.
#'   }
#' }
#'
#' The \code{eventPred} package uses days as its primary time unit.
#' If you need to convert enrollment or event rates per month to
#' rates per day, simply divide by 30.4375.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Emilia Bagiella and Daniel F. Heitjan. Predicting analysis times in
#' randomized clinical trials. Stat in Med. 2001; 20:2055-2063.
#'
#' Gui-shuang Ying and Daniel F. Heitjan. Weibull prediction of event
#' times in clinical trials. Pharm Stat. 2008; 7:107-120.
#'
#' Xiaoxi Zhang and Qi Long. Stochastic modeling and prediction for
#' accrual in clinical trials. Stat in Med. 2010; 29:649-658.
#'
#' Patrick Royston and Mahesh K. B. Parmar. Flexible parametric
#' proportional-hazards and proportional-odds models for censored
#' survival data, with application to prognostic modelling and
#' estimation of treatment effects. Stat in Med. 2002; 21:2175-2197.
#'
#'
#' @importFrom dplyr %>% arrange as_tibble bind_rows cross_join filter
#'   group_by mutate n rename rename_all row_number select slice
#'   summarize tibble
#' @importFrom plotly plot_ly add_lines add_ribbons hide_legend layout
#' @importFrom survival Surv survfit survreg
#' @importFrom splines bs
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm rmvnorm
#' @importFrom rstpm2 vuniroot
#' @importFrom numDeriv grad
#' @importFrom purrr list_c map map_dbl
#' @importFrom stats as.formula dweibull dlnorm loess.smooth model.matrix
#'   optim optimHess pexp plnorm plogis pnorm pweibull qlogis qlnorm qnorm
#'   quantile rbinom reorder rexp rlnorm rlogis rmultinom rnorm runif
#'   rweibull uniroot var
#' @importFrom flexsurv flexsurvspline psurvspline qsurvspline rsurvspline
#'   dllogis pllogis
#' @importFrom erify check_bool check_class check_content check_n
#'   check_positive
#' @importFrom rlang .data
#' @importFrom shiny shinyAppDir
#'
NULL
