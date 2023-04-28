#' @title Enrollment and event prediction
#' @description Performs enrollment and event prediction by utilizing
#'   observed data and specified enrollment and event models.
#'
#' @param df The subject-level enrollment and event data, including
#'   \code{trialsdt}, \code{randdt}, \code{cutoffdt}, \code{time},
#'   \code{event}, and \code{dropout}. By default, it is set to
#'   \code{NULL} for enrollment and event prediction at the design stage.
#' @param to_predict Specifies what to predict: "enrollment only", "event
#'   only", or "enrollment and event". By default, it is set to
#'   "enrollment and event".
#' @param target_n The target number of subjects to enroll in the study.
#' @param target_d The target number of events to reach in the study.
#' @param enroll_model The enrollment model which can be specified as
#'   "Poisson", "Time-decay", "B-spline", or
#'   "Piecewise Poisson". By default, it is set to "B-spline".
#' @param nknots The number of inner knots for the B-spline enrollment
#'   model. By default, it is set to 0.
#' @param lags The day lags to compute the average enrollment rate to
#'   carry forward for the B-spline enrollment model. By default,
#'   it is set to 30.
#' @param accrualTime The accrual time intervals for the piecewise
#'   Poisson model. Must start with 0, e.g., c(0, 30) breaks the
#'   time axis into 2 accrual intervals: [0, 30) and [30, Inf).
#'   By default, it is set to 0.
#' @param enroll_model_parameter The enrollment model parameters for
#'   design-stage enrollment prediction.
#' @param event_model The event model used to analyze the event data
#'   which can be set to one of the
#'   following options: "exponential", "Weibull", "log-normal",
#'   "piecewise exponential", or "model averaging". The model averaging
#'   uses the \code{exp(-bic/2)} weighting and combines Weibull and
#'   log-normal models. By default, it is set to "model
#'   averaging".
#' @param piecewiseSurvivalTime A vector that specifies the time
#'   intervals for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param event_model_parameter The event model parameters for
#'   design-stage event prediction.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options: "exponential",
#'   "Weibull", "log-normal", or "piecewise exponential". By default,
#'   it is set to "exponential".
#' @param piecewiseDropoutTime A vector that specifies the time
#'   intervals for the piecewise exponential dropout distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param dropout_model_parameter The dropout model parameters for
#'   design-stage event prediction.
#' @param fixedFollowup A Boolean variable indicating whether a fixed
#'   follow-up design is used. By default, it is set to \code{FALSE}
#'   for a variable follow-up design.
#' @param followupTime The follow-up time for a fixed
#'   follow-up design, in days. By default, it is set to 365.
#' @param pilevel The prediction interval level. By default,
#'   it is set to 0.90.
#' @param nyears The number of years after the data cut for prediction.
#'   By default, it is set to 4.
#' @param nreps The number of replications for simulation. By default,
#'   it is set to 500.
#' @param showEnrollment A Boolean variable to control whether or not to
#'   show the number of enrolled subjects. By default, it is set to
#'   \code{TRUE}.
#' @param showEvent A Boolean variable to control whether or not to
#'   show the number of events. By default, it is set to
#'   \code{TRUE}.
#' @param showDropout A Boolean variable to control whether or not to
#'   show the number of dropouts. By default, it is set to
#'   \code{FALSE}.
#' @param showOngoing A Boolean variable to control whether or not to
#'   show the number of ongoing subjects. By default, it is set to
#'   \code{FALSE}.
#' @param showsummary A Boolean variable to control whether or not to
#'   show the prediction summary. By default, it is set to \code{TRUE}.
#' @param showplot A Boolean variable to control whether or not to
#'   show the plots. By default, it is set to \code{TRUE}.
#'
#' @details
#' For the time-decay model, the mean function is
#' \code{mu(t) = mu/delta*(t - 1/delta*(1 - exp(-delta*t)))}
#' and the rate function is
#' \code{lambda(t) = mu/delta*(1 - exp(-delta*t))}.
#' For the B-spline model, the daily enrollment rate is approximated as
#' \code{lambda(t) = exp(B(t)*theta)},
#' where \code{B(t)} represents the B-spline basis functions.
#'
#' The \code{enroll_model_parameter} variable can be used for
#' enrollment prediction at the design stage. A piecewise Poisson
#' can be parameterized through the time intervals,
#' \code{accrualTime}, which is treated as fixed, and the
#' enrollment rates in the intervals, \code{accrualIntensity},
#' the log of which is used as the model parameter. For the
#' homogeneous Poisson, time-decay, and piecewise Poisson models,
#' \code{enroll_model_parameter} is used to specify the prior
#' distribution of model parameters, with a very small variance
#' being used to fix the parameter values. It should be noted
#' that the B-spline model is not appropriate for use during
#' the design stage.
#'
#' The \code{event_model_parameter} variable should be a list that
#' includes \code{model} to specify the event process,
#' \code{ngroups} to indicate the number of treatment groups,
#' \code{alloc} to indicate the treatment allocation within a
#' randomization block, \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix, both of which
#' have \code{ngroups} blocks with the \code{j}-th block specifying
#' the prior distribution of model parameters for the \code{j}-th
#' treatment group. For the piecewise exponential event model,
#' this should also include \code{piecewiseSurvivalTime} to indicate
#' the location of knots. It should be noted that the model averaging
#' option is not appropriate for use during the design stage.
#'
#' The \code{dropout_model_parameter} should be a list that
#' includes \code{model} to specify the dropout process,
#' \code{ngroups} to indicate the number of treatment groups,
#' \code{alloc} to indicate the treatment allocation within a
#' randomization block, \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix, both of which
#' have \code{ngroups} blocks with the \code{j}-th block specifying
#' the prior distribution of model parameters for the \code{j}-th
#' treatment group. For the piecewise exponential dropout model,
#' this should also include \code{piecewiseDropoutTime} to indicate
#' the location of knots.
#'
#' For analysis-stage enrollment and event prediction, the
#' \code{enroll_model_parameter}, \code{event_model_parameter}, and
#' \code{dropout_model_parameter} are either set to \code{NULL} to
#' use the observed data only, or specify the prior distribution
#' of model parameters to be combined with observed data likelihood
#' for enhanced modeling flexibility.
#'
#' @return A list that includes the fits of observed data models,
#' as well as simulated enrollment data for new subjects and
#' simulated event data for ongoing and new subjects.
#'
#' @examples
#'
#' # Event prediction after enrollment completion
#'
#' pred <- getPrediction(
#'   df = interimData2, to_predict = "event only",
#'   target_d = 200,
#'   event_model = "weibull",
#'   dropout_model = "exponential",
#'   pilevel = 0.90, nreps = 100)
#'
#' @export
#'
getPrediction <- function(
    df = NULL, to_predict = "enrollment and event",
    target_n = NA, target_d = NA,
    enroll_model = "b-spline", nknots = 0, lags = 30, accrualTime = 0,
    enroll_model_parameter = NULL,
    event_model = "model averaging", piecewiseSurvivalTime = 0,
    event_model_parameter = NULL,
    dropout_model = "exponential", piecewiseDropoutTime = 0,
    dropout_model_parameter = NULL,
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nyears = 4, nreps = 500,
    showEnrollment = TRUE, showEvent = TRUE,
    showDropout = FALSE, showOngoing = FALSE,
    showsummary = TRUE, showplot = TRUE) {

  if (!is.null(df)) erify::check_class(df, "data.frame")

  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))

  if (!is.na(target_n)) erify::check_n(target_n)
  if (!is.na(target_d)) erify::check_n(target_d)
  if (is.na(target_n) && is.na(target_d))
    stop("At least one of target_n and target_d must be specified.")
  if (!is.na(target_n) && !is.na(target_d) && target_d > target_n)
    stop("target_d cannot exceed target_n.")


  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay", "b-spline",
                         "piecewise poisson"))

  erify::check_n(nknots, zero = TRUE)
  erify::check_n(lags, zero = TRUE)

  if (accrualTime[1] != 0) {
    stop("accrualTime must start with 0")
  }
  if (length(accrualTime) > 1 && any(diff(accrualTime) <= 0)) {
    stop("accrualTime should be increasing")
  }

  if (!is.null(enroll_model_parameter)) {
    erify::check_class(enroll_model_parameter, "list")
    erify::check_content(tolower(enroll_model_parameter$model),
                         c("poisson", "time-decay", "piecewise poisson"))

    if (!is.null(df)) {
      if (tolower(enroll_model_parameter$model) != tolower(enroll_model)) {
        stop("Prior and likelihood must use the same enrollment model.")
      }

      if (tolower(enroll_model_parameter$model) == "piecewise poisson" &&
          !all.equal(enroll_model_parameter$accrualTime, accrualTime)) {
        stop("Knots for piecewise Poisson must be the same for prior.")
      }
    }
  }


  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-normal",
                         "piecewise exponential", "model averaging"))

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }
  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (!is.null(event_model_parameter)) {
    erify::check_class(event_model_parameter, "list")
    erify::check_content(tolower(event_model_parameter$model),
                         c("exponential", "weibull", "log-normal",
                           "piecewise exponential"))

    if (length(event_model_parameter$alloc) !=
        event_model_parameter$ngroups) {
      stop(paste("Number of treatments in alloc must be equal to ngroups",
                 "in event_model_parameter"))
    }

    if (nrow(event_model_parameter$vtheta) !=
        length(event_model_parameter$theta) ||
        ncol(event_model_parameter$vtheta) !=
        length(event_model_parameter$theta)) {
      stop(paste("Dimensions of vtheta must be compatible with the length",
                 "of theta in event_model_parameter"))
    }

    if ((tolower(event_model_parameter$model) == "exponential" &&
         length(event_model_parameter$theta) !=
         event_model_parameter$ngroups) ||
        (tolower(event_model_parameter$model) == "weibull" &&
         length(event_model_parameter$theta) !=
         2*event_model_parameter$ngroups) ||
        (tolower(event_model_parameter$model) == "log-normal" &&
         length(event_model_parameter$theta) !=
         2*event_model_parameter$ngroups) ||
        (tolower(event_model_parameter$model) == "piecewise_exponential" &&
         length(event_model_parameter$theta) !=
         length(event_model_parameter$piecewiseSurvivalTime)*
         event_model_parameter$ngroups)) {
      stop(paste("Length of theta must be compatible with ngroups",
                 "in event_model_parameter"))
    }

    if (tolower(event_model_parameter$model) == "piecewise_exponential") {
      if (event_model_parameter$piecewiseSurvivalTime[1] != 0) {
        stop(paste("piecewiseSurvivalTime must start with 0",
                   "in event_model_parameter"))
      }
      if (length(event_model_parameter$piecewiseSurvivalTime) > 1 &&
          any(diff(event_model_parameter$piecewiseSurvivalTime) <= 0)) {
        stop(paste("piecewiseSurvivalTime should be increasing",
                   "in event_model_parameter"))
      }
    }

    if (!is.null(df)) {
      if (tolower(event_model_parameter$model) != tolower(event_model)) {
        stop("Prior and likelihood must use the same event model.")
      }

      if (tolower(event_model_parameter$model) == "piecewise exponential"
          && !all.equal(event_model_parameter$piecewiseSurvivalTime,
                        piecewiseSurvivalTime)) {
        stop(paste("Knots for piecewise exponential survival must be",
        "the same for prior."))
      }
    }
  }

  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-normal",
                         "piecewise exponential"))

  if (piecewiseDropoutTime[1] != 0) {
    stop("piecewiseDropoutTime must start with 0")
  }
  if (length(piecewiseDropoutTime) > 1 &&
      any(diff(piecewiseDropoutTime) <= 0)) {
    stop("piecewiseDropoutTime should be increasing")
  }


  if (!is.null(dropout_model_parameter)) {
    erify::check_class(dropout_model_parameter, "list")
    erify::check_content(tolower(dropout_model_parameter$model),
                         c("exponential", "weibull", "log-normal",
                           "piecewise exponential"))

    if (length(dropout_model_parameter$alloc) !=
        dropout_model_parameter$ngroups) {
      stop(paste("Number of treatments in alloc must be equal to ngroups",
                 "in dropout_model_parameter"))
    }

    if (nrow(dropout_model_parameter$vtheta) !=
        length(dropout_model_parameter$theta) ||
        ncol(dropout_model_parameter$vtheta) !=
        length(dropout_model_parameter$theta)) {
      stop(paste("Dimensions of vtheta must be compatible with the length",
                 "of theta in dropout_model_parameter"))
    }

    if ((tolower(dropout_model_parameter$model) == "exponential" &&
         length(dropout_model_parameter$theta) !=
         dropout_model_parameter$ngroups) ||
        (tolower(dropout_model_parameter$model) == "weibull" &&
         length(dropout_model_parameter$theta) !=
         2*dropout_model_parameter$ngroups) ||
        (tolower(dropout_model_parameter$model) == "log-normal" &&
         length(dropout_model_parameter$theta) !=
         2*dropout_model_parameter$ngroups) ||
        (tolower(dropout_model_parameter$model) == "piecewise_exponential" &&
         length(dropout_model_parameter$theta) !=
         length(dropout_model_parameter$piecewiseDropoutTime)*
         dropout_model_parameter$ngroups)) {
      stop(paste("Length of theta must be compatible with ngroups",
                 "in dropout_model_parameter"))
    }

    if (tolower(dropout_model_parameter$model) == "piecewise_exponential") {
      if (dropout_model_parameter$piecewiseDropoutTime[1] != 0) {
        stop(paste("piecewiseDropoutTime must start with 0",
                   "in dropout_model_parameter"))
      }
      if (length(dropout_model_parameter$piecewiseDropoutTime) > 1 &&
          any(diff(dropout_model_parameter$piecewiseDropoutTime) <= 0)) {
        stop(paste("piecewiseDropoutTime should be increasing",
                   "in dropout_model_parameter"))
      }
    }


    if (!is.null(df)) {
      if (tolower(dropout_model_parameter$model) != tolower(dropout_model)) {
        stop("Prior and likelihood must use the same dropout model.")
      }

      if (tolower(dropout_model_parameter$model) == "piecewise exponential"
          && !all.equal(dropout_model_parameter$piecewiseDropoutTime,
                        piecewiseDropoutTime)) {
        stop(paste("Knots for piecewise exponential dropout must be",
                   "the same for prior."))
      }
    }


    if (!is.null(event_model_parameter)) {
      if (event_model_parameter$ngroups != dropout_model_parameter$ngroups) {
        stop(paste("Number of treatments must match between",
                   "event_model_parameter and dropout_model_parameter"))
      }

      if (!all.equal(event_model_parameter$alloc,
                     dropout_model_parameter$alloc)) {
        stop(paste("Treatment allocation must match between",
                   "event_model_parameter and dropout_model_parameter"))
      }
    }
  }


  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_positive(nyears)
  erify::check_n(nreps)
  erify::check_bool(showEnrollment)
  erify::check_bool(showEvent)
  erify::check_bool(showDropout)
  erify::check_bool(showOngoing)
  erify::check_bool(showsummary)
  erify::check_bool(showplot)


  if (!is.null(df)) {
    df <- dplyr::as_tibble(df)
    names(df) <- tolower(names(df))
    df$trialsdt <- as.Date(df$trialsdt)
    df$randdt <- as.Date(df$randdt)
    df$cutoffdt <- as.Date(df$cutoffdt)

    trialsdt = df$trialsdt[1]
    cutoffdt = df$cutoffdt[1]

    # summarize observed data
    observed <- summarizeObserved(df, to_predict, showplot)
  }


  # fit and predict enrollment
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) {
      erify::check_n(target_n - observed$n0,
                     supplement = "Enrollment target reached.")

      enroll_fit <- fitEnrollment(df = observed$adsl, enroll_model,
                                  nknots, accrualTime, showplot)
      enroll_fit1 <- enroll_fit$enroll_fit

      # combine prior and likelihood to yield posterior
      if (!is.null(enroll_model_parameter)) {
        if (length(enroll_model$theta) == 1) {
          enroll_fit1$theta <-
            1/(1/enroll_fit1$vtheta +
                 1/enroll_model_parameter$vtheta)*
            (1/enroll_fit1$vtheta*enroll_fit1$theta +
               1/enroll_model_parameter$vtheta*enroll_model_parameter$theta)
          enroll_fit1$vtheta <-
            1/(1/enroll_fit1$vtheta +
                 1/enroll_model_parameter$vtheta)
        } else {
          enroll_fit1$theta <-
            solve(solve(enroll_fit1$vtheta) +
                    solve(enroll_model_parameter$vtheta),
                  solve(enroll_fit1$vtheta, enroll_fit1$theta) +
                    solve(enroll_model_parameter$vtheta,
                          enroll_model_parameter$theta))
          enroll_fit1$vtheta <- solve(solve(enroll_fit1$vtheta) +
                                        solve(enroll_model_parameter$vtheta))
        }
      }

      # enrollment prediction at the analysis stage
      enroll_pred <- predictEnrollment(
        df = observed$adsl, target_n,
        enroll_fit = enroll_fit1,
        lags, pilevel, nyears, nreps,
        showsummary, showplot = FALSE)
    } else {
      # enrollment prediction at the design stage
      enroll_pred <- predictEnrollment(
        df = NULL, target_n,
        enroll_fit = enroll_model_parameter,
        lags, pilevel, nyears, nreps,
        showsummary, showplot = FALSE)
    }
  }


  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) { # event prediction at analysis stage
      erify::check_n(target_d - observed$d0,
                     supplement = "Event target reached.")

      # convert prior by treatment to prior overall
      if (!is.null(event_model_parameter)) {
        k = event_model_parameter$ngroups
        alloc = event_model_parameter$alloc
        w = alloc/sum(alloc)

        if (tolower(event_model_parameter$model) == "exponential") {
          # match the overall mean
          lambda = exp(event_model_parameter$theta)
          lambda1 = 1/sum(w/lambda)  # hazard rate for pooled
          theta1 = log(lambda1)

          # use delta-method to obtain the variance
          vtheta1 = 0
          for (i in 1:k) {
            vtheta1 = vtheta1 + (w[i]/lambda[i])^2 *
              event_model_parameter$vtheta[i,i]
          }
          vtheta1 = vtheta1*lambda1^2

          event_prior <- list(
            model = event_model_parameter$model,
            theta = theta1,
            vtheta = vtheta1)
        } else if (tolower(event_model_parameter$model) == "weibull") {
          # match the overall mean and variance

          # mean and variance of weibull as a function of theta
          fmweibull <- function(theta) {
            k = length(theta)/2
            shape = exp(theta[2*(1:k)-1])
            scale = exp(theta[2*(1:k)])
            list(mean = scale*gamma(1+1/shape),
                 var = scale^2*(gamma(1+2/shape) - (gamma(1+1/shape))^2))
          }

          # gradient vector
          gmweibull <- function(theta) {
            g1 = numDeriv::grad(function(theta) fmweibull(theta)$mean, theta)
            g2 = numDeriv::grad(function(theta) fmweibull(theta)$var, theta)
            matrix(c(g1, g2), nrow=2, byrow=TRUE)
          }

          # mean and variance by treatment group
          m1 = fmweibull(event_model_parameter$theta)

          # mean and variance for pooled
          m2 = list(mean = sum(w*m1$mean),
                    var = sum(w*m1$var) +
                      sum(w*m1$mean^2) - (sum(w*m1$mean))^2)

          # solve for theta given the mean and variance for pooled
          theta11 = uniroot(function(x)
            lgamma(1+2/exp(x)) - 2*lgamma(1+1/exp(x)) -
              log(m2$var/m2$mean^2 + 1),
            c(min(event_model_parameter$theta[2*(1:k)-1]) - 1,
              max(event_model_parameter$theta[2*(1:k)-1]) + 1),
            extendInt = "yes")$root

          theta12 = log(m2$mean) - lgamma(1+1/exp(theta11))
          theta1 = c(theta11, theta12)

          # gradient of theta with respect to mean and variance for pooled
          ig = solve(gmweibull(theta1))

          # variance of theta for pooled
          vtheta1 = 0
          for (i in 1:k) {
            index = c(2*i-1, 2*i)
            gi = gmweibull(event_model_parameter$theta[index])
            vm1i = gi * event_model_parameter$vtheta[index,index] * t(gi)
            li = w[i]*matrix(c(1, 2*(m1$mean[i] - m2$mean), 0, 1), ncol=2)
            vtheta1 = vtheta1 + li %*% vm1i%*% t(li)
          }
          vtheta1 = ig %*% vtheta1 %*% t(ig)

          event_prior <- list(
            model = event_model_parameter$model,
            theta = theta1,
            vtheta = vtheta1)
        } else if (tolower(event_model_parameter$model) == "log-normal") {
          # match the overall mean and variance

          # mean and variance of log-normal as a function of theta
          fmlnorm <- function(theta) {
            k = length(theta)/2
            meanlog = theta[2*(1:k)-1]
            sdlog = exp(theta[2*(1:k)])
            list(mean = exp(meanlog + sdlog^2/2),
                 var = (exp(sdlog^2) - 1)*exp(2*meanlog + sdlog^2))
          }

          # gradient vector
          gmlnorm <- function(theta) {
            g1 = numDeriv::grad(function(theta) fmlnorm(theta)$mean, theta)
            g2 = numDeriv::grad(function(theta) fmlnorm(theta)$var, theta)
            matrix(c(g1, g2), nrow=2, byrow=TRUE)
          }

          # mean and variance by treatment group
          m1 = fmlnorm(event_model_parameter$theta)

          # mean and variance for pooled
          m2 = list(mean = sum(w*m1$mean),
                    var = sum(w*m1$var) +
                      sum(w*m1$mean^2) - (sum(w*m1$mean))^2)

          # solve for theta given the mean and variance for pooled
          theta12 = 0.5*log(log(m2$var/m2$mean^2 + 1))
          theta11 = log(m2$mean) - 0.5*exp(2*theta12)
          theta1 = c(theta11, theta12)

          # gradient of theta with respect to mean and variance for pooled
          ig = solve(gmlnorm(theta1))

          # variance of theta for pooled
          vtheta1 = 0
          for (i in 1:k) {
            index = c(2*i-1, 2*i)
            gi = gmlnorm(event_model_parameter$theta[index])
            vm1i = gi * event_model_parameter$vtheta[index,index] * t(gi)
            li = w[i]*matrix(c(1, 2*(m1$mean[i] - m2$mean), 0, 1), ncol=2)
            vtheta1 = vtheta1 + li %*% vm1i%*% t(li)
          }
          vtheta1 = ig %*% vtheta1 %*% t(ig)

          event_prior <- list(
            model = event_model_parameter$model,
            theta = theta1,
            vtheta = vtheta1)
        } else if (tolower(event_model_parameter$model) ==
                   "piecewise exponential") {
          # match the mean within each interval
          lambda = exp(event_model_parameter$theta)
          npieces = length(event_model_parameter$piecewiseSurvivalTime)

          # construct theta and vtheta piece by piece
          theta1 = rep(NA, npieces)
          vtheta1 = 0*diag(npieces)
          for (j in 1:npieces) {
            lambdaj = lambda[seq(0, k-1)*npieces + j]
            lambda1j = 1/sum(w/lambdaj)
            theta1[j] = log(lambda1j)
            for (i in 1:k) {
              vtheta1[j,j] = vtheta1[j,j] + (w[i]/lambdaj[i])^2 *
                event_model_parameter$vtheta[(i-1)*npieces + j,
                                             (i-1)*npieces + j]
            }
            vtheta1[j,j] = vtheta1[j,j]*lambda1j^2
          }

          event_prior <- list(
            model = event_model_parameter$model,
            theta = theta1,
            vtheta = vtheta1,
            piecewiseSurvivalTime =
              event_model_parameter$piecewiseSurvivalTime)
        }
      }

      if (sum(observed$adtte$event == 1) > 0) {
        event_fit <- fitEvent(df = observed$adtte, event_model,
                              piecewiseSurvivalTime, showplot)
        event_fit1 <- event_fit$event_fit
      } else {
        if (is.null(event_model_parameter)) {
          stop("Prior must be specified if there is no event observed.")
        }
        event_fit <- list()
        event_fit1 <- event_prior
        event_fit1$vtheta = event_prior$vtheta*1e8
        event_fit1$bic = NA
      }

      # combine prior and likelihood to yield posterior
      if (!is.null(event_model_parameter)) {
        if (length(event_fit1$theta) == 1) {
          event_fit1$theta <- 1/(1/event_fit1$vtheta + 1/event_prior$vtheta)*
            (1/event_fit1$vtheta*event_fit1$theta +
               1/event_prior$vtheta*event_prior$theta)
          event_fit1$vtheta <- 1/(1/event_fit1$vtheta + 1/event_prior$vtheta)
        } else {
          event_fit1$theta <-
            solve(solve(event_fit1$vtheta) + solve(event_prior$vtheta),
                  solve(event_fit1$vtheta, event_fit1$theta) +
                    solve(event_prior$vtheta, event_prior$theta))
          event_fit1$vtheta <-
            solve(solve(event_fit1$vtheta) + solve(event_prior$vtheta))
        }
      }


      # whether to include dropout model
      if (tolower(dropout_model) != "none") {

        # convert prior by treatment to prior overall
        if (!is.null(dropout_model_parameter)) {
          k = dropout_model_parameter$ngroups
          alloc = dropout_model_parameter$alloc

          w = alloc/sum(alloc)

          if (tolower(dropout_model_parameter$model) == "exponential") {
            # match the overall mean
            lambda = exp(dropout_model_parameter$theta)
            lambda1 = 1/sum(w/lambda)  # hazard rate for pooled
            theta1 = log(lambda1)

            # use delta-method to obtain the variance
            vtheta1 = 0
            for (i in 1:k) {
              vtheta1 = vtheta1 + (w[i]/lambda[i])^2 *
                dropout_model_parameter$vtheta[i,i]
            }
            vtheta1 = vtheta1*lambda1^2

            dropout_prior <- list(
              model = dropout_model_parameter$model,
              theta = theta1,
              vtheta = vtheta1)
            } else if (tolower(dropout_model_parameter$model) == "weibull") {
            # match the overall mean and variance

            # mean and variance of weibull as a function of theta
            fmweibull <- function(theta) {
              k = length(theta)/2
              shape = exp(theta[2*(1:k)-1])
              scale = exp(theta[2*(1:k)])
              list(mean = scale*gamma(1+1/shape),
                   var = scale^2*(gamma(1+2/shape) - (gamma(1+1/shape))^2))
            }

            # gradient vector
            gmweibull <- function(theta) {
              g1 = numDeriv::grad(function(theta) fmweibull(theta)$mean, theta)
              g2 = numDeriv::grad(function(theta) fmweibull(theta)$var, theta)
              matrix(c(g1, g2), nrow=2, byrow=TRUE)
            }

            # mean and variance by treatment group
            m1 = fmweibull(dropout_model_parameter$theta)

            # mean and variance for pooled
            m2 = list(mean = sum(w*m1$mean),
                      var = sum(w*m1$var) +
                        sum(w*m1$mean^2) - (sum(w*m1$mean))^2)

            # solve for theta given the mean and variance for pooled
            theta11 = uniroot(function(x)
              lgamma(1+2/exp(x)) - 2*lgamma(1+1/exp(x)) -
                log(m2$var/m2$mean^2 + 1),
              c(min(dropout_model_parameter$theta[2*(1:k)-1]) - 1,
                max(dropout_model_parameter$theta[2*(1:k)-1]) + 1),
              extendInt = "yes")$root

            theta12 = log(m2$mean) - lgamma(1+1/exp(theta11))
            theta1 = c(theta11, theta12)

            # gradient of theta with respect to mean and variance for pooled
            ig = solve(gmweibull(theta1))

            # variance of theta for pooled
            vtheta1 = 0
            for (i in 1:k) {
              index = c(2*i-1, 2*i)
              gi = gmweibull(dropout_model_parameter$theta[index])
              vm1i = gi * dropout_model_parameter$vtheta[index,index] * t(gi)
              li = w[i]*matrix(c(1, 2*(m1$mean[i] - m2$mean), 0, 1), ncol=2)
              vtheta1 = vtheta1 + li %*% vm1i%*% t(li)
            }
            vtheta1 = ig %*% vtheta1 %*% t(ig)

            dropout_prior <- list(
              model = dropout_model_parameter$model,
              theta = theta1,
              vtheta = vtheta1)
          } else if (tolower(dropout_model_parameter$model) == "log-normal") {
            # match the overall mean and variance

            # mean and variance of log-normal as a function of theta
            fmlnorm <- function(theta) {
              k = length(theta)/2
              meanlog = theta[2*(1:k)-1]
              sdlog = exp(theta[2*(1:k)])
              list(mean = exp(meanlog + sdlog^2/2),
                   var = (exp(sdlog^2) - 1)*exp(2*meanlog + sdlog^2))
            }

            # gradient vector
            gmlnorm <- function(theta) {
              g1 = numDeriv::grad(function(theta) fmlnorm(theta)$mean, theta)
              g2 = numDeriv::grad(function(theta) fmlnorm(theta)$var, theta)
              matrix(c(g1, g2), nrow=2, byrow=TRUE)
            }

            # mean and variance by treatment group
            m1 = fmlnorm(dropout_model_parameter$theta)

            # mean and variance for pooled
            m2 = list(mean = sum(w*m1$mean),
                      var = sum(w*m1$var) +
                        sum(w*m1$mean^2) - (sum(w*m1$mean))^2)

            # solve for theta given the mean and variance for pooled
            theta12 = 0.5*log(log(m2$var/m2$mean^2 + 1))
            theta11 = log(m2$mean) - 0.5*exp(2*theta12)
            theta1 = c(theta11, theta12)

            # gradient of theta with respect to mean and variance for pooled
            ig = solve(gmlnorm(theta1))

            # variance of theta for pooled
            vtheta1 = 0
            for (i in 1:k) {
              index = c(2*i-1, 2*i)
              gi = gmlnorm(dropout_model_parameter$theta[index])
              vm1i = gi * dropout_model_parameter$vtheta[index,index] * t(gi)
              li = w[i]*matrix(c(1, 2*(m1$mean[i] - m2$mean), 0, 1), ncol=2)
              vtheta1 = vtheta1 + li %*% vm1i%*% t(li)
            }
            vtheta1 = ig %*% vtheta1 %*% t(ig)

            dropout_prior <- list(
              model = dropout_model_parameter$model,
              theta = theta1,
              vtheta = vtheta1)
          } else if (tolower(dropout_model_parameter$model) ==
                     "piecewise exponential") {
            # match the mean within each interval
            lambda = exp(dropout_model_parameter$theta)
            npieces = length(dropout_model_parameter$piecewiseDropoutTime)

            # construct theta and vtheta piece by piece
            theta1 = rep(NA, npieces)
            vtheta1 = 0*diag(npieces)
            for (j in 1:npieces) {
              lambdaj = lambda[seq(0, k-1)*npieces + j]
              lambda1j = 1/sum(w/lambdaj)
              theta1[j] = log(lambda1j)
              for (i in 1:k) {
                vtheta1[j,j] = vtheta1[j,j] + (w[i]/lambdaj[i])^2 *
                  dropout_model_parameter$vtheta[(i-1)*npieces + j,
                                               (i-1)*npieces + j]
              }
              vtheta1[j,j] = vtheta1[j,j]*lambda1j^2
            }

            dropout_prior <- list(
              model = dropout_model_parameter$model,
              theta = theta1,
              vtheta = vtheta1,
              piecewiseDropoutTime =
                dropout_model_parameter$piecewiseDropoutTime)
          }
        }


        if (sum(observed$adtte$dropout == 1) > 0) {
          dropout_fit <- fitDropout(df = observed$adtte, dropout_model,
                                    piecewiseDropoutTime, showplot)
          dropout_fit1 <- dropout_fit$dropout_fit
        } else {
          if (is.null(dropout_model_parameter)) {
            stop("Prior must be specified if there is no dropout observed.")
          }
          dropout_fit <- list()
          dropout_fit1 <- dropout_prior
          dropout_fit1$vtheta <- dropout_prior$vtheta*1e8
          dropout_fit1$bic = NA
        }

        # combine prior and likelihood to yield posterior
        if (!is.null(dropout_model_parameter)) {
          if (length(dropout_fit1$theta) == 1) {
            dropout_fit1$theta <-
              1/(1/dropout_fit1$vtheta + 1/dropout_prior$vtheta)*
              (1/dropout_fit1$vtheta*dropout_fit1$theta +
                 1/dropout_prior$vtheta*dropout_prior$theta)
            dropout_fit1$vtheta <-
              1/(1/dropout_fit1$vtheta + 1/dropout_prior$vtheta)
          } else {
            dropout_fit1$theta <-
              solve(solve(dropout_fit1$vtheta) + solve(dropout_prior$vtheta),
                    solve(dropout_fit1$vtheta, dropout_fit1$theta) +
                      solve(dropout_prior$vtheta, dropout_prior$theta))
            dropout_fit1$vtheta <-
              solve(solve(dropout_fit1$vtheta) + solve(dropout_prior$vtheta))
          }
        }


        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1,
            dropout_fit = dropout_fit1,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            dropout_fit = dropout_fit1,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE)
        }
      } else {
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1,
            dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE)
        }
      }
    } else { # event prediction at design stage
      if (!is.null(dropout_model_parameter)) {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_model_parameter,
          dropout_fit = dropout_model_parameter,
          fixedFollowup, followupTime, pilevel, nyears, nreps,
          showEnrollment, showEvent, showDropout, showOngoing,
          showsummary, showplot = FALSE)
      } else {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_model_parameter,
          dropout_fit = NULL,
          fixedFollowup, followupTime, pilevel, nyears, nreps,
          showEnrollment, showEvent, showDropout, showOngoing,
          showsummary, showplot = FALSE)
      }
    }
  }


  # output results
  if (is.null(df)) { # design stage prediction
    if (tolower(to_predict) == "enrollment only") {
      if (showplot) print(enroll_pred$enroll_pred_plot)

      list(stage = "Design stage",
           to_predict = "Enrollment only",
           enroll_fit = enroll_model_parameter, enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "enrollment and event") {
      if (showplot) print(event_pred$event_pred_plot)

      if (!is.null(dropout_model_parameter)) {
        list(stage = "Design stage",
             to_predict = "Enrollment and event",
             enroll_fit = enroll_model_parameter, enroll_pred = enroll_pred,
             event_fit = event_model_parameter,
             dropout_fit = dropout_model_parameter, event_pred = event_pred)
      } else {
        list(stage = "Design stage",
             to_predict = "Enrollment and event",
             enroll_fit = enroll_model_parameter, enroll_pred = enroll_pred,
             event_fit = event_model_parameter, event_pred = event_pred)
      }
    }
  } else { # analysis stage prediction
    if (tolower(to_predict) == "enrollment only") {
      if (showplot) print(enroll_pred$enroll_pred_plot)

      list(stage = "Real-time before enrollment completion",
           to_predict = "Enrollment only",
           observed = observed, enroll_fit = enroll_fit,
           enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "enrollment and event") {
      if (showplot) print(event_pred$event_pred_plot)

      if (tolower(dropout_model) != "none") {
        list(stage = "Real-time before enrollment completion",
             to_predict = "Enrollment and event",
             observed = observed, enroll_fit = enroll_fit,
             enroll_pred = enroll_pred, event_fit = event_fit,
             dropout_fit = dropout_fit, event_pred = event_pred)
      } else {
        list(stage = "Real-time before enrollment completion",
             to_predict = "Enrollment and event",
             observed = observed, enroll_fit = enroll_fit,
             enroll_pred = enroll_pred, event_fit = event_fit,
             event_pred = event_pred)
      }
    } else if (tolower(to_predict) == "event only") {
      if (showplot) print(event_pred$event_pred_plot)

      if (tolower(dropout_model) != "none") {
        list(stage = "Real-time after enrollment completion",
             to_predict = "Event only",
             observed = observed, event_fit = event_fit,
             dropout_fit = dropout_fit, event_pred = event_pred)
      } else {
        list(stage = "Real-time after enrollment completion",
             to_predict = "Event only",
             observed = observed, event_fit = event_fit,
             event_pred = event_pred)
      }
    }
  }
}


