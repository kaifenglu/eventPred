#' @title Enrollment and event prediction
#' @description Performs enrollment and event prediction by utilizing
#'   observed data and specified enrollment and event models.
#'
#' @param df The subject-level enrollment and event data,
#'   including \code{randdt}, \code{cutoffdt}, \code{time}, \code{event},
#'   and \code{dropout}. By default, it is set to \code{NULL} for
#'   enrollment and event prediction at the design stage.
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
#'   Poisson model. Must start with 0, e.g., c(0, 3) breaks the
#'   time axis into 2 accrual intervals: [0, 3) and [3, Inf).
#'   By default, it is set to 0.
#' @param enroll_model_parameter The enrollment model parameters for
#'   design-stage enrollment prediction.
#' @param event_model The event model used to analyze the event data
#'   which can be set to one of the
#'   following options: "exponential", "Weibull", "log-normal",
#'   "piecewise exponential", or "model averaging". The model averaging
#'   uses the \code{exp(-bic)} weighting and combines Weibull and
#'   log-normal models. By default, it is set to "model
#'   averaging".
#' @param piecewiseSurvivalTime A vector that specifies the time
#'   intervals for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 6) breaks the time axis into 2
#'   event intervals: [0, 6) and [6, Inf). By default, it is set to 0.
#' @param event_model_parameter The event model parameters for
#'   design-stage event prediction.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options: "exponential",
#'   "Weibull", or "log-normal". By default, it is set to "Weibull".
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
#' @param showplot A Boolean variable to control whether or not to
#'   show the plots. By default, it is set to \code{FALSE}.
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
#' \code{prob} to indicate the randomization probabilities
#' for each group, \code{theta} and \code{vtheta} to indicate
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
#' \code{prob} to indicate the randomization probabilities
#' for each group, \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix, both of which
#' have \code{ngroups} blocks with the \code{j}-th block specifying
#' the prior distribution of model parameters for the \code{j}-th
#' treatment group.
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
    dropout_model = "weibull",
    dropout_model_parameter = NULL,
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nyears = 4, nreps = 500,
    showEnrollment = TRUE, showEvent = TRUE,
    showDropout = FALSE, showOngoing = FALSE,
    showplot = TRUE) {

  if (!is.null(df)) erify::check_class(df, "data.frame")

  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))

  if (!is.na(target_n)) erify::check_n(target_n)
  if (!is.na(target_d)) erify::check_n(target_d)
  if (is.na(target_n) & is.na(target_d))
    stop("At least one of target_n and target_d must be specified.")

  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay", "b-spline",
                         "piecewise poisson"))

  erify::check_n(nknots, zero = TRUE)

  if (accrualTime[1] != 0) {
    stop("accrualTime must start with 0");
  }
  if (length(accrualTime) > 1 & any(diff(accrualTime) <= 0)) {
    stop("accrualTime should be increasing")
  }

  if (!is.null(enroll_model_parameter)) {
    erify::check_class(enroll_model_parameter, "list")
    if (!is.null(df)) {
      if (tolower(enroll_model_parameter$model) != tolower(enroll_model)) {
        stop("Prior and likelihood must use the same enrollment model.")
      }

      if (tolower(enroll_model_parameter$model) == "piecewise poisson" &
          !all.equal(enroll_model_parameter$accrualTime, accrualTime)) {
        stop("Knots for piecewise Poisson must be the same for prior.")
      }

      if (tolower(enroll_model_parameter$model) == "b-spline") {
        stop("B-spline enrollment model cannot be used as prior.")
      }
    }
  }

  erify::check_n(lags, zero=TRUE)

  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-normal",
                         "piecewise exponential", "model averaging"))

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  if (length(piecewiseSurvivalTime) > 1 &
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (!is.null(event_model_parameter)) {
    erify::check_class(event_model_parameter, "list")
    if (!is.null(df)) {
      if (tolower(event_model_parameter$model) != tolower(event_model)) {
        stop("Prior and likelihood must use the same event model.")
      }

      if (tolower(event_model_parameter$model) == "piecewise exponential" &
          !all.equal(event_model_parameter$piecewiseSurvivalTime,
                     piecewiseSurvivalTime)) {
        stop("Knots for piecewise exponential must be the same for prior.")
      }

      if (tolower(event_model_parameter$model) == "model averaging") {
        stop("Model averaging event model cannot be used as prior.")
      }
    }
  }

  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-normal"))

  if (!is.null(dropout_model_parameter)) {
    erify::check_class(dropout_model_parameter, "list")
    if (!is.null(df)) {
      if (tolower(dropout_model_parameter$model) != tolower(dropout_model)) {
        stop("Prior and likelihood must use the same dropout model.")
      }
      if (tolower(dropout_model_parameter$model) == "none") {
        stop("None dropout model cannot be used as prior.")
      }
    }
  }

  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_positive(nyears)
  erify::check_n(nreps)
  erify::check_bool(showOngoing)


  if (!is.null(df)) {
    df <- dplyr::as_tibble(df)
    names(df) <- tolower(names(df))
    df$randdt <- as.Date(df$randdt)
    df$cutoffdt <- as.Date(df$cutoffdt)

    trialsdt = min(df$randdt)
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
        if (tolower(enroll_model) == "poisson") {
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
        lags, pilevel, nyears, nreps, showplot = FALSE)
    } else {
      # enrollment prediction at the design stage
      enroll_pred <- predictEnrollment(
        df = NULL, target_n,
        enroll_fit = enroll_model_parameter,
        lags, pilevel, nyears, nreps, showplot = FALSE)
    }
  }

  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) { # event prediction at analysis stage
      erify::check_n(target_d - observed$d0,
                     supplement = "Event target reached.")

      # convert prior by treatment to prior overall
      if (!is.null(event_model_parameter)) {
        if (tolower(event_model_parameter$model) == "exponential") {
          # match the overall mean
          w = event_model_parameter$prob
          lambda = exp(event_model_parameter$theta)

          event_prior <- list(
            model = event_model_parameter$model,
            theta = log(1/sum(w/lambda)),
            vtheta = event_model_parameter$vtheta[1,1])
        } else if (tolower(event_model_parameter$model) == "weibull") {
          # match the overall mean and variance
          k = event_model_parameter$ngroups
          w = event_model_parameter$prob
          shape = exp(event_model_parameter$theta[2*(1:k)-1])
          scale = exp(event_model_parameter$theta[2*(1:k)])

          fmweibull <- function(shape, scale) {
            list(mean = scale*gamma(1+1/shape),
                 var = scale^2*(gamma(1+2/shape) - (gamma(1+1/shape))^2))
          }

          sol <- rootSolve::multiroot(function(x) {
            m1 = fmweibull(shape, scale)
            m2 = fmweibull(exp(x[1]), exp(x[2]))
            y1 = sum(w*m1$mean) - m2$mean
            y2 = sum(w*m1$var) + as.numeric(
              m1$mean %*% (diag(w) - w %*% t(w)) %*% m1$mean) - m2$var
            c(y1, y2)}, log(c(mean(shape), mean(scale))))$root

          event_prior <- list(
            model = event_model_parameter$model,
            theta = sol,
            vtheta = event_model_parameter$vtheta[1:2,1:2])
        } else if (tolower(event_model_parameter$model) == "log-normal") {
          # match the overall mean and variance
          k = event_model_parameter$ngroups
          w = event_model_parameter$prob
          meanlog = event_model_parameter$theta[2*(1:k)-1]
          sdlog = exp(event_model_parameter$theta[2*(1:k)])

          fmlnorm <- function(meanlog, sdlog) {
            list(mean = exp(meanlog + sdlog^2/2),
                 var = (exp(sdlog^2) - 1)*exp(2*meanlog + sdlog^2))
          }

          sol <- rootSolve::multiroot(function(x) {
            m1 = fmlnorm(meanlog, sdlog)
            m2 = fmweibull(x[1], exp(x[2]))
            y1 = sum(w*m1$mean) - m2$mean
            y2 = sum(w*m1$var) + as.numeric(
              m1$mean %*% (diag(w) - w %*% t(w)) %*% m1$mean) - m2$var
            c(y1, y2)}, c(mean(meanlog), mean(log(sdlog))))$root

          event_prior <- list(
            model = event_model_parameter$model,
            theta = sol,
            vtheta = event_model_parameter$vtheta[1:2,1:2])
        } else if (tolower(event_model_parameter$model) ==
                   "piecewise exponential") {
          # match within each interval
          k = event_model_parameter$ngroups
          w = event_model_parameter$prob
          lambda = exp(event_model_parameter$theta)
          npieces = length(lambda)/k

          event_prior <- list(
            model = event_model_parameter$model,
            theta = log(1/as.numeric(w %*% matrix(
              1/lambda, nrow=k, ncol=npieces, byrow=TRUE))),
            vtheta = event_model_parameter$vtheta[1:npieces,1:npieces],
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
          if (tolower(dropout_model_parameter$model) == "exponential") {
            # match the overall mean
            w = dropout_model_parameter$prob
            lambda = exp(dropout_model_parameter$theta)

            dropout_prior <- list(
              model = dropout_model_parameter$model,
              theta = log(1/sum(w/lambda)),
              vtheta = dropout_model_parameter$vtheta[1,1])
          } else if (tolower(dropout_model_parameter$model) == "weibull") {
            # match the overall mean and variance
            k = dropout_model_parameter$ngroups
            w = dropout_model_parameter$prob
            shape = exp(dropout_model_parameter$theta[2*(1:k)-1])
            scale = exp(dropout_model_parameter$theta[2*(1:k)])

            fmweibull <- function(shape, scale) {
              list(mean = scale*gamma(1+1/shape),
                   var = scale^2*(gamma(1+2/shape) - (gamma(1+1/shape))^2))
            }

            sol <- rootSolve::multiroot(function(x) {
              m1 = fmweibull(shape, scale)
              m2 = fmweibull(exp(x[1]), exp(x[2]))
              y1 = sum(w*m1$mean) - m2$mean
              y2 = sum(w*m1$var) + as.numeric(
                m1$mean %*% (diag(w) - w %*% t(w)) %*% m1$mean) - m2$var
              c(y1, y2)}, log(c(mean(shape), mean(scale))))$root

            dropout_prior <- list(
              model = dropout_model_parameter$model,
              theta = sol,
              vtheta = dropout_model_parameter$vtheta[1:2,1:2])
          } else if (tolower(dropout_model_parameter$model) == "log-normal") {
            # match the overall mean and variance
            k = dropout_model_parameter$ngroups
            w = dropout_model_parameter$prob
            meanlog = dropout_model_parameter$theta[2*(1:k)-1]
            sdlog = exp(dropout_model_parameter$theta[2*(1:k)])

            fmlnorm <- function(meanlog, sdlog) {
              list(mean = exp(meanlog + sdlog^2/2),
                   var = (exp(sdlog^2) - 1)*exp(2*meanlog + sdlog^2))
            }

            sol <- rootSolve::multiroot(function(x) {
              m1 = fmlnorm(meanlog, sdlog)
              m2 = fmweibull(x[1], exp(x[2]))
              y1 = sum(w*m1$mean) - m2$mean
              y2 = sum(w*m1$var) + as.numeric(
                m1$mean %*% (diag(w) - w %*% t(w)) %*% m1$mean) - m2$var
              c(y1, y2)}, c(mean(meanlog), mean(log(sdlog))))$root

            dropout_prior <- list(
              model = dropout_model_parameter$model,
              theta = sol,
              vtheta = dropout_model_parameter$vtheta[1:2,1:2])
          }
        }


        if (sum(observed$adtte$dropout == 1) > 0) {
          dropout_fit <- fitDropout(df = observed$adtte, dropout_model,
                                    showplot)
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
            showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            dropout_fit = dropout_fit1,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showplot = FALSE)
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
            showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showplot = FALSE)
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
          showplot = FALSE)
      } else {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_model_parameter,
          dropout_fit = NULL,
          fixedFollowup, followupTime, pilevel, nyears, nreps,
          showEnrollment, showEvent, showDropout, showOngoing,
          showplot = FALSE)
      }
    }
  }


  # output results
  if (!is.null(df)) { # analysis stage prediction
    if (tolower(to_predict) == "enrollment only") {
      if (showplot) print(enroll_pred$enroll_pred_plot)
      list(observed = observed, enroll_fit = enroll_fit,
           enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "event only") {
      if (showplot) print(event_pred$event_pred_plot)

      if (tolower(dropout_model) != "none") {
        list(observed = observed, event_fit = event_fit,
             dropout_fit = dropout_fit, event_pred = event_pred)
      } else {
        list(observed = observed, event_fit = event_fit,
             event_pred = event_pred)
      }
    } else if (tolower(to_predict) == "enrollment and event") {
      if (showplot) print(event_pred$event_pred_plot)

      if (tolower(dropout_model) != "none") {
        list(observed = observed, enroll_fit = enroll_fit,
             enroll_pred = enroll_pred, event_fit = event_fit,
             dropout_fit = dropout_fit, event_pred = event_pred)
      } else {
        list(observed = observed, enroll_fit = enroll_fit,
             enroll_pred = enroll_pred, event_fit = event_fit,
             event_pred = event_pred)
      }
    }
  } else { # design stage prediction
    if (tolower(to_predict) == "enrollment only") {
      if (showplot) print(enroll_pred$enroll_pred_plot)
      list(enroll_fit = enroll_model_parameter, enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "event only") {
      if (showplot) print(event_pred$event_pred_plot)

      if (!is.null(dropout_model_parameter)) {
        list(event_fit = event_model_parameter,
             dropout_fit = dropout_model_parameter, event_pred = event_pred)
      } else {
        list(event_fit = event_model_parameter, event_pred = event_pred)
      }
    } else if (tolower(to_predict) == "enrollment and event") {
      if (showplot) print(event_pred$event_pred_plot)

      if (!is.null(dropout_model_parameter)) {
        list(enroll_fit = enroll_model_parameter, enroll_pred = enroll_pred,
             event_fit = event_model_parameter,
             dropout_fit = dropout_model_parameter, event_pred = event_pred)
      } else {
        list(enroll_fit = enroll_model_parameter, enroll_pred = enroll_pred,
             event_fit = event_model_parameter, event_pred = event_pred)
      }
    }
  }

}


