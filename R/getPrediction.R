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
#' @param accrualTime The accrual time intervals for the piecewise
#'   Poisson model. Must start with 0, e.g., c(0, 3) breaks the
#'   time axis into 2 accrual intervals: [0, 3) and [3, Inf).
#'   By default, it is set to 0.
#' @param parameter_enroll_model The enrollment model parameters for
#'   design-stage enrollment prediction.
#' @param lags The day lags to compute the average enrollment rate to
#'   carry forward for the B-spline enrollment model. By default,
#'   it is set to 30.
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
#' @param parameter_event_model The event model parameters for
#'   design-stage event prediction.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options: "exponential",
#'   "Weibull", or "log-normal". By default, it is set to "Weibull".
#' @param parameter_dropout_model The dropout model parameters for
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
#' @param showDropout A Boolean variable to control whether or not to
#'   show the number of dropouts. By default, it is set to
#'   \code{FALSE}.
#' @param showOngoing A Boolean variable to control whether or not to
#'   show the number of ongoing subjects. By default, it is set to
#'   \code{FALSE}.
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
#' The \code{parameter_enroll_model} variable can be used for
#' enrollment prediction at the design stage. A piecewise Poisson
#' can be parameterized through the time intervals,
#' \code{accrualTime}, which is treated as fixed, and the
#' enrollment rates in the intervals, \code{accrualIntensity},
#' the log of which is used as the model parameter. For the
#' homogeneous Poisson, time-decay, and piecewise Poisson models,
#' \code{parameter_enroll_model} is used to specify the prior
#' distribution of model parameters, with a very small variance
#' being used to fix the parameter values. It should be noted
#' that the B-spline model is not appropriate for use during
#' the design stage.
#'
#' The \code{parameter_event_model} variable should be a list that
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
#' The \code{parameter_dropout_model} should be a list that
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
#' \code{parameter_enroll_model}, \code{parameter_event_model}, and
#' \code{parameter_dropout_model} are either set to \code{NULL} to
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
#'   target_d = 244,
#'   event_model = "weibull",
#'   dropout_model = "exponential",
#'   pilevel = 0.90, nreps = 200)
#'
#' @export
#'
getPrediction <- function(
    df = NULL, to_predict = "enrollment and event",
    target_n = NA, target_d = NA,
    enroll_model = "b-spline", nknots = 0, accrualTime = 0,
    parameter_enroll_model = NULL, lags = 30,
    event_model = "model averaging", piecewiseSurvivalTime = 0,
    parameter_event_model = NULL,
    dropout_model = "weibull",
    parameter_dropout_model = NULL,
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nyears = 4, nreps = 500,
    showDropout = FALSE, showOngoing = FALSE) {

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

  if (!is.null(parameter_enroll_model)) {
    erify::check_class(parameter_enroll_model, "list")
    if (!is.null(df)) {
      if (tolower(parameter_enroll_model$model) != tolower(enroll_model)) {
        stop("Prior and likelihood must use the same enrollment model.")
      }

      if (tolower(parameter_enroll_model$model) == "piecewise poisson" &
          !all.equal(parameter_enroll_model$accrualTime, accrualTime)) {
        stop("Knots for piecewise Poisson must be the same for prior.")
      }

      if (tolower(parameter_enroll_model$model) == "b-spline") {
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

  if (!is.null(parameter_event_model)) {
    erify::check_class(parameter_event_model, "list")
    if (!is.null(df)) {
      if (tolower(parameter_event_model$model) != tolower(event_model)) {
        stop("Prior and likelihood must use the same event model.")
      }

      if (tolower(parameter_event_model$model) == "piecewise exponential" &
          !all.equal(parameter_event_model$piecewiseSurvivalTime,
          piecewiseSurvivalTime)) {
        stop("Knots for piecewise exponential must be the same for prior.")
      }

      if (tolower(parameter_event_model$model) == "model averaging") {
        stop("Model averaging event model cannot be used as prior.")
      }
    }
  }

  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-normal"))

  if (!is.null(parameter_dropout_model)) {
    erify::check_class(parameter_dropout_model, "list")
    if (!is.null(df)) {
      if (tolower(parameter_dropout_model$model) != tolower(dropout_model)) {
        stop("Prior and likelihood must use the same dropout model.")
      }
      if (tolower(parameter_dropout_model$model) == "none") {
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
    trialsdt = min(df$randdt)
    cutoffdt = df$cutoffdt[1]

    # summarize observed data
    observed <- summarizeObserved(df, to_predict, dropout_model)
  }

  # fit and predict enrollment
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) {
      erify::check_n(target_n - observed$n0,
                     supplement = "Enrollment target reached.")

      enroll_fit <- fitEnrollment(df = observed$adsl, enroll_model,
                                  nknots, accrualTime)

      # combine prior and likelihood to yield posterior
      if (!is.null(parameter_enroll_model)) {
        if (tolower(enroll_model) == "poisson") {
          enroll_fit$theta <-
            1/(1/enroll_fit$vtheta +
                 1/parameter_enroll_model$vtheta)*
            (1/enroll_fit$vtheta*enroll_fit$theta +
               1/parameter_enroll_model$vtheta*parameter_enroll_model$theta)
          enroll_fit$vtheta <-
            1/(1/enroll_fit$vtheta +
                 1/parameter_enroll_model$vtheta)
        } else {
          enroll_fit$theta <-
            solve(solve(enroll_fit$vtheta) +
                    solve(parameter_enroll_model$vtheta),
                  solve(enroll_fit$vtheta, enroll_fit$theta) +
                    solve(parameter_enroll_model$vtheta,
                          parameter_enroll_model$theta))
          enroll_fit$vtheta <- solve(solve(enroll_fit$vtheta) +
                                       solve(parameter_enroll_model$vtheta))
        }
      }

      # enrollment prediction at the analysis stage
      enroll_pred <- predictEnrollment(
        df = observed$adsl, target_n,
        enroll_fit = enroll_fit,
        lags, pilevel, nreps, showplot = FALSE)
    } else {
      # enrollment prediction at the design stage
      enroll_pred <- predictEnrollment(
        df = NULL, target_n,
        enroll_fit = parameter_enroll_model,
        lags, pilevel, nreps, showplot = FALSE)
    }
  }

  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) { # event prediction at analysis stage
      erify::check_n(target_d - observed$d0,
                     supplement = "Event target reached.")

      # convert prior by treatment to prior overall
      if (!is.null(parameter_event_model)) {
        if (tolower(parameter_event_model$model) == "exponential") {
          # match the overall mean
          w = parameter_event_model$prob
          lambda = exp(parameter_event_model$theta)

          event_prior <- list(
            model = parameter_event_model$model,
            theta = log(1/sum(w/lambda)),
            vtheta = parameter_event_model$vtheta[1,1])
        } else if (tolower(parameter_event_model$model) == "weibull") {
          # match the overall mean and variance
          k = parameter_event_model$ngroups
          w = parameter_event_model$prob
          shape = exp(parameter_event_model$theta[2*(1:k)-1])
          scale = exp(parameter_event_model$theta[2*(1:k)])

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
            model = parameter_event_model$model,
            theta = sol,
            vtheta = parameter_event_model$vtheta[1:2,1:2])
        } else if (tolower(parameter_event_model$model) == "log-normal") {
          # match the overall mean and variance
          k = parameter_event_model$ngroups
          w = parameter_event_model$prob
          meanlog = parameter_event_model$theta[2*(1:k)-1]
          sdlog = exp(parameter_event_model$theta[2*(1:k)])

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
            model = parameter_event_model$model,
            theta = sol,
            vtheta = parameter_event_model$vtheta[1:2,1:2])
        } else if (tolower(parameter_event_model$model) ==
                   "piecewise exponential") {
          # match within each interval
          k = parameter_event_model$ngroups
          w = parameter_event_model$prob
          lambda = exp(parameter_event_model$theta)
          npieces = length(lambda)/k

          event_prior <- list(
            model = parameter_event_model$model,
            theta = log(1/as.numeric(w %*% matrix(
              1/lambda, nrow=k, ncol=npieces, byrow=TRUE))),
            vtheta = parameter_event_model$vtheta[1:npieces,1:npieces],
            piecewiseSurvivalTime =
              parameter_event_model$piecewiseSurvivalTime)
        }
      }

      if (sum(observed$adtte$event == 1) > 0) {
        event_fit <- fitEvent(df = observed$adtte, event_model,
                              piecewiseSurvivalTime)
      } else {
        if (is.null(parameter_event_model)) {
          stop("Prior must be specified if there is no event observed.")
        }
        event_fit <- event_prior
        event_fit$vtheta = event_prior$vtheta*1e8
        event_fit$bic = NA
      }

      # combine prior and likelihood to yield posterior
      if (!is.null(parameter_event_model)) {
        if (tolower(event_model) == "exponential") {
          event_fit$theta <- 1/(1/event_fit$vtheta + 1/event_prior$vtheta)*
            (1/event_fit$vtheta*event_fit$theta +
               1/event_prior$vtheta*event_prior$theta)
          event_fit$vtheta <- 1/(1/event_fit$vtheta + 1/event_prior$vtheta)
        } else {
          event_fit$theta <-
            solve(solve(event_fit$vtheta) + solve(event_prior$vtheta),
                  solve(event_fit$vtheta, event_fit$theta) +
                    solve(event_prior$vtheta, event_prior$theta))
          event_fit$vtheta <-
            solve(solve(event_fit$vtheta) + solve(event_prior$vtheta))
        }
      }


      # whether to include dropout model
      if (tolower(dropout_model) != "none") {

        # convert prior by treatment to prior overall
        if (!is.null(parameter_dropout_model)) {
          if (tolower(parameter_dropout_model$model) == "exponential") {
            # match the overall mean
            w = parameter_dropout_model$prob
            lambda = exp(parameter_dropout_model$theta)

            dropout_prior <- list(
              model = parameter_dropout_model$model,
              theta = log(1/sum(w/lambda)),
              vtheta = parameter_dropout_model$vtheta[1,1])
          } else if (tolower(parameter_dropout_model$model) == "weibull") {
            # match the overall mean and variance
            k = parameter_dropout_model$ngroups
            w = parameter_dropout_model$prob
            shape = exp(parameter_dropout_model$theta[2*(1:k)-1])
            scale = exp(parameter_dropout_model$theta[2*(1:k)])

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
              model = parameter_dropout_model$model,
              theta = sol,
              vtheta = parameter_dropout_model$vtheta[1:2,1:2])
          } else if (tolower(parameter_dropout_model$model) == "log-normal") {
            # match the overall mean and variance
            k = parameter_dropout_model$ngroups
            w = parameter_dropout_model$prob
            meanlog = parameter_dropout_model$theta[2*(1:k)-1]
            sdlog = exp(parameter_dropout_model$theta[2*(1:k)])

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
              model = parameter_dropout_model$model,
              theta = sol,
              vtheta = parameter_dropout_model$vtheta[1:2,1:2])
          }
        }


        if (sum(observed$adtte$dropout == 1) > 0) {
          dropout_fit <- fitDropout(df = observed$adtte, dropout_model)
        } else {
          if (is.null(parameter_dropout_model)) {
            stop("Prior must be specified if there is no dropout observed.")
          }
          dropout_fit <- dropout_prior
          dropout_fit$vtheta <- dropout_prior$vtheta*1e8
          dropout_fit$bic = NA
        }

        # combine prior and likelihood to yield posterior
        if (!is.null(parameter_dropout_model)) {
          if (tolower(dropout_model) == "exponential") {
            dropout_fit$theta <-
              1/(1/dropout_fit$vtheta + 1/dropout_prior$vtheta)*
              (1/dropout_fit$vtheta*dropout_fit$theta +
                 1/dropout_prior$vtheta*dropout_prior$theta)
            dropout_fit$vtheta <-
              1/(1/dropout_fit$vtheta + 1/dropout_prior$vtheta)
          } else {
            dropout_fit$theta <-
              solve(solve(dropout_fit$vtheta) + solve(dropout_prior$vtheta),
                    solve(dropout_fit$vtheta, dropout_fit$theta) +
                      solve(dropout_prior$vtheta, dropout_prior$theta))
            dropout_fit$vtheta <-
              solve(solve(dropout_fit$vtheta) + solve(dropout_prior$vtheta))
          }
        }


        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit,
            dropout_fit = dropout_fit,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showDropout, showOngoing, showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = NULL,
            event_fit = event_fit,
            dropout_fit = dropout_fit,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showDropout, showOngoing, showplot = FALSE)
        }
      } else {
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit,
            dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showDropout, showOngoing, showplot = FALSE)
        } else {
          event_pred <- predictEvent(
            df = observed$adtte, target_d,
            newSubjects = NULL,
            event_fit = event_fit,
            dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showDropout, showOngoing, showplot = FALSE)
        }
      }
    } else { # event prediction at design stage
      if (!is.null(parameter_dropout_model)) {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = parameter_event_model,
          dropout_fit = parameter_dropout_model,
          fixedFollowup, followupTime, pilevel, nyears, nreps,
          showDropout, showOngoing, showplot = FALSE)
      } else {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = parameter_event_model,
          dropout_fit = NULL,
          fixedFollowup, followupTime, pilevel, nyears, nreps,
          showDropout, showOngoing, showplot = FALSE)
      }
    }
  }


  # output results
  if (!is.null(df)) { # analysis stage prediction
    if (tolower(to_predict) == "enrollment only") {
      dfs <- enroll_pred$plotEnrollment

      # separate data into observed and predicted
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      # plot the enrollment data with month as x-axis label
      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$date,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_step(data=dfa, ggplot2::aes(x=.data$date, y=.data$n),
                           color="black") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$date, y=.data$n),
                           color="blue") +
        ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
        ggplot2::scale_x_date(name = NULL,
                              labels = scales::date_format("%b"),
                              breaks = scales::breaks_width(bw),
                              minor_breaks = NULL,
                              expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Subjects", title = "Predicted subjects") +
        ggplot2::theme_bw()

      g2 <- flabel(dfs, trialsdt)

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      print(p1)

      list(observed = observed, enroll_fit = enroll_fit,
           enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "event only") {
      dfs <- event_pred$plotEvent

      # separate data into observed and predicted
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      # plot the enrollment and time to event data with month as x-axis label
      # generate the plot
      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$date,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_step(data=dfa, ggplot2::aes(x=.data$date, y=.data$n),
                           color="black") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$date, y=.data$n),
                           color="blue") +
        ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_date(name = NULL,
                              labels = scales::date_format("%b"),
                              breaks = scales::breaks_width(bw),
                              minor_breaks = NULL,
                              expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Events", title = "Predicted events") +
        ggplot2::theme_bw()

      g2 <- flabel(dfs, trialsdt)

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      print(p1)



      # whether to show dropout subjects
      if (showDropout) {
        dft <- event_pred$plotDropout

        # separate data into observed and predicted
        dfc <- dft %>% dplyr::filter(is.na(.data$lower))
        dfd <- dft %>% dplyr::filter(!is.na(.data$lower))

        n_months = lubridate::interval(min(dft$date),
                                       max(dft$date)) %/% months(1)
        bw = fbw(n_months)

        # generate the plot
        g3 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dfd, ggplot2::aes(x=.data$date,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_step(data=dfc, ggplot2::aes(x=.data$date, y=.data$n),
                             color="black") +
          ggplot2::geom_line(data=dfd, ggplot2::aes(x=.data$date, y=.data$n),
                             color="blue") +
          ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
          ggplot2::scale_x_date(name = NULL,
                                labels = scales::date_format("%b"),
                                breaks = scales::breaks_width(bw),
                                minor_breaks = NULL,
                                expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Dropouts", title = "Predicted dropouts") +
          ggplot2::theme_bw()

        g4 <- flabel(dft, trialsdt)

        # stack them together
        p2 <- g3 + g4 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
        print(p2)
      }



      if (showOngoing) {
        dfu <- event_pred$plotOngoing

        # separate data into observed and predicted
        dfe <- dfu %>% dplyr::filter(is.na(.data$lower))
        dff <- dfu %>% dplyr::filter(!is.na(.data$lower))

        n_months = lubridate::interval(min(dfu$date),
                                       max(dfu$date)) %/% months(1)
        bw = fbw(n_months)

        # generate the plot
        g5 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dff, ggplot2::aes(x=.data$date,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_step(data=dfe, ggplot2::aes(x=.data$date, y=.data$n),
                             color="black") +
          ggplot2::geom_line(data=dff, ggplot2::aes(x=.data$date, y=.data$n),
                             color="blue") +
          ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
          ggplot2::scale_x_date(name = NULL,
                                labels = scales::date_format("%b"),
                                breaks = scales::breaks_width(bw),
                                minor_breaks = NULL,
                                expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Ongoing subjects",
                        title = "Predicted ongoing subjects") +
          ggplot2::theme_bw()

        g6 <- flabel(dfu, trialsdt)

        # stack them together
        p3 <- g5 + g6 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
        print(p3)
      }


      if (tolower(dropout_model) != "none") {
        list(observed = observed, event_fit = event_fit,
             dropout_fit = dropout_fit, event_pred = event_pred)
      } else {
        list(observed = observed, event_fit = event_fit,
             event_pred = event_pred)
      }
    } else if (tolower(to_predict) == "enrollment and event") {
      df1 <- enroll_pred$plotEnrollment %>%
        dplyr::mutate(parameter = "subjects")
      df1last <- df1 %>% dplyr::slice(dplyr::n())

      df2 <- event_pred$plotEvent %>%
        dplyr::mutate(parameter = "events")
      df2last <- df2 %>% dplyr::slice(dplyr::n())

      # extend enrollment prediction time to event prediction time
      df12 <- df2last %>%
        dplyr::mutate(n = df1last$n,
                        lower = df1last$lower,
                        upper = df1last$upper,
                        parameter = df1last$parameter)

      dfs <- df1 %>%
        dplyr::bind_rows(df12) %>%
        dplyr::bind_rows(df2)

      # separate data into observed and predicted
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$date,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper,
                                                    group=.data$parameter),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_step(data=dfa, ggplot2::aes(x=.data$date, y=.data$n,
                                                  group=.data$parameter),
                           color="black") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$date, y=.data$n,
                                                  group=.data$parameter),
                           color="blue") +
        ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_date(name = NULL,
                              labels = scales::date_format("%b"),
                              breaks = scales::breaks_width(bw),
                              minor_breaks = NULL,
                              expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Subjects / Events",
                      title = "Predicted subjects and events") +
        ggplot2::theme_bw()

      g2 <- flabel(dfs, trialsdt)

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      print(p1)


      # whether to show dropout subjects
      if (showDropout) {
        dft <- event_pred$plotDropout

        # separate data into observed and predicted
        dfc <- dft %>% dplyr::filter(is.na(.data$lower))
        dfd <- dft %>% dplyr::filter(!is.na(.data$lower))

        n_months = lubridate::interval(min(dft$date),
                                       max(dft$date)) %/% months(1)
        bw = fbw(n_months)

        # generate the plot
        g3 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dfd, ggplot2::aes(x=.data$date,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_step(data=dfc, ggplot2::aes(x=.data$date, y=.data$n),
                             color="black") +
          ggplot2::geom_line(data=dfd, ggplot2::aes(x=.data$date, y=.data$n),
                             color="blue") +
          ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
          ggplot2::scale_x_date(name = NULL,
                                labels = scales::date_format("%b"),
                                breaks = scales::breaks_width(bw),
                                minor_breaks = NULL,
                                expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Dropouts", title = "Predicted dropouts") +
          ggplot2::theme_bw()

        g4 <- flabel(dft, trialsdt)

        # stack them together
        p2 <- g3 + g4 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
        print(p2)
      }


      if (showOngoing) {
        dfu <- event_pred$plotOngoing

        # separate data into observed and predicted
        dfe <- dfu %>% dplyr::filter(is.na(.data$lower))
        dff <- dfu %>% dplyr::filter(!is.na(.data$lower))

        n_months = lubridate::interval(min(dfu$date),
                                       max(dfu$date)) %/% months(1)
        bw = fbw(n_months)

        # generate the plot
        g5 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dff, ggplot2::aes(x=.data$date,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_step(data=dfe, ggplot2::aes(x=.data$date, y=.data$n),
                             color="black") +
          ggplot2::geom_line(data=dff, ggplot2::aes(x=.data$date, y=.data$n),
                             color="blue") +
          ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
          ggplot2::scale_x_date(name = NULL,
                                labels = scales::date_format("%b"),
                                breaks = scales::breaks_width(bw),
                                minor_breaks = NULL,
                                expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Ongoing subjects",
                        title = "Predicted ongoing subjects") +
          ggplot2::theme_bw()

        g6 <- flabel(dfu, trialsdt)

        # stack them together
        p3 <- g5 + g6 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
        print(p3)
      }

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
      dfb <- enroll_pred$plotEnrollment

      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$t,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$t, y=.data$n),
                           color="blue") +
        ggplot2::scale_x_continuous(name = "Days since randomization",
                                    expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Subjects",
                      title = "Predicted subjects") +
        ggplot2::theme_bw()

      print(g1)

      list(enroll_fit = parameter_enroll_model, enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "event only") {
      dfb <- event_pred$plotEvent

      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$t,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$t, y=.data$n),
                           color="blue") +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_continuous(name = "Days since randomization",
                                    expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Events", title = "Predicted events") +
        ggplot2::theme_bw()

      print(g1)


      if (showDropout) {
        dfd <- event_pred$plotDropout

        g3 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dfd, ggplot2::aes(x=.data$t,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_line(data=dfd, ggplot2::aes(x=.data$t, y=.data$n),
                             color="blue") +
          ggplot2::scale_x_continuous(name = "Days since randomization",
                                      expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Dropouts", title = "Predicted dropouts") +
          ggplot2::theme_bw()

        print(g3)
      }


      if (showOngoing) {
        dff <- event_pred$plotOngoing

        g5 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dff, ggplot2::aes(x=.data$t,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_line(data=dff, ggplot2::aes(x=.data$t, y=.data$n),
                             color="blue") +
          ggplot2::scale_x_continuous(name = "Days since randomization",
                                      expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Ongoing subjects",
                        title = "Predicted ongoing subjects") +
          ggplot2::theme_bw()

        print(g5)
      }

      if (!is.null(parameter_dropout_model)) {
        list(event_fit = parameter_event_model,
             dropout_fit = parameter_dropout_model, event_pred = event_pred)
      } else {
        list(event_fit = parameter_event_model, event_pred = event_pred)
      }
    } else if (tolower(to_predict) == "enrollment and event") {
      df1 <- enroll_pred$plotEnrollment %>%
        dplyr::mutate(parameter = "subjects")
      df1last <- df1 %>% dplyr::slice(dplyr::n())

      df2 <- event_pred$plotEvent %>%
        dplyr::mutate(parameter = "events")
      df2last <- df2 %>% dplyr::slice(dplyr::n())

      # extend enrollment prediction time to event prediction time
      df12 <- df2last %>%
        dplyr::mutate(n = df1last$n,
                      lower = df1last$lower,
                      upper = df1last$upper,
                      parameter = df1last$parameter)

      dfb <- df1 %>%
        dplyr::bind_rows(df12) %>%
        dplyr::bind_rows(df2)

      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$t,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper,
                                                    group=.data$parameter),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$t, y=.data$n,
                                                  group=.data$parameter),
                           color="blue") +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_continuous(name = "Days since randomization",
                                    expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Subjects / Events",
                      title = "Predicted subjects and events") +
        ggplot2::theme_bw()

      print(g1)


      if (showDropout) {
        dfd <- event_pred$plotDropout

        g3 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dfd, ggplot2::aes(x=.data$t,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_line(data=dfd, ggplot2::aes(x=.data$t, y=.data$n),
                             color="blue") +
          ggplot2::scale_x_continuous(name = "Days since randomization",
                                      expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Dropouts", title = "Predicted dropouts") +
          ggplot2::theme_bw()

        print(g3)
      }


      if (showOngoing) {
        dff <- event_pred$plotOngoing

        g5 <- ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=dff, ggplot2::aes(x=.data$t,
                                                      ymin=.data$lower,
                                                      ymax=.data$upper),
                               alpha=0.5, fill="lightblue") +
          ggplot2::geom_line(data=dff, ggplot2::aes(x=.data$t, y=.data$n),
                             color="blue") +
          ggplot2::scale_x_continuous(name = "Days since randomization",
                                      expand = c(0.01, 0.01)) +
          ggplot2::labs(y = "Ongoing subjects",
                        title = "Predicted ongoing subjects") +
          ggplot2::theme_bw()

        print(g5)
      }

      if (!is.null(parameter_dropout_model)) {
        list(enroll_fit = parameter_enroll_model, enroll_pred = enroll_pred,
             event_fit = parameter_event_model,
             dropout_fit = parameter_dropout_model, event_pred = event_pred)
      } else {
        list(enroll_fit = parameter_enroll_model, enroll_pred = enroll_pred,
             event_fit = parameter_event_model, event_pred = event_pred)
      }
    }
  }

}


