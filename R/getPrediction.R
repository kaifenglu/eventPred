#' @title Enrollment and event prediction
#' @description Performs enrollment and event prediction by utilizing
#'   observed data and specified enrollment and event models.
#'
#' @param df The subject-level enrollment and event data, including
#'   \code{trialsdt}, \code{usubjid}, \code{randdt}, and \code{cutoffdt} for
#'   enrollment prediction, and, additionally, \code{time}, \code{event},
#'   and \code{dropout} for event prediction. The data should also include
#'   \code{treatment} coded as 1, 2, and so on, and
#'   \code{treatment_description} for enrollment and
#'   event prediction by treatment. By default, it is set to
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
#' @param enroll_prior The prior of enrollment model parameters.
#' @param event_model The event model used to analyze the event data
#'   which can be set to one of the following options:
#'   "exponential", "Weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", or "model averaging". The model averaging
#'   uses the \code{exp(-bic/2)} weighting and combines Weibull and
#'   log-normal models. By default, it is set to "model
#'   averaging".
#' @param piecewiseSurvivalTime A vector that specifies the time
#'   intervals for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param k The number of inner knots of the spline event model of
#'   Royston and Parmar (2002). The default
#'   \code{k=0} gives a Weibull, log-logistic or log-normal model,
#'   if \code{scale} is "hazard", "odds", or "normal", respectively.
#'   The knots are chosen as equally-spaced quantiles of the log
#'   uncensored survival times. The boundary knots are chosen as the
#'   minimum and maximum log uncensored survival times.
#' @param scale If "hazard", the log cumulative hazard is modeled
#'   as a spline function. If "odds", the log cumulative odds is
#'   modeled as a spline function. If "normal", -qnorm(S(t)) is
#'   modeled as a spline function.
#' @param event_prior The prior of event model parameters.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options: "exponential",
#'   "Weibull", "log-logistic", "log-normal", or "piecewise exponential".
#'   By default, it is set to "exponential".
#' @param piecewiseDropoutTime A vector that specifies the time
#'   intervals for the piecewise exponential dropout distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param dropout_prior The prior of dropout model parameters.
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
#' @param by_treatment A Boolean variable to control whether or not to
#'   predict by treatment group. By default, it is set to \code{FALSE}.
#' @param ngroups The number of treatment groups for enrollment prediction
#'   at the design stage. By default, it is set to 1.
#'   It is replaced with the actual number of
#'   treatment groups in the observed data if \code{df} is not \code{NULL}.
#' @param alloc The treatment allocation in a randomization block.
#'   By default, it is set to \code{NULL}, which yields equal allocation
#'   among the treatment groups.
#' @param treatment_label The treatment labels for treatments in a
#'   randomization block for design stage prediction.
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
#' The \code{enroll_prior} variable should be a list that
#' includes \code{model} to specify the enrollment model
#' (poisson, time-decay, or piecewise poisson),
#' \code{theta} and \code{vtheta} to indicate the parameter
#' values and the covariance matrix. One can use a very small
#' value of \code{vtheta} to fix the parameter values.
#' For the piecewise Poisson enrollment model, the list
#' should also include \code{accrualTime}. It should be noted
#' that the B-spline model is not appropriate for use as prior.
#'
#' The \code{event_prior} variable should be a list with one element
#' per treatment. For each treatment, the element should include \code{w}
#' to specify the weight of the treatment in a randomization block,
#' \code{model} to specify the event model
#' (exponential, weibull, log-logistic, log-normal,
#' or piecewise exponential), \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix.
#' For the piecewise exponential event model, the list
#' should also include \code{piecewiseSurvivalTime} to indicate
#' the location of knots. It should be noted that the model averaging
#' and spline options are not appropriate for use as prior.
#'
#' The \code{dropout_prior} should be a list with one element
#' per treatment. For each treatment, the element should include \code{w}
#' to specify the weight of the treatment in a randomization block,
#' \code{model} to specify the dropout model
#' (exponential, weibull, log-logistic, log-normal,
#' or piecewise exponential), \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix.
#' For the piecewise exponential dropout model, the list
#' should also include \code{piecewiseDropoutTime} to indicate
#' the location of knots.
#'
#' For analysis-stage enrollment and event prediction, the
#' \code{enroll_prior}, \code{event_prior}, and
#' \code{dropout_prior} are either set to \code{NULL} to
#' use the observed data only, or specify the prior distribution
#' of model parameters to be combined with observed data likelihood
#' for enhanced modeling flexibility.
#'
#' @return A list that includes the fits of observed data models,
#' as well as simulated enrollment data for new subjects and
#' simulated event data for ongoing and new subjects.
#'
#' @examples
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
    enroll_prior = NULL,
    event_model = "model averaging", piecewiseSurvivalTime = 0,
    k = 0, scale = "hazard",
    event_prior = NULL,
    dropout_model = "exponential", piecewiseDropoutTime = 0,
    dropout_prior = NULL,
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nyears = 4, nreps = 500,
    showEnrollment = TRUE, showEvent = TRUE,
    showDropout = FALSE, showOngoing = FALSE,
    showsummary = TRUE, showplot = TRUE,
    by_treatment = FALSE, ngroups = 1, alloc = NULL,
    treatment_label = NULL) {

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


  # check by_treatment, ngroups, and alloc
  erify::check_bool(by_treatment)

  if (is.null(df)) by_treatment = TRUE

  if (by_treatment) {
    if (!is.null(df)) {
      ngroups = length(table(df$treatment))
    }

    if (is.null(alloc)) {
      alloc = rep(1, ngroups)
    } else {
      if (length(alloc) != ngroups) {
        stop("length of alloc must be equal to the number of treatments.")
      }

      if (any(alloc <= 0 | alloc != round(alloc))) {
        stop("elements of alloc must be positive integers.")
      }
    }
  } else {
    ngroups = 1
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }


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

  # check enrollment model prior
  if (!is.null(enroll_prior)) {
    erify::check_class(enroll_prior, "list")
    erify::check_content(tolower(enroll_prior$model), c(
      "poisson", "time-decay", "piecewise poisson"))

    model = tolower(enroll_prior$model)
    p = length(enroll_prior$theta)
    vtheta = enroll_prior$vtheta

    if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                   ncol(vtheta) != p)) ||
        (p == 1 && length(c(vtheta)) != 1)) {
      stop(paste("Dimensions of vtheta must be compatible with the length",
                 "of theta in enroll_prior"))
    }

    if ((model == "poisson" && p != 1) ||
        (model == "time-decay" && p != 2) ||
        (model == "piecewise poisson" &&
         p != length(enroll_prior$accrualTime))) {
      stop(paste("Length of theta must be compatible with model",
                 "in enroll_prior"))
    }

    if (model == "piecewise poisson") {
      if (enroll_prior$accrualTime[1] != 0) {
        stop("accrualTime must start with 0 in enroll_prior")
      }
      if (length(enroll_prior$accrualTime) > 1 &&
          any(diff(enroll_prior$accrualTime) <= 0)) {
        stop("accrualTime should be increasing in enroll_prior")
      }
    }

    if (!is.null(df)) {
      if (tolower(enroll_prior$model) != tolower(enroll_model)) {
        stop("Prior and likelihood must use the same enrollment model.")
      }

      if (tolower(enroll_prior$model) == "piecewise poisson" &&
          (length(enroll_prior$accrualTime) < length(accrualTime) ||
           !all.equal(enroll_prior$accrualTime[1:length(accrualTime)],
                      accrualTime))) {
        stop(paste("accrualTime of piecewise Poisson must be a subset of",
                   "that in enroll_prior."))
      }
    }
  }


  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "model averaging", "spline"))

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }
  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  erify::check_n(k, zero = TRUE)
  erify::check_content(tolower(scale), c("hazard", "odds", "normal"))


  # check event model prior
  if (!is.null(event_prior)) {
    erify::check_class(event_prior, "list")

    if (by_treatment) {
      if (length(event_prior) != ngroups) {
        stop("event_prior must be a list with one element per treatment.")
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(event_prior)) {
      event_prior2 <- list()
      event_prior2[[1]] <- event_prior
    } else {
      event_prior2 <- event_prior
    }

    for (j in 1:length(event_prior2)) {
      erify::check_content(tolower(event_prior2[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential"))

      model = tolower(event_prior2[[j]]$model)
      p = length(event_prior2[[j]]$theta)
      vtheta = event_prior2[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in event_prior"))
      }

      if ((model == "exponential" && p != 1) ||
          (model == "weibull" && p != 2) ||
          (model == "log-logistic" && p != 2) ||
          (model == "log-normal" && p != 2) ||
          (model == "piecewise exponential" &&
           p != length(event_prior2[[j]]$piecewiseSurvivalTime))) {
        stop(paste("Length of theta must be compatible with model",
                   "in event_prior"))
      }

      if (model == "piecewise exponential") {
        if (event_prior2[[j]]$piecewiseSurvivalTime[1] != 0) {
          stop(paste("piecewiseSurvivalTime must start with 0",
                     "in event_prior"))
        }
        if (length(event_prior2[[j]]$piecewiseSurvivalTime) > 1 &&
            any(diff(event_prior2[[j]]$piecewiseSurvivalTime) <= 0)) {
          stop(paste("piecewiseSurvivalTime should be increasing",
                     "in event_prior"))
        }
      }


      if (!is.null(df)) {
        if (tolower(event_prior2[[j]]$model) != tolower(event_model)) {
          stop("Prior and likelihood must use the same event model.")
        }

        if (tolower(event_prior2[[j]]$model) == "piecewise exponential" &&
            (length(event_prior2[[j]]$piecewiseSurvivalTime) <
             length(piecewiseSurvivalTime) ||
             !all.equal(event_prior2[[j]]$piecewiseSurvivalTime[
               1:length(piecewiseSurvivalTime)], piecewiseSurvivalTime))) {
          stop(paste("piecewiseSurvivalTime of piecewise exponential model",
                     "must be a subset of that in event_prior."))
        }
      }
    }
  }

  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential"))

  if (piecewiseDropoutTime[1] != 0) {
    stop("piecewiseDropoutTime must start with 0")
  }
  if (length(piecewiseDropoutTime) > 1 &&
      any(diff(piecewiseDropoutTime) <= 0)) {
    stop("piecewiseDropoutTime should be increasing")
  }

  # check dropout model prior
  if (!is.null(dropout_prior)) {
    erify::check_class(dropout_prior, "list")

    if (by_treatment) {
      if (length(dropout_prior) != ngroups) {
        stop("dropout_prior must be a list with one element per treatment.")
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(dropout_prior)) {
      dropout_prior2 <- list()
      dropout_prior2[[1]] <- dropout_prior
    } else {
      dropout_prior2 <- dropout_prior
    }

    for (j in 1:length(dropout_prior2)) {
      erify::check_content(tolower(dropout_prior2[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential"))

      model = tolower(dropout_prior2[[j]]$model)
      p = length(dropout_prior2[[j]]$theta)
      vtheta = dropout_prior2[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in dropout_prior"))
      }

      if ((model == "exponential" && p != 1) ||
          (model == "weibull" && p != 2) ||
          (model == "log-logistic" && p != 2) ||
          (model == "log-normal" && p != 2) ||
          (model == "piecewise exponential" &&
           p != length(dropout_prior2[[j]]$piecewiseDropoutTime))) {
        stop(paste("Length of theta must be compatible with model",
                   "in event_prior"))
      }

      if (model == "piecewise exponential") {
        if (dropout_prior2[[j]]$piecewiseDropoutTime[1] != 0) {
          stop(paste("piecewiseDropoutTime must start with 0",
                     "in dropout_prior"))
        }
        if (length(dropout_prior2[[j]]$piecewiseDropoutTime) > 1 &&
            any(diff(dropout_prior2[[j]]$piecewiseDropoutTime) <= 0)) {
          stop(paste("piecewiseDropoutTime should be increasing",
                     "in dropout_prior"))
        }
      }


      if (!is.null(df)) {
        if (tolower(dropout_prior2[[j]]$model) != tolower(dropout_model)) {
          stop("Prior and likelihood must use the same dropout model.")
        }

        if (tolower(dropout_prior2[[j]]$model) == "piecewise exponential" &&
            (length(dropout_prior2[[j]]$piecewiseDropoutTime) <
             length(piecewiseDropoutTime) ||
             !all.equal(dropout_prior2[[j]]$piecewiseDropoutTime[
               1:length(piecewiseDropoutTime)], piecewiseDropoutTime))) {
          stop(paste("piecewiseDropoutTime of piecewise exponential model",
                     "must be a subset of that in dropout_prior."))
        }
      }


      if (!is.null(event_prior) && "w" %in% names(event_prior2[[j]])) {
        if (event_prior2[[j]]$w != dropout_prior2[[j]]$w) {
          stop("w must be equal between event prior and dropout prior.")
        }
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
    observed <- summarizeObserved(df, to_predict, showplot, by_treatment)
  }


  # fit and predict enrollment
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) {
      erify::check_n(target_n - observed$n0,
                     supplement = "Enrollment target reached.")

      enroll_fit <- fitEnrollment(df = df, enroll_model,
                                  nknots, accrualTime, showplot)
      enroll_fit1 <- enroll_fit$enroll_fit

      # combine prior and likelihood to yield posterior
      if (!is.null(enroll_prior)) {
        # pad additional pieces with prior parameter values
        if (tolower(enroll_model) == "piecewise poisson" &&
            length(enroll_prior$accrualTime) > length(accrualTime)) {
          i = (length(accrualTime) + 1):length(enroll_prior$accrualTime)
          enroll_fit1$theta = c(enroll_fit1$theta, enroll_prior$theta[i])
          enroll_fit1$vtheta = as.matrix(Matrix::bdiag(
            enroll_fit1$vtheta, enroll_prior$vtheta[i,i]*1e8))
          enroll_fit1$accrualTime = enroll_prior$accrualTime
        }

        enroll_fit1$theta <-
          solve(solve(enroll_fit1$vtheta) + solve(enroll_prior$vtheta),
                solve(enroll_fit1$vtheta, enroll_fit1$theta) +
                  solve(enroll_prior$vtheta, enroll_prior$theta))
        enroll_fit1$vtheta <-
          solve(solve(enroll_fit1$vtheta) + solve(enroll_prior$vtheta))
      }

      # enrollment prediction at the analysis stage
      enroll_pred <- predictEnrollment(
        df = df, target_n, enroll_fit = enroll_fit1,
        lags, pilevel, nyears, nreps, showsummary, showplot = FALSE,
        by_treatment, ngroups, alloc, treatment_label)
    } else {
      # enrollment prediction at the design stage
      enroll_pred <- predictEnrollment(
        df = NULL, target_n, enroll_fit = enroll_prior,
        lags, pilevel, nyears, nreps, showsummary, showplot = FALSE,
        by_treatment, ngroups, alloc, treatment_label)
    }
  }


  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) { # event prediction at analysis stage
      erify::check_n(target_d - observed$d0,
                     supplement = "Event target reached.")

      if (by_treatment) {
        sum_by_trt <- df %>%
          dplyr::group_by(.data$treatment) %>%
          dplyr::summarise(n0 = dplyr::n(),
                           d0 = sum(.data$event),
                           c0 = sum(.data$dropout),
                           r0 = sum(!(.data$event | .data$dropout)))
      }

      # convert prior by treatment to prior overall
      if (!is.null(event_prior) && !by_treatment &&
          !("model" %in% names(event_prior))) {

        m = length(event_prior)
        w = sapply(event_prior, function(sub_list) sub_list$w)
        w = w/sum(w)

        # check model consistency across treatments
        model = tolower(event_prior[[1]]$model)
        if (m > 1) {
          for (j in 2:m) {
            if (tolower(event_prior[[j]]$model) != model) {
              stop("Prior event model must be equal across treatments.")
            }
          }

          if (model == "piecewise exponential") {
            for (j in 2:m) {
              if (!all.equal(event_prior[[j]]$piecewiseSurvivalTime,
                             event_prior[[1]]$piecewiseSurvivalTime)) {
                stop(paste("piecewiseSurvivalTime must be equal",
                           "across treatments."))
              }
            }
          }
        }


        if (model == "exponential") {
          # match the overall mean
          theta = sapply(event_prior, function(sub_list) sub_list$theta)
          lambda = exp(theta)
          lambda1 = 1/sum(w/lambda)  # hazard rate for pooled
          theta1 = log(lambda1)

          # use delta-method to obtain the variance
          vtheta1 = 0
          for (i in 1:m) {
            vtheta1 = vtheta1 + (w[i]/lambda[i])^2*event_prior[[i]]$vtheta
          }
          vtheta1 = vtheta1*lambda1^2

          event_prior1 <- list(
            model = model, theta = theta1, vtheta = vtheta1)
        } else if (model == "weibull") {
          # match the overall mean and variance

          # mean and variance of weibull as a function of theta
          fmweibull <- function(theta) {
            shape = exp(-theta[2])
            scale = exp(theta[1])
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
          theta = lapply(event_prior, function(sub_list) sub_list$theta)
          m1 = lapply(theta, fmweibull)
          m1mean = sapply(m1, function(sub_list) sub_list$mean)
          m1var = sapply(m1, function(sub_list) sub_list$var)

          # mean and variance for pooled
          m2 = list(mean = sum(w*m1mean),
                    var = sum(w*m1var) +
                      sum(w*m1mean^2) - (sum(w*m1mean))^2)

          # solve for theta given the mean and variance for pooled
          theta2s = sapply(theta, function(sub_list) sub_list[2])

          theta12 = uniroot(function(x)
            lgamma(1+2/exp(-x)) - 2*lgamma(1+1/exp(-x)) -
              log(m2$var/m2$mean^2 + 1),
            c(min(theta2s) - 1, max(theta2s) + 1), extendInt = "yes")$root

          theta11 = log(m2$mean) - lgamma(1+1/exp(-theta12))
          theta1 = c(theta11, theta12)

          # gradient of theta with respect to mean and variance for pooled
          ig = solve(gmweibull(theta1))

          # variance of theta for pooled
          vtheta1 = 0
          for (i in 1:m) {
            gi = gmweibull(event_prior[[i]]$theta)
            vm1i = gi * event_prior[[i]]$vtheta * t(gi)
            li = w[i]*matrix(c(1, 2*(m1[[i]]$mean - m2$mean), 0, 1), ncol=2)
            vtheta1 = vtheta1 + li %*% vm1i %*% t(li)
          }
          vtheta1 = ig %*% vtheta1 %*% t(ig)

          event_prior1 <- list(
            model = model, theta = theta1, vtheta = vtheta1)
        } else if (model == "log-logistic") {
          # since the mean and variance of log-logistic distribution
          # may not exist, match the cdf at the weighted average of
          # treatment-specific 97.5% percentiles and medians

          fllogis <- function(theta) {
            k = length(theta)/2
            mus = theta[seq(1,2*k-1,2)]
            sigmas = exp(theta[seq(2,2*k,2)])
            t1 = sum(w*exp(mus + qlogis(0.975)*sigmas))
            t2 = sum(w*exp(mus))
            a1 = log(1/sum(w*plogis(-(log(t1) - mus)/sigmas)) - 1)
            a2 = log(1/sum(w*plogis(-(log(t2) - mus)/sigmas)) - 1)
            gamma = (a1 - a2)/(log(t1) - log(t2))
            mu = log(t1) - 1/gamma*a1
            c(mu, -log(gamma))
          }

          # gradient vector
          gllogis <- function(theta) {
            g1 = numDeriv::grad(function(theta) fllogis(theta)[1], theta)
            g2 = numDeriv::grad(function(theta) fllogis(theta)[2], theta)
            matrix(c(g1, g2), nrow=2, byrow=TRUE)
          }

          # concatenating treatment-specific model parameters
          theta = lapply(event_prior, function(sub_list) sub_list$theta)

          # parameter and variance for the overall population
          theta1 = fllogis(theta)
          g = gllogis(theta)
          vtheta1 = 0
          for (i in 1:m) {
            gi = g[,(2*i-1):(2*i)]
            vtheta1 = vtheta1 + gi * event_prior[[i]]$vtheta * t(gi)
          }

          event_prior1 <- list(
            model = model, theta = theta1, vtheta = vtheta1)
        } else if (model == "log-normal") {
          # match the overall mean and variance

          # mean and variance of log-normal as a function of theta
          fmlnorm <- function(theta) {
            meanlog = theta[1]
            sdlog = exp(theta[2])
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
          theta = lapply(event_prior, function(sub_list) sub_list$theta)
          m1 = lapply(theta, fmlnorm)
          m1mean = sapply(m1, function(sub_list) sub_list$mean)
          m1var = sapply(m1, function(sub_list) sub_list$var)

          # mean and variance for pooled
          m2 = list(mean = sum(w*m1mean),
                    var = sum(w*m1var) +
                      sum(w*m1mean^2) - (sum(w*m1mean))^2)

          # solve for theta given the mean and variance for pooled
          theta12 = 0.5*log(log(m2$var/m2$mean^2 + 1))
          theta11 = log(m2$mean) - 0.5*exp(2*theta12)
          theta1 = c(theta11, theta12)

          # gradient of theta with respect to mean and variance for pooled
          ig = solve(gmlnorm(theta1))

          # variance of theta for pooled
          vtheta1 = 0
          for (i in 1:m) {
            gi = gmlnorm(event_prior[[i]]$theta)
            vm1i = gi * event_prior[[i]]$vtheta * t(gi)
            li = w[i]*matrix(c(1, 2*(m1[[i]]$mean - m2$mean), 0, 1), ncol=2)
            vtheta1 = vtheta1 + li %*% vm1i %*% t(li)
          }
          vtheta1 = ig %*% vtheta1 %*% t(ig)

          event_prior1 <- list(
            model = model, theta = theta1, vtheta = vtheta1)
        } else if (model == "piecewise exponential") {
          # match 1/lambda within each interval
          theta = sapply(event_prior, function(sub_list) sub_list$theta)
          lambda = exp(theta)
          npieces = length(event_prior[[1]]$piecewiseSurvivalTime)

          # construct theta and vtheta piece by piece
          theta1 = rep(NA, npieces)
          vtheta1 = 0*diag(npieces)
          for (j in 1:npieces) {
            lambdaj = lambda[j,]
            lambda1j = 1/sum(w/lambdaj)
            theta1[j] = log(lambda1j)
            for (i in 1:ngroups) {
              vtheta1[j,j] = vtheta1[j,j] + (w[i]/lambdaj[i])^2 *
                event_prior[[i]]$vtheta[j,j]
            }
            vtheta1[j,j] = vtheta1[j,j]*lambda1j^2
          }

          event_prior1 <- list(
            model = model, theta = theta1, vtheta = vtheta1,
            piecewiseSurvivalTime = event_prior[[1]]$piecewiseSurvivalTime)
        }
      } else if (!is.null(event_prior)) {
        event_prior1 <- event_prior
      }

      # fit the event model
      if ((!by_treatment && observed$d0 > 0) ||
          (by_treatment && all(sum_by_trt$d0 > 0))) {
        event_fit <- fitEvent(df = df, event_model,
                              piecewiseSurvivalTime, k, scale, showplot,
                              by_treatment)
        event_fit1 <- event_fit$event_fit
      } else {
        if (is.null(event_prior)) {
          stop("Prior must be specified if there is no event observed.")
        }

        event_fit <- list()
        event_fit1 <- event_prior1

        # inflate the variance
        if (!by_treatment) {
          event_fit1$vtheta <- event_prior1$vtheta*1e8
          event_fit1$bic <- NA
        } else {
          for (i in 1:ngroups) {
            event_fit1[[i]]$vtheta <- event_prior1[[i]]$vtheta*1e8
            event_fit1[[i]]$bic <- NA
          }
        }
      }

      # combine prior and likelihood to yield posterior
      if (!is.null(event_prior)) {
        if (!by_treatment) {
          # pad additional pieces with prior parameter values
          if (tolower(event_model) == "piecewise exponential" &&
              length(event_prior1$piecewiseSurvivalTime) >
              length(piecewiseSurvivalTime)) {
            i = (length(piecewiseSurvivalTime) + 1):
              length(event_prior1$piecewiseSurvivalTime)
            event_fit1$theta = c(event_fit1$theta, event_prior1$theta[i])
            event_fit1$vtheta = as.matrix(Matrix::bdiag(
              event_fit1$vtheta, event_prior1$vtheta[i,i]*1e8))
            event_fit1$piecewiseSurvivalTime =
              event_prior1$piecewiseSurvivalTime
          }

          event_fit1$theta <-
            solve(solve(event_fit1$vtheta) + solve(event_prior1$vtheta),
                  solve(event_fit1$vtheta, event_fit1$theta) +
                    solve(event_prior1$vtheta, event_prior1$theta))
          event_fit1$vtheta <-
            solve(solve(event_fit1$vtheta) + solve(event_prior1$vtheta))
        } else {
          for (j in 1:ngroups) {
            if (tolower(event_model) == "piecewise exponential" &&
                length(event_prior1[[j]]$piecewiseSurvivalTime) >
                length(piecewiseSurvivalTime)) {
              i = (length(piecewiseSurvivalTime) + 1):
                length(event_prior1[[j]]$piecewiseSurvivalTime)
              event_fit1[[j]]$theta =
                c(event_fit1[[j]]$theta, event_prior1[[j]]$theta[i])
              event_fit1[[j]]$vtheta = as.matrix(Matrix::bdiag(
                event_fit1[[j]]$vtheta, event_prior1[[j]]$vtheta[i,i]*1e8))
              event_fit1[[j]]$piecewiseSurvivalTime =
                event_prior1[[j]]$piecewiseSurvivalTime
            }

            event_fit1[[j]]$theta <-
              solve(solve(event_fit1[[j]]$vtheta) +
                      solve(event_prior1[[j]]$vtheta),
                    solve(event_fit1[[j]]$vtheta, event_fit1[[j]]$theta) +
                      solve(event_prior1[[j]]$vtheta,
                            event_prior1[[j]]$theta))
            event_fit1[[j]]$vtheta <-
              solve(solve(event_fit1[[j]]$vtheta) +
                      solve(event_prior1[[j]]$vtheta))
          }
        }
      }


      # whether to include dropout model
      if (tolower(dropout_model) != "none") {

        # convert prior by treatment to prior overall
        if (!is.null(dropout_prior) && !by_treatment &&
            !("model" %in% names(dropout_prior))) {

          m = length(dropout_prior)
          w = sapply(dropout_prior, function(sub_list) sub_list$w)
          w = w/sum(w)

          # check model consistency across treatments
          model = tolower(dropout_prior[[1]]$model)
          if (m > 1) {
            for (j in 2:m) {
              if (tolower(dropout_prior[[j]]$model) != model) {
                stop("Prior dropout model must be equal across treatments.")
              }
            }

            if (model == "piecewise exponential") {
              for (j in 2:m) {
                if (!all.equal(dropout_prior[[j]]$piecewiseDropoutTime,
                               dropout_prior[[1]]$piecewiseDropoutTime)) {
                  stop(paste("piecewiseDropoutTime must be equal",
                             "across treatments."))
                }
              }
            }
          }


          if (model == "exponential") {
            # match the overall mean
            theta = sapply(dropout_prior, function(sub_list) sub_list$theta)
            lambda = exp(theta)
            lambda1 = 1/sum(w/lambda)  # hazard rate for pooled
            theta1 = log(lambda1)

            # use delta-method to obtain the variance
            vtheta1 = 0
            for (i in 1:m) {
              vtheta1 = vtheta1 + (w[i]/lambda[i])^2*dropout_prior[[i]]$vtheta
            }
            vtheta1 = vtheta1*lambda1^2

            dropout_prior1 <- list(
              model = model, theta = theta1, vtheta = vtheta1)
          } else if (model == "weibull") {
            # match the overall mean and variance

            # mean and variance of weibull as a function of theta
            fmweibull <- function(theta) {
              shape = exp(-theta[2])
              scale = exp(theta[1])
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
            theta = lapply(dropout_prior, function(sub_list) sub_list$theta)
            m1 = lapply(theta, fmweibull)
            m1mean = sapply(m1, function(sub_list) sub_list$mean)
            m1var = sapply(m1, function(sub_list) sub_list$var)

            # mean and variance for pooled
            m2 = list(mean = sum(w*m1mean),
                      var = sum(w*m1var) +
                        sum(w*m1mean^2) - (sum(w*m1mean))^2)

            # solve for theta given the mean and variance for pooled
            theta2s = sapply(theta, function(sub_list) sub_list[2])

            theta12 = uniroot(function(x)
              lgamma(1+2/exp(-x)) - 2*lgamma(1+1/exp(-x)) -
                log(m2$var/m2$mean^2 + 1),
              c(min(theta2s) - 1, max(theta2s) + 1), extendInt = "yes")$root

            theta11 = log(m2$mean) - lgamma(1+1/exp(-theta12))
            theta1 = c(theta11, theta12)

            # gradient of theta with respect to mean and variance for pooled
            ig = solve(gmweibull(theta1))

            # variance of theta for pooled
            vtheta1 = 0
            for (i in 1:m) {
              gi = gmweibull(dropout_prior[[i]]$theta)
              vm1i = gi * dropout_prior[[i]]$vtheta * t(gi)
              li = w[i]*matrix(c(1, 2*(m1[[i]]$mean - m2$mean), 0, 1), ncol=2)
              vtheta1 = vtheta1 + li %*% vm1i %*% t(li)
            }
            vtheta1 = ig %*% vtheta1 %*% t(ig)

            dropout_prior1 <- list(
              model = model, theta = theta1, vtheta = vtheta1)
          } else if (model == "log-logistic") {
            # since the mean and variance of log-logistic distribution
            # may not exist, match the cdf at the weighted average of
            # treatment-specific 97.5% percentiles and medians

            fllogis <- function(theta) {
              k = length(theta)/2
              mus = theta[seq(1,2*k-1,2)]
              sigmas = exp(theta[seq(2,2*k,2)])
              t1 = sum(w*exp(mus + qlogis(0.975)*sigmas))
              t2 = sum(w*exp(mus))
              a1 = log(1/sum(w*plogis(-(log(t1) - mus)/sigmas)) - 1)
              a2 = log(1/sum(w*plogis(-(log(t2) - mus)/sigmas)) - 1)
              gamma = (a1 - a2)/(log(t1) - log(t2))
              mu = log(t1) - 1/gamma*a1
              c(mu, -log(gamma))
            }

            # gradient vector
            gllogis <- function(theta) {
              g1 = numDeriv::grad(function(theta) fllogis(theta)[1], theta)
              g2 = numDeriv::grad(function(theta) fllogis(theta)[2], theta)
              matrix(c(g1, g2), nrow=2, byrow=TRUE)
            }

            # concatenating treatment-specific model parameters
            theta = lapply(dropout_prior, function(sub_list) sub_list$theta)

            # parameter and variance for the overall population
            theta1 = fllogis(theta)
            g = gllogis(theta)
            vtheta1 = 0
            for (i in 1:m) {
              gi = g[,(2*i-1):(2*i)]
              vtheta1 = vtheta1 + gi * dropout_prior[[i]]$vtheta * t(gi)
            }

            dropout_prior1 <- list(
              model = model, theta = theta1, vtheta = vtheta1)
          } else if (model == "log-normal") {
            # match the overall mean and variance

            # mean and variance of log-normal as a function of theta
            fmlnorm <- function(theta) {
              meanlog = theta[1]
              sdlog = exp(theta[2])
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
            theta = lapply(dropout_prior, function(sub_list) sub_list$theta)
            m1 = lapply(theta, fmlnorm)
            m1mean = sapply(m1, function(sub_list) sub_list$mean)
            m1var = sapply(m1, function(sub_list) sub_list$var)

            # mean and variance for pooled
            m2 = list(mean = sum(w*m1mean),
                      var = sum(w*m1var) +
                        sum(w*m1mean^2) - (sum(w*m1mean))^2)

            # solve for theta given the mean and variance for pooled
            theta12 = 0.5*log(log(m2$var/m2$mean^2 + 1))
            theta11 = log(m2$mean) - 0.5*exp(2*theta12)
            theta1 = c(theta11, theta12)

            # gradient of theta with respect to mean and variance for pooled
            ig = solve(gmlnorm(theta1))

            # variance of theta for pooled
            vtheta1 = 0
            for (i in 1:m) {
              gi = gmlnorm(dropout_prior[[i]]$theta)
              vm1i = gi * dropout_prior[[i]]$vtheta * t(gi)
              li = w[i]*matrix(c(1, 2*(m1[[i]]$mean - m2$mean), 0, 1), ncol=2)
              vtheta1 = vtheta1 + li %*% vm1i %*% t(li)
            }
            vtheta1 = ig %*% vtheta1 %*% t(ig)

            dropout_prior1 <- list(
              model = model, theta = theta1, vtheta = vtheta1)
          } else if (model == "piecewise exponential") {
            # match 1/lambda within each interval
            theta = sapply(dropout_prior, function(sub_list) sub_list$theta)
            lambda = exp(theta)
            npieces = length(dropout_prior[[1]]$piecewiseDropoutTime)

            # construct theta and vtheta piece by piece
            theta1 = rep(NA, npieces)
            vtheta1 = 0*diag(npieces)
            for (j in 1:npieces) {
              lambdaj = lambda[j,]
              lambda1j = 1/sum(w/lambdaj)
              theta1[j] = log(lambda1j)
              for (i in 1:ngroups) {
                vtheta1[j,j] = vtheta1[j,j] + (w[i]/lambdaj[i])^2 *
                  dropout_prior[[i]]$vtheta[j,j]
              }
              vtheta1[j,j] = vtheta1[j,j]*lambda1j^2
            }

            dropout_prior1 <- list(
              model = model, theta = theta1, vtheta = vtheta1,
              piecewiseDropoutTime = dropout_prior[[1]]$piecewiseDropoutTime)
          }
        } else if (!is.null(dropout_prior)) {
          dropout_prior1 <- dropout_prior
        }


        # fit the dropout model
        if ((!by_treatment && observed$c0 > 0) ||
            (by_treatment && all(sum_by_trt$c0 > 0))) {
          dropout_fit <- fitDropout(df = df, dropout_model,
                                    piecewiseDropoutTime, showplot,
                                    by_treatment)
          dropout_fit1 <- dropout_fit$dropout_fit
        } else {
          if (is.null(dropout_prior)) {
            stop("Prior must be specified if there is no dropout observed.")
          }

          dropout_fit <- list()
          dropout_fit1 <- dropout_prior1

          # inflate the variance
          if (!by_treatment) {
            dropout_fit1$vtheta <- dropout_prior1$vtheta*1e8
            dropout_fit1$bic <- NA
          } else {
            for (i in 1:ngroups) {
              dropout_fit1[[i]]$vtheta <- dropout_prior1[[i]]$vtheta*1e8
              dropout_fit1[[i]]$bic <- NA
            }
          }
        }

        # combine prior and likelihood to yield posterior
        if (!is.null(dropout_prior)) {
          if (!by_treatment) {
            # pad additional pieces with prior parameter values
            if (tolower(dropout_model) == "piecewise exponential" &&
                length(dropout_prior1$piecewiseDropoutTime) >
                length(piecewiseDropoutTime)) {
              i = (length(piecewiseDropoutTime) + 1):
                length(dropout_prior1$piecewiseDropoutTime)
              dropout_fit1$theta = c(
                dropout_fit1$theta, dropout_prior1$theta[i])
              dropout_fit1$vtheta = as.matrix(Matrix::bdiag(
                dropout_fit1$vtheta, dropout_prior1$vtheta[i,i]*1e8))
              dropout_fit1$piecewiseDropoutTime =
                dropout_prior1$piecewiseDropoutTime
            }

            dropout_fit1$theta <-
              solve(solve(dropout_fit1$vtheta) +
                      solve(dropout_prior1$vtheta),
                    solve(dropout_fit1$vtheta, dropout_fit1$theta) +
                      solve(dropout_prior1$vtheta, dropout_prior1$theta))
            dropout_fit1$vtheta <-
              solve(solve(dropout_fit1$vtheta) +
                      solve(dropout_prior1$vtheta))
          } else {
            for (j in 1:ngroups) {
              if (tolower(dropout_model) == "piecewise exponential" &&
                  length(dropout_prior1[[j]]$piecewiseDropoutTime) >
                  length(piecewiseDropoutTime)) {
                i = (length(piecewiseDropoutTime) + 1):
                  length(dropout_prior1[[j]]$piecewiseDropoutTime)
                dropout_fit1[[j]]$theta =
                  c(dropout_fit1[[j]]$theta, dropout_prior1[[j]]$theta[i])
                dropout_fit1[[j]]$vtheta = as.matrix(Matrix::bdiag(
                  dropout_fit1[[j]]$vtheta,
                  dropout_prior1[[j]]$vtheta[i,i]*1e8))
                dropout_fit1[[j]]$piecewiseDropoutTime =
                  dropout_prior1[[j]]$piecewiseDropoutTime
              }

              dropout_fit1[[j]]$theta <-
                solve(solve(dropout_fit1[[j]]$vtheta) +
                        solve(dropout_prior1[[j]]$vtheta),
                      solve(dropout_fit1[[j]]$vtheta,
                            dropout_fit1[[j]]$theta) +
                        solve(dropout_prior1[[j]]$vtheta,
                              dropout_prior1[[j]]$theta))
              dropout_fit1[[j]]$vtheta <-
                solve(solve(dropout_fit1[[j]]$vtheta) +
                        solve(dropout_prior1[[j]]$vtheta))
            }
          }
        }


        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = df, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1,
            dropout_fit = dropout_fit1,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment)
        } else {
          event_pred <- predictEvent(
            df = df, target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            dropout_fit = dropout_fit1,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment)
        }
      } else {  # no dropout model
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = df, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1,
            dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment)
        } else {
          event_pred <- predictEvent(
            df = df, target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            dropout_fit = NULL,
            fixedFollowup, followupTime, pilevel, nyears, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment)
        }
      }
    } else { # event prediction at design stage
      if (!is.null(dropout_prior)) {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_prior,
          dropout_fit = dropout_prior,
          fixedFollowup, followupTime, pilevel, nyears, nreps,
          showEnrollment, showEvent, showDropout, showOngoing,
          showsummary, showplot = FALSE, by_treatment)
      } else {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_prior,
          dropout_fit = NULL,
          fixedFollowup, followupTime, pilevel, nyears, nreps,
          showEnrollment, showEvent, showDropout, showOngoing,
          showsummary, showplot = FALSE, by_treatment)
      }
    }
  }


  # output results
  if (is.null(df)) { # design stage prediction
    if (tolower(to_predict) == "enrollment only") {
      if (showplot) print(enroll_pred$enroll_pred_plot)

      list(stage = "Design stage",
           to_predict = "Enrollment only",
           enroll_fit = enroll_prior, enroll_pred = enroll_pred)
    } else if (tolower(to_predict) == "enrollment and event") {
      if (showplot) print(event_pred$event_pred_plot)

      if (!is.null(dropout_prior)) {
        list(stage = "Design stage",
             to_predict = "Enrollment and event",
             enroll_fit = enroll_prior, enroll_pred = enroll_pred,
             event_fit = event_prior,
             dropout_fit = dropout_prior, event_pred = event_pred)
      } else {
        list(stage = "Design stage",
             to_predict = "Enrollment and event",
             enroll_fit = enroll_prior, enroll_pred = enroll_pred,
             event_fit = event_prior, event_pred = event_pred)
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


