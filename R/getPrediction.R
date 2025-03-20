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
#'   "piecewise exponential", "model averaging", "spline", or "cox".
#'   The model averaging uses the \code{exp(-bic/2)} weighting and
#'   combines Weibull and log-normal models. By default, it is set to
#'   "model averaging".
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
#' @param m The number of event time intervals to extrapolate the hazard
#'   function beyond the last observed event time.
#' @param event_prior The prior of event model parameters.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options:
#'   "none", "exponential", "Weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", "model averaging", "spline", or "cox".
#'   The model averaging uses the \code{exp(-bic/2)} weighting and
#'   combines Weibull and log-normal models. By default, it is set to
#'   "exponential".
#' @param piecewiseDropoutTime A vector that specifies the time
#'   intervals for the piecewise exponential dropout distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param k_dropout The number of inner knots of the spline dropout model of
#'   Royston and Parmar (2002). The default
#'   \code{k_dropout=0} gives a Weibull, log-logistic or log-normal model,
#'   if \code{scale_dropout} is "hazard", "odds", or "normal", respectively.
#'   The knots are chosen as equally-spaced quantiles of the log
#'   uncensored survival times. The boundary knots are chosen as the
#'   minimum and maximum log uncensored survival times.
#' @param scale_dropout If "hazard", the log cumulative hazard for dropout
#'   is modeled as a spline function. If "odds", the log cumulative odds is
#'   modeled as a spline function. If "normal", -qnorm(S(t)) is
#'   modeled as a spline function.
#' @param m_dropout The number of dropout time intervals to extrapolate
#'   the hazard function beyond the last observed dropout time.
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
#' @param target_t The target number of days after the data cutoff
#'   used to predict both the number of events and the probability
#'   of achieving the target event count.
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
#'   It is replaced with the treatment_description
#'   in the observed data if \code{df} is not \code{NULL}.
#' @param covariates_event The names of baseline covariates from the input
#'   data frame to include in the event model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#' @param event_prior_with_covariates The prior of event model parameters
#'   in the presence of covariates.
#' @param covariates_dropout The names of baseline covariates from the input
#'   data frame to include in the dropout model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#' @param dropout_prior_with_covariates The prior of dropout model
#'   parameters in the presence of covariates.
#' @param fix_parameter Whether to fix parameters at the maximum
#'   likelihood estimates when generating new data for prediction.
#'   Defaults to \code{FALSE}, in which case, parameters will be drawn from
#'   their approximate posterior distribution.
#'
#' @details
#' For the time-decay model, the mean function is
#' \eqn{\mu(t) = (\mu/\delta)(t - (1/\delta)(1 - \exp(-\delta t)))}
#' and the rate function is
#' \eqn{\lambda(t) = (\mu/\delta)(1 - \exp(-\delta t))}.
#' For the B-spline model, the daily enrollment rate is approximated as
#' \eqn{\lambda(t) = \exp(B(t)' \theta)},
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
#' For event prediction by treatment with prior information,
#' the \code{event_prior} (\code{dropout_prior}) variable should be
#' a list with one element per treatment. For each treatment, the
#' element should include \code{model} to specify the event (dropout)
#' model (exponential, weibull, log-logistic, log-normal,
#' or piecewise exponential), and \code{theta} and \code{vtheta} to
#' indicate the parameter values and the covariance matrix.
#' For the piecewise exponential event (dropout) model, the list
#' should also include \code{piecewiseSurvivalTime}
#' (\code{piecewiseDropoutTime}) to indicate the location of knots.
#' It should be noted that the model averaging, spline, and
#' cox options are not appropriate for use as prior.
#'
#' If the event prediction is not by treatment while the prior
#' information is given by treatment, then each element of
#' \code{event_prior} (\code{dropout_prior}) should also include
#' \code{w} to specify the weight of the treatment in a
#' randomization block. If the prediction is not by treatment and
#' the prior is given for the overall study, then \code{event_prior}
#' (\code{dropout_prior}) is a flat list with \code{model},
#' \code{theta}, and \code{vtheta}. For the piecewise exponential
#' event (dropout) model, it should also include
#' \code{piecewiseSurvivalTime} (\code{piecewiseDropoutTime}) to
#' indicate the location of knots.
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
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' # Event prediction after enrollment completion
#' set.seed(3000)
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
    enroll_model = "b-spline", nknots = 0, lags = 30,
    accrualTime = 0,
    enroll_prior = NULL,
    event_model = "model averaging", piecewiseSurvivalTime = 0,
    k = 0, scale = "hazard", m = 5,
    event_prior = NULL,
    dropout_model = "exponential", piecewiseDropoutTime = 0,
    k_dropout = 0, scale_dropout = "hazard", m_dropout = 5,
    dropout_prior = NULL,
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nyears = 4,
    target_t = NA, nreps = 500,
    showEnrollment = TRUE, showEvent = TRUE,
    showDropout = FALSE, showOngoing = FALSE,
    showsummary = TRUE, showplot = TRUE,
    by_treatment = FALSE, ngroups = 1, alloc = NULL,
    treatment_label = NULL,
    covariates_event = NULL,
    event_prior_with_covariates = NULL,
    covariates_dropout = NULL,
    dropout_prior_with_covariates = NULL,
    fix_parameter = FALSE) {

  if (!is.null(df)) erify::check_class(df, "data.frame")

  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))

  if (!is.na(target_n)) erify::check_n(target_n)
  if (!is.na(target_d)) erify::check_n(target_d)
  if (is.na(target_n) && is.na(target_d))
    stop("At least one of target_n and target_d must be specified")
  if (!is.na(target_n) && !is.na(target_d) && target_d > target_n)
    stop("target_d cannot exceed target_n")


  # check by_treatment, ngroups, and alloc
  erify::check_bool(by_treatment)

  if (is.null(df) && by_treatment) {
    erify::check_n(ngroups)
  }

  if (is.null(df)) by_treatment = TRUE
  if (!is.null(df)) dt <- data.table::setDT(data.table::copy(df))

  if (by_treatment) {
    if (!is.null(df)) {
      ngroups = dt[, data.table::uniqueN(get("treatment"))]
    }

    if (is.null(alloc)) {
      alloc = rep(1, ngroups)
    } else {
      if (length(alloc) != ngroups) {
        stop("length of alloc must be equal to the number of treatments")
      }

      if (any(alloc <= 0 | alloc != round(alloc))) {
        stop("elements of alloc must be positive integers")
      }
    }
  } else {
    ngroups = 1
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }

  if (!is.null(treatment_label) && length(treatment_label) != ngroups) {
    stop(paste("length of treatment_label must be equal to",
               "the number of treatments"))
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
        stop("Prior and likelihood must use the same enrollment model")
      }

      if (tolower(enroll_prior$model) == "piecewise poisson" &&
          (length(enroll_prior$accrualTime) < length(accrualTime) ||
           !all.equal(enroll_prior$accrualTime[1:length(accrualTime)],
                      accrualTime))) {
        stop(paste("accrualTime of piecewise Poisson must be a subset of",
                   "that in enroll_prior"))
      }
    }
  }


  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "model averaging", "spline", "cox"))

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
        stop("event_prior must be a list with one element per treatment")
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
          stop("Prior and likelihood must use the same event model")
        }

        if (tolower(event_prior2[[j]]$model) == "piecewise exponential" &&
            (length(event_prior2[[j]]$piecewiseSurvivalTime) <
             length(piecewiseSurvivalTime) ||
             !all.equal(event_prior2[[j]]$piecewiseSurvivalTime[
               1:length(piecewiseSurvivalTime)], piecewiseSurvivalTime))) {
          stop(paste("piecewiseSurvivalTime of piecewise exponential model",
                     "must be a subset of that in event_prior"))
        }
      }
    }
  }

  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "model averaging", "spline", "cox"))

  if (piecewiseDropoutTime[1] != 0) {
    stop("piecewiseDropoutTime must start with 0")
  }
  if (length(piecewiseDropoutTime) > 1 &&
      any(diff(piecewiseDropoutTime) <= 0)) {
    stop("piecewiseDropoutTime should be increasing")
  }

  erify::check_n(k_dropout, zero = TRUE)
  erify::check_content(tolower(scale_dropout), c("hazard", "odds", "normal"))

  # check dropout model prior
  if (!is.null(dropout_prior)) {
    erify::check_class(dropout_prior, "list")

    if (by_treatment) {
      if (length(dropout_prior) != ngroups) {
        stop("dropout_prior must be a list with one element per treatment")
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
                   "in dropout_prior"))
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
          stop("Prior and likelihood must use the same dropout model")
        }

        if (tolower(dropout_prior2[[j]]$model) == "piecewise exponential" &&
            (length(dropout_prior2[[j]]$piecewiseDropoutTime) <
             length(piecewiseDropoutTime) ||
             !all.equal(dropout_prior2[[j]]$piecewiseDropoutTime[
               1:length(piecewiseDropoutTime)], piecewiseDropoutTime))) {
          stop(paste("piecewiseDropoutTime of piecewise exponential model",
                     "must be a subset of that in dropout_prior"))
        }
      }


      if (!is.null(event_prior) && "w" %in% names(event_prior2[[j]])) {
        if (!("w" %in% names(dropout_prior2[[j]])) ||
            event_prior2[[j]]$w != dropout_prior2[[j]]$w) {
          stop("w must be equal between event_prior and dropout_prior")
        }
      }
    }
  }


  erify::check_n(m)
  erify::check_n(m_dropout)
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
  erify::check_bool(fix_parameter)

  if (!all(is.na(target_t))) {
    target_t <- target_t[!is.na(target_t)]
    if (any(target_t <= 0 | target_t > nyears*365)) {
      stop("target_t must be positive and less than nyears*365")
    }
  }


  if (!is.null(covariates_event)) {
    if (is.null(df)) {
      stop("df must be provided in the presence of covariates_event")
    }

    if (!all(covariates_event %in% colnames(dt))) {
      stop("All covariates_event must exist in df")
    }

    xnames = paste(covariates_event, collapse = "+")
    formula = as.formula(paste("~", xnames))
    x_event = model.matrix(formula, dt)  # design matrix with intercept
    q_event = ncol(x_event) - 1
  }


  # check event_prior_with_covariates
  event_prior_w_x <- event_prior_with_covariates
  if (!is.null(event_prior_w_x)) {
    if (is.null(covariates_event)) {
      stop(paste("covariates_event must be provided",
                 "for event_prior_with_covariates"))
    }

    erify::check_class(event_prior_w_x, "list")

    if (by_treatment) {
      if (length(event_prior_w_x) != ngroups) {
        stop(paste("event_prior_with_covariates must be a list with",
                   "one element per treatment"))
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(event_prior_w_x)) {
      event_prior2_w_x <- list()
      event_prior2_w_x[[1]] <- event_prior_w_x
    } else {
      event_prior2_w_x <- event_prior_w_x
    }

    for (j in 1:length(event_prior2_w_x)) {
      erify::check_content(tolower(event_prior2_w_x[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential"))

      model = tolower(event_prior2_w_x[[j]]$model)
      p = length(event_prior2_w_x[[j]]$theta)
      vtheta = event_prior2_w_x[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in event_prior_with_covariates"))
      }

      if ((model == "exponential" && p != 1 + q_event) ||
          (model == "weibull" && p != 2 + q_event) ||
          (model == "log-logistic" && p != 2 + q_event) ||
          (model == "log-normal" && p != 2 + q_event) ||
          (model == "piecewise exponential" &&
           p != length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) +
           q_event)) {
        stop(paste("Length of theta must be compatible with model",
                   "in event_prior_with_covariates"))
      }

      if (model == "piecewise exponential") {
        if (event_prior2_w_x[[j]]$piecewiseSurvivalTime[1] != 0) {
          stop(paste("piecewiseSurvivalTime must start with 0",
                     "in event_prior_with_covariates"))
        }
        if (length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) > 1 &&
            any(diff(event_prior2_w_x[[j]]$piecewiseSurvivalTime) <= 0)) {
          stop(paste("piecewiseSurvivalTime should be increasing",
                     "in event_prior_with_covariates"))
        }
      }


      if (!is.null(df)) {
        if (tolower(event_prior2_w_x[[j]]$model) != tolower(event_model)) {
          stop("Prior and likelihood must use the same event model")
        }

        if (tolower(event_prior2_w_x[[j]]$model) == "piecewise exponential" &&
            (length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) <
             length(piecewiseSurvivalTime) ||
             !all.equal(event_prior2_w_x[[j]]$piecewiseSurvivalTime[
               1:length(piecewiseSurvivalTime)], piecewiseSurvivalTime))) {
          stop(paste("piecewiseSurvivalTime of piecewise exponential model",
                     "must be a subset of that in",
                     "event_prior_with_covariates"))
        }
      }
    }
  }


  if (!is.null(covariates_dropout)) {
    if (is.null(df)) {
      stop("df must be provided in the presence of covariates_dropout")
    }

    if (!all(covariates_dropout %in% colnames(dt))) {
      stop("All covariates_dropout must exist in df")
    }

    xnames = paste(covariates_dropout, collapse = "+")
    formula = as.formula(paste("~", xnames))
    x_dropout = model.matrix(formula, dt)  # design matrix with intercept
    q_dropout = ncol(x_dropout) - 1
  }


  # check event_prior_with_covariates
  dropout_prior_w_x <- dropout_prior_with_covariates
  if (!is.null(dropout_prior_w_x)) {
    if (is.null(covariates_dropout)) {
      stop(paste("covariates_dropout must be provided for",
                 "dropout_prior_with_covariates"))
    }

    erify::check_class(dropout_prior_w_x, "list")

    if (by_treatment) {
      if (length(dropout_prior_w_x) != ngroups) {
        stop(paste("dropout_prior_with_covariates must be a list with",
                   "one element per treatment"))
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(dropout_prior_w_x)) {
      dropout_prior2_w_x <- list()
      dropout_prior2_w_x[[1]] <- dropout_prior_w_x
    } else {
      dropout_prior2_w_x <- dropout_prior_w_x
    }

    for (j in 1:length(dropout_prior2_w_x)) {
      erify::check_content(tolower(dropout_prior2_w_x[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential"))

      model = tolower(dropout_prior2_w_x[[j]]$model)
      p = length(dropout_prior2_w_x[[j]]$theta)
      vtheta = dropout_prior2_w_x[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in dropout_prior_with_covariates"))
      }

      if ((model == "exponential" && p != 1 + q_dropout) ||
          (model == "weibull" && p != 2 + q_dropout) ||
          (model == "log-logistic" && p != 2 + q_dropout) ||
          (model == "log-normal" && p != 2 + q_dropout) ||
          (model == "piecewise exponential" &&
           p != length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) +
           q_dropout)) {
        stop(paste("Length of theta must be compatible with model",
                   "in dropout_prior_with_covariates"))
      }

      if (model == "piecewise exponential") {
        if (dropout_prior2_w_x[[j]]$piecewiseDropoutTime[1] != 0) {
          stop(paste("piecewiseDropoutTime must start with 0",
                     "in dropout_prior_with_covariates"))
        }
        if (length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) > 1 &&
            any(diff(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) <= 0)) {
          stop(paste("piecewiseDropoutTime should be increasing",
                     "in dropout_prior_with_covariates"))
        }
      }


      if (!is.null(df)) {
        if (tolower(dropout_prior2_w_x[[j]]$model) !=
            tolower(dropout_model)) {
          stop("Prior and likelihood must use the same dropout model")
        }

        if (tolower(dropout_prior2_w_x[[j]]$model) == "piecewise exponential"
            && (length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) <
                length(piecewiseDropoutTime) ||
                !all.equal(dropout_prior2_w_x[[j]]$piecewiseDropoutTime[
                  1:length(piecewiseDropoutTime)], piecewiseDropoutTime))) {
          stop(paste("piecewiseDropoutTime of piecewise exponential model",
                     "must be a subset of that in",
                     "dropout_prior_with_covariates"))
        }
      }
    }
  }


  # check input data set to ensure it has all the required columns
  if (!is.null(df)) {
    cols = colnames(dt)

    if (tolower(to_predict) == "enrollment only") {
      req_cols = c("trialsdt", "usubjid", "randdt", "cutoffdt")
    } else {
      req_cols = c("trialsdt", "usubjid", "randdt", "time", "event",
                   "dropout", "cutoffdt")
    }

    if (by_treatment) {
      req_cols <- c(req_cols, "treatment")
    }

    if (!all(req_cols %in% cols)) {
      stop(paste("The following columns are missing from df:",
                 paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
    }

    if (any(is.na(dt[, mget(req_cols)]))) {
      stop(paste("The following columns of df have missing values:",
                 paste(req_cols[sapply(dt, function(x) any(is.na(x)))],
                       collapse = ", ")))
    }

    if ("treatment" %in% cols && !("treatment_description" %in% cols)) {
      dt[, `:=`(treatment_description = paste("Treatment", get("treatment")))]
    }
  }


  if (!is.null(df)) {
    dt$trialsdt <- as.Date(dt$trialsdt)
    dt$randdt <- as.Date(dt$randdt)
    dt$cutoffdt <- as.Date(dt$cutoffdt)

    trialsdt = dt[1, get("trialsdt")]
    cutoffdt = dt[1, get("cutoffdt")]

    # summarize observed data
    observed <- summarizeObserved(dt, to_predict, showplot, by_treatment)
  }

  if (!is.null(covariates_event)) {
    covariates_event <- tolower(covariates_event)
  }

  if (!is.null(covariates_dropout)) {
    covariates_dropout <- tolower(covariates_dropout)
  }


  # fit and predict enrollment
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) {
      erify::check_n(target_n - observed$n0,
                     supplement = "Enrollment target reached.")

      enroll_fit <- fitEnrollment(df = dt, enroll_model,
                                  nknots, accrualTime, showplot)
      enroll_fit1 <- enroll_fit$fit

      # combine prior and likelihood to yield posterior
      if (!is.null(enroll_prior)) {
        if (tolower(enroll_model) == "piecewise poisson" &&
            length(enroll_prior$accrualTime) > length(accrualTime)) {
          l1 <- length(accrualTime)

          # information from prior
          info = solve(enroll_prior$vtheta)
          # add information from data
          info[1:l1,1:l1] = info[1:l1,1:l1] +
            solve(enroll_fit1$vtheta)

          mu = solve(enroll_prior$vtheta,
                     enroll_prior$theta)
          mu[1:l1] = mu[1:l1] + solve(enroll_fit1$vtheta,
                                      enroll_fit1$theta)
        } else {
          info = solve(enroll_prior$vtheta) +
            solve(enroll_fit1$vtheta)
          mu = solve(enroll_prior$vtheta,
                     enroll_prior$theta) +
            solve(enroll_fit1$vtheta,
                  enroll_fit1$theta)
        }

        enroll_fit1 <- enroll_prior
        enroll_fit1$vtheta = solve(info)
        enroll_fit1$theta = as.numeric(enroll_fit1$vtheta %*% mu)
      }

      # enrollment prediction at the analysis stage
      enroll_pred <- predictEnrollment(
        df = dt, target_n, enroll_fit = enroll_fit1,
        lags, pilevel, nyears, nreps, showsummary, showplot = FALSE,
        by_treatment, ngroups, alloc, treatment_label,
        fix_parameter)
    } else {
      # enrollment prediction at the design stage
      enroll_pred <- predictEnrollment(
        df = NULL, target_n, enroll_fit = enroll_prior,
        lags, pilevel, nyears, nreps, showsummary, showplot = FALSE,
        by_treatment, ngroups, alloc, treatment_label,
        fix_parameter)
    }
  }


  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) { # event prediction at analysis stage
      erify::check_n(target_d - observed$d0,
                     supplement = "Event target reached.")

      if (by_treatment) {
        sum_by_trt <- dt[, list(
          n0 = .N, d0 = sum(get("event")), c0 = sum(get("dropout")),
          r0 = sum(!(get("event") | get("dropout")))),
          by = "treatment"]
      }

      # convert prior by treatment to prior overall
      if (!is.null(event_prior) && !by_treatment &&
          !("model" %in% names(event_prior))) {

        m0 = length(event_prior)
        w = sapply(event_prior, function(sub_list) sub_list$w)
        if (any(w <= 0)) {
          stop("w must be positive in event_prior")
        }
        w = w/sum(w)

        # check model consistency across treatments
        model = tolower(event_prior[[1]]$model)
        if (m0 > 1) {
          for (j in 2:m0) {
            if (tolower(event_prior[[j]]$model) != model) {
              stop(paste("Event model must be the same across",
                         "treatments in event_prior"))
            }
          }

          if (model == "piecewise exponential") {
            for (j in 2:m0) {
              if (!all.equal(event_prior[[j]]$piecewiseSurvivalTime,
                             event_prior[[1]]$piecewiseSurvivalTime)) {
                stop(paste("piecewiseSurvivalTime must be equal across",
                           "treatments in event_prior"))
              }
            }
          }
        }


        # prior mean
        theta = 0
        for (j in 1:m0) {
          theta = theta + w[j]*event_prior[[j]]$theta
        }

        # prior variance
        vtheta = 0
        for (j in 1:m0) {
          vtheta = vtheta + w[j]*(event_prior[[j]]$vtheta +
                                    event_prior[[j]]$theta %*%
                                    t(event_prior[[j]]$theta))
        }
        vtheta = vtheta - theta %*% t(theta)
        if (length(theta) == 1) vtheta = as.numeric(vtheta)

        if (model %in% c("exponential", "weibull", "log-logistic",
                         "log-normal")) {
          event_prior1 <- list(
            model = model, theta = theta, vtheta = vtheta)
        } else if (model == "piecewise exponential") {
          event_prior1 <- list(
            model = model, theta = theta, vtheta = vtheta,
            piecewiseSurvivalTime =
              event_prior[[1]]$piecewiseSurvivalTime)
        }
      } else if (!is.null(event_prior)) {
        event_prior1 <- event_prior
      }


      if (!is.null(event_prior_w_x) && !by_treatment &&
          !("model" %in% names(event_prior_w_x))) {

        m0 = length(event_prior_w_x)
        w = sapply(event_prior_w_x, function(sub_list) sub_list$w)
        if (any(w <= 0)) {
          stop("w must be positive in event_prior_with_covariates")
        }
        w = w/sum(w)

        # check model consistency across treatments
        model = tolower(event_prior_w_x[[1]]$model)
        if (m0 > 1) {
          for (j in 2:m0) {
            if (tolower(event_prior_w_x[[j]]$model) != model) {
              stop(paste("Event model must be equal across",
                         "treatments in event_prior_with_covariates"))
            }
          }

          if (model == "piecewise exponential") {
            for (j in 2:m0) {
              if (!all.equal(event_prior_w_x[[j]]$piecewiseSurvivalTime,
                             event_prior_w_x[[1]]$piecewiseSurvivalTime)) {
                stop(paste("piecewiseSurvivalTime must be equal across",
                           "treatments in event_prior_with_covariates"))
              }
            }
          }
        }


        # prior mean
        theta = 0
        for (j in 1:m0) {
          theta = theta + w[j]*event_prior_w_x[[j]]$theta
        }

        # prior variance
        vtheta = 0
        for (j in 1:m0) {
          vtheta = vtheta + w[j]*(event_prior_w_x[[j]]$vtheta +
                                    event_prior_w_x[[j]]$theta %*%
                                    t(event_prior_w_x[[j]]$theta))
        }
        vtheta = vtheta - theta %*% t(theta)

        if (model %in% c("exponential", "weibull", "log-logistic",
                         "log-normal")) {
          event_prior1_w_x <- list(
            model = model, theta = theta, vtheta = vtheta)
        } else if (model == "piecewise exponential") {
          event_prior1_w_x <- list(
            model = model, theta = theta, vtheta = vtheta,
            piecewiseSurvivalTime =
              event_prior_w_x[[1]]$piecewiseSurvivalTime)
        }
      } else if (!is.null(event_prior_w_x)) {
        event_prior1_w_x <- event_prior_w_x
      }


      # fit the event model without covariates
      if (!(to_predict == "event only" && !is.null(covariates_event))) {
        event_fit <- fitEvent(dt, event_model,
                              piecewiseSurvivalTime,
                              k, scale, m,
                              showplot, by_treatment)

        if (!by_treatment) {
          event_fit1 <- list()
          event_fit1[[1]] <- event_fit$fit
        } else {
          event_fit1 <- purrr::map(event_fit, function(fit) fit$fit)
        }

        # combine prior and likelihood to yield posterior
        if (!is.null(event_prior)) {
          if (!by_treatment) {
            event_prior2 <- list()
            event_prior2[[1]] <- event_prior1
          } else {
            event_prior2 <- event_prior1
          }

          for (j in 1:ngroups) {
            if (tolower(event_model) == "piecewise exponential" &&
                length(event_prior2[[j]]$piecewiseSurvivalTime) >
                length(piecewiseSurvivalTime)) {
              l1 <- length(piecewiseSurvivalTime)

              # information from prior
              info = solve(event_prior2[[j]]$vtheta)
              # add information from data
              info[1:l1,1:l1] = info[1:l1,1:l1] +
                solve(event_fit1[[j]]$vtheta)

              mu = solve(event_prior2[[j]]$vtheta,
                         event_prior2[[j]]$theta)
              mu[1:l1] = mu[1:l1] + solve(event_fit1[[j]]$vtheta,
                                          event_fit1[[j]]$theta)
            } else {
              info = solve(event_prior2[[j]]$vtheta) +
                solve(event_fit1[[j]]$vtheta)
              mu = solve(event_prior2[[j]]$vtheta,
                         event_prior2[[j]]$theta) +
                solve(event_fit1[[j]]$vtheta,
                      event_fit1[[j]]$theta)
            }

            event_fit2 <- event_prior2[[j]]
            event_fit2$vtheta = solve(info)
            event_fit2$theta = as.numeric(event_fit2$vtheta %*% mu)

            event_fit1[[j]] <- event_fit2
          }
        }

        if (!by_treatment) {
          event_fit1 <- event_fit1[[1]]
        }
      } else {
        event_fit1 <- NULL
      }


      # fit the event model with covariates
      if (!is.null(covariates_event)) {
        event_fit_w_x <- fitEvent(dt, event_model,
                                  piecewiseSurvivalTime,
                                  k, scale, m,
                                  showplot, by_treatment,
                                  covariates_event)

        if (!by_treatment) {
          event_fit1_w_x <- list()
          event_fit1_w_x[[1]] <- event_fit_w_x$fit
        } else {
          event_fit1_w_x <- purrr::map(event_fit_w_x,
                                       function(fit) fit$fit)
        }

        # combine prior and likelihood to yield posterior
        if (!is.null(event_prior_w_x)) {
          if (!by_treatment) {
            event_prior2_w_x <- list()
            event_prior2_w_x[[1]] <- event_prior1_w_x
          } else {
            event_prior2_w_x <- event_prior1_w_x
          }

          for (j in 1:ngroups) {
            if (tolower(event_model) == "piecewise exponential" &&
                length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) >
                length(piecewiseSurvivalTime)) {
              l1 <- length(piecewiseSurvivalTime)

              # information from prior
              info = solve(event_prior2_w_x[[j]]$vtheta)
              # add information from data
              info[1:l1,1:l1] = info[1:l1,1:l1] +
                solve(event_fit1_w_x[[j]]$vtheta)

              mu = solve(event_prior2_w_x[[j]]$vtheta,
                         event_prior2_w_x[[j]]$theta)
              mu[1:l1] = mu[1:l1] + solve(event_fit1_w_x[[j]]$vtheta,
                                          event_fit1_w_x[[j]]$theta)
            } else {
              info = solve(event_prior2_w_x[[j]]$vtheta) +
                solve(event_fit1_w_x[[j]]$vtheta)
              mu = solve(event_prior2_w_x[[j]]$vtheta,
                         event_prior2_w_x[[j]]$theta) +
                solve(event_fit1_w_x[[j]]$vtheta,
                      event_fit1_w_x[[j]]$theta)
            }

            event_fit2_w_x <- event_prior2_w_x[[j]]
            event_fit2_w_x$vtheta = solve(info)
            event_fit2_w_x$theta = as.numeric(event_fit2_w_x$vtheta
                                              %*% mu)

            event_fit1_w_x[[j]] <- event_fit2_w_x
          }
        }

        if (!by_treatment) {
          event_fit1_w_x <- event_fit1_w_x[[1]]
        }
      } else {
        event_fit1_w_x <- NULL
      }


      # whether to include dropout model
      if (tolower(dropout_model) != "none") {

        # convert prior by treatment to prior overall
        if (!is.null(dropout_prior) && !by_treatment &&
            !("model" %in% names(dropout_prior))) {

          m0 = length(dropout_prior)
          w = sapply(dropout_prior, function(sub_list) sub_list$w)
          if (any(w <= 0)) {
            stop("w must be positive in dropout_prior")
          }
          w = w/sum(w)

          # check model consistency across treatments
          model = tolower(dropout_prior[[1]]$model)
          if (m0 > 1) {
            for (j in 2:m0) {
              if (tolower(dropout_prior[[j]]$model) != model) {
                stop(paste("Dropout model must be the same across",
                           "treatments in dropout_prior"))
              }
            }

            if (model == "piecewise exponential") {
              for (j in 2:m0) {
                if (!all.equal(dropout_prior[[j]]$piecewiseDropoutTime,
                               dropout_prior[[1]]$piecewiseDropoutTime)) {
                  stop(paste("piecewiseDropoutTime must be equal across",
                             "treatments in dropout_prior"))
                }
              }
            }
          }


          # prior mean
          theta = 0
          for (j in 1:m0) {
            theta = theta + w[j]*dropout_prior[[j]]$theta
          }

          # prior variance
          vtheta = 0
          for (j in 1:m0) {
            vtheta = vtheta + w[j]*(dropout_prior[[j]]$vtheta +
                                      dropout_prior[[j]]$theta %*%
                                      t(dropout_prior[[j]]$theta))
          }
          vtheta = vtheta - theta %*% t(theta)
          if (length(theta) == 1) vtheta = as.numeric(vtheta)

          if (model %in% c("exponential", "weibull", "log-logistic",
                           "log-normal")) {
            dropout_prior1 <- list(
              model = model, theta = theta, vtheta = vtheta)
          } else if (model == "piecewise exponential") {
            dropout_prior1 <- list(
              model = model, theta = theta, vtheta = vtheta,
              piecewiseDropoutTime =
                dropout_prior[[1]]$piecewiseDropoutTime)
          }
        } else if (!is.null(dropout_prior)) {
          dropout_prior1 <- dropout_prior
        }


        if (!is.null(dropout_prior_w_x) && !by_treatment &&
            !("model" %in% names(dropout_prior_w_x))) {

          m0 = length(dropout_prior_w_x)
          w = sapply(dropout_prior_w_x, function(sub_list) sub_list$w)
          if (any(w <= 0)) {
            stop("w must be positive in dropout_prior_with_covariates")
          }
          w = w/sum(w)

          # check model consistency across treatments
          model = tolower(dropout_prior_w_x[[1]]$model)
          if (m0 > 1) {
            for (j in 2:m0) {
              if (tolower(dropout_prior_w_x[[j]]$model) != model) {
                stop(paste("Dropout model must be the same across",
                           "treatments in dropout_prior_with_covariates"))
              }
            }

            if (model == "piecewise exponential") {
              for (j in 2:m0) {
                if (!all.equal(dropout_prior_w_x[[j]]$piecewiseDropoutTime,
                               dropout_prior_w_x[[1]]$piecewiseDropoutTime)) {
                  stop(paste("piecewiseDropoutTime must be equal across",
                             "treatments in dropout_prior_with_covariates"))
                }
              }
            }
          }


          # prior mean
          theta = 0
          for (j in 1:m0) {
            theta = theta + w[j]*dropout_prior_w_x[[j]]$theta
          }

          # prior variance
          vtheta = 0
          for (j in 1:m0) {
            vtheta = vtheta + w[j]*(dropout_prior_w_x[[j]]$vtheta +
                                      dropout_prior_w_x[[j]]$theta %*%
                                      t(dropout_prior_w_x[[j]]$theta))
          }
          vtheta = vtheta - theta %*% t(theta)

          if (model %in% c("exponential", "weibull", "log-logistic",
                           "log-normal")) {
            dropout_prior1_w_x <- list(
              model = model, theta = theta, vtheta = vtheta)
          } else if (model == "piecewise exponential") {
            dropout_prior1_w_x <- list(
              model = model, theta = theta, vtheta = vtheta,
              piecewiseDropoutTime =
                dropout_prior_w_x[[1]]$piecewiseDropoutTime)
          }
        } else if (!is.null(dropout_prior_w_x)) {
          dropout_prior1_w_x <- dropout_prior_w_x
        }


        # fit the dropout model without covariates
        if (!(to_predict == "event only" && !is.null(covariates_dropout))) {
          dropout_fit <- fitDropout(dt, dropout_model,
                                    piecewiseDropoutTime,
                                    k_dropout, scale_dropout, m_dropout,
                                    showplot, by_treatment)

          if (!by_treatment) {
            dropout_fit1 <- list()
            dropout_fit1[[1]] <- dropout_fit$fit
          } else {
            dropout_fit1 <- purrr::map(dropout_fit, function(fit) fit$fit)
          }

          # combine prior and likelihood to yield posterior
          if (!is.null(dropout_prior)) {
            if (!by_treatment) {
              dropout_prior2 <- list()
              dropout_prior2[[1]] <- dropout_prior1
            } else {
              dropout_prior2 <- dropout_prior1
            }

            for (j in 1:ngroups) {
              if (tolower(dropout_model) == "piecewise exponential" &&
                  length(dropout_prior2[[j]]$piecewiseDropoutTime) >
                  length(piecewiseDropoutTime)) {
                l1 <- length(piecewiseDropoutTime)

                # information from prior
                info = solve(dropout_prior2[[j]]$vtheta)
                # add information from data
                info[1:l1,1:l1] = info[1:l1,1:l1] +
                  solve(dropout_fit1[[j]]$vtheta)

                mu = solve(dropout_prior2[[j]]$vtheta,
                           dropout_prior2[[j]]$theta)
                mu[1:l1] = mu[1:l1] + solve(dropout_fit1[[j]]$vtheta,
                                            dropout_fit1[[j]]$theta)
              } else {
                info = solve(dropout_prior2[[j]]$vtheta) +
                  solve(dropout_fit1[[j]]$vtheta)
                mu = solve(dropout_prior2[[j]]$vtheta,
                           dropout_prior2[[j]]$theta) +
                  solve(dropout_fit1[[j]]$vtheta,
                        dropout_fit1[[j]]$theta)
              }

              dropout_fit2 <- dropout_prior2[[j]]
              dropout_fit2$vtheta = solve(info)
              dropout_fit2$theta = as.numeric(dropout_fit2$vtheta %*% mu)

              dropout_fit1[[j]] <- dropout_fit2
            }
          }

          if (!by_treatment) {
            dropout_fit1 <- dropout_fit1[[1]]
          }
        } else {
          dropout_fit1 <- NULL
        }


        # fit the dropout model with covariates
        if (!is.null(covariates_dropout)) {
          dropout_fit_w_x <- fitDropout(dt, dropout_model,
                                        piecewiseDropoutTime,
                                        k_dropout, scale_dropout, m_dropout,
                                        showplot, by_treatment,
                                        covariates_dropout)

          if (!by_treatment) {
            dropout_fit1_w_x <- list()
            dropout_fit1_w_x[[1]] <- dropout_fit_w_x$fit
          } else {
            dropout_fit1_w_x <- purrr::map(dropout_fit_w_x,
                                           function(fit) fit$fit)
          }

          # combine prior and likelihood to yield posterior
          if (!is.null(dropout_prior_w_x)) {
            if (!by_treatment) {
              dropout_prior2_w_x <- list()
              dropout_prior2_w_x[[1]] <- dropout_prior1_w_x
            } else {
              dropout_prior2_w_x = dropout_prior1_w_x
            }

            for (j in 1:ngroups) {
              if (tolower(dropout_model) == "piecewise exponential" &&
                  length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) >
                  length(piecewiseDropoutTime)) {
                l1 <- length(piecewiseDropoutTime)

                # information from prior
                info = solve(dropout_prior2_w_x[[j]]$vtheta)
                # add information from data
                info[1:l1,1:l1] = info[1:l1,1:l1] +
                  solve(dropout_fit1_w_x[[j]]$vtheta)

                mu = solve(dropout_prior2_w_x[[j]]$vtheta,
                           dropout_prior2_w_x[[j]]$theta)
                mu[1:l1] = mu[1:l1] + solve(dropout_fit1_w_x[[j]]$vtheta,
                                            dropout_fit1_w_x[[j]]$theta)
              } else {
                info = solve(dropout_prior2_w_x[[j]]$vtheta) +
                  solve(dropout_fit1_w_x[[j]]$vtheta)
                mu = solve(dropout_prior2_w_x[[j]]$vtheta,
                           dropout_prior2_w_x[[j]]$theta) +
                  solve(dropout_fit1_w_x[[j]]$vtheta,
                        dropout_fit1_w_x[[j]]$theta)
              }

              dropout_fit2_w_x <- dropout_prior2_w_x[[j]]
              dropout_fit2_w_x$vtheta = solve(info)
              dropout_fit2_w_x$theta = as.numeric(dropout_fit2_w_x$vtheta
                                                  %*% mu)

              dropout_fit1_w_x[[j]] <- dropout_fit2_w_x
            }
          }

          if (!by_treatment) {
            dropout_fit1_w_x <- dropout_fit1_w_x[[1]]
          }
        } else {
          dropout_fit1_w_x <- NULL
        }


        # event prediction with a dropout model
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = dt, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1, m,
            dropout_fit = dropout_fit1, m_dropout,
            fixedFollowup, followupTime, pilevel,
            nyears, target_t, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment,
            covariates_event, event_fit1_w_x,
            covariates_dropout, dropout_fit1_w_x,
            fix_parameter)
        } else {
          event_pred <- predictEvent(
            df = dt, target_d,
            newSubjects = NULL,
            event_fit = event_fit1, m,
            dropout_fit = dropout_fit1, m_dropout,
            fixedFollowup, followupTime, pilevel,
            nyears, target_t, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment,
            covariates_event, event_fit1_w_x,
            covariates_dropout, dropout_fit1_w_x,
            fix_parameter)
        }
      } else {  # no dropout model
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = dt, target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1, m,
            dropout_fit = NULL, m_dropout,
            fixedFollowup, followupTime, pilevel,
            nyears, target_t, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment,
            covariates_event, event_fit1_w_x,
            covariates_dropout,
            dropout_fit_with_covariates = NULL,
            fix_parameter)
        } else {
          event_pred <- predictEvent(
            df = dt, target_d,
            newSubjects = NULL,
            event_fit = event_fit1, m,
            dropout_fit = NULL, m_dropout,
            fixedFollowup, followupTime, pilevel,
            nyears, target_t, nreps,
            showEnrollment, showEvent, showDropout, showOngoing,
            showsummary, showplot = FALSE, by_treatment,
            covariates_event, event_fit1_w_x,
            covariates_dropout,
            dropout_fit_with_covariates = NULL,
            fix_parameter)
        }
      }
    } else { # event prediction at design stage
      if (!is.null(dropout_prior)) {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_prior, m,
          dropout_fit = dropout_prior, m_dropout,
          fixedFollowup, followupTime, pilevel,
          nyears, target_t, nreps,
          showEnrollment, showEvent, showDropout, showOngoing,
          showsummary, showplot = FALSE, by_treatment,
          covariates_event,
          event_fit_with_covariates = event_prior_with_covariates,
          covariates_dropout,
          dropout_fit_with_covariates = dropout_prior_with_covariates,
          fix_parameter)
      } else {
        event_pred <- predictEvent(
          df = NULL, target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_prior, m,
          dropout_fit = NULL, m_dropout,
          fixedFollowup, followupTime, pilevel,
          nyears, target_t, nreps,
          showEnrollment, showEvent, showDropout, showOngoing,
          showsummary, showplot = FALSE, by_treatment,
          covariates_event,
          event_fit_with_covariates = event_prior_with_covariates,
          covariates_dropout,
          dropout_fit_with_covariates = dropout_prior_with_covariates,
          fix_parameter)
      }
    }
  }


  # obtain subject-level data from all subjects
  if (tolower(to_predict) == "enrollment only") {
    subject_data <- enroll_pred$newSubjects
    if (!is.null(df)) {
      dt[, `:=`(arrivalTime = as.numeric(get("randdt") - get("trialsdt")+1))]

      if (by_treatment) {
        subject_data <- data.table::rbindlist(list(
          dt[, `:=`(draw = 0)][
            , mget(c("draw", "usubjid", "arrivalTime",
                     "treatment", "treatment_description"))],
          subject_data), use.names = TRUE)
      } else {
        subject_data <- data.table::rbindlist(list(
          dt[, `:=`(draw = 0)][
            , mget(c("draw", "usubjid", "arrivalTime"))],
          subject_data), use.names = TRUE)
      }
    }
  } else {
    subject_data <- event_pred$newEvents
    if (!is.null(df)) {
      dt[, `:=`(
        arrivalTime = as.numeric(get("randdt") - get("trialsdt") + 1),
        totalTime = as.numeric(get("randdt") - get("trialsdt")) +
          get("time"))]

      if (by_treatment) {
        subject_data <- data.table::rbindlist(list(
          dt[get("event") | get("dropout"), `:=`(draw = 0)][
            , mget(c("draw", "usubjid", "arrivalTime", "treatment",
                     "treatment_description", "time", "event",
                     "dropout", "totalTime"))],
          subject_data), use.names = TRUE)
      } else {
        subject_data <- data.table::rbindlist(list(
          dt[get("event") | get("dropout"), `:=`(draw = 0)][
            , mget(c("draw", "usubjid", "arrivalTime", "time",
                     "event", "dropout", "totalTime"))],
          subject_data), use.names = TRUE)
      }
    }
  }

  # merge in other information such as covariates from raw data
  if (!is.null(df)) {
    varnames <- c(setdiff(names(dt), names(subject_data)), "usubjid")
    subject_data <- merge(dt[, mget(varnames)], subject_data,
                          by = "usubjid", all.y = TRUE)[order(get("draw"))]
  }


  # output results
  if (is.null(df)) { # design stage prediction
    if (tolower(to_predict) == "enrollment only") {
      if (showplot) print(enroll_pred$enroll_pred_plot)

      list(stage = "Design stage",
           to_predict = "Enrollment only",
           enroll_fit = enroll_prior, enroll_pred = enroll_pred,
           subject_data = subject_data)
    } else if (tolower(to_predict) == "enrollment and event") {
      if (showplot) print(event_pred$event_pred_plot)

      if (!is.null(dropout_prior)) {
        list(stage = "Design stage",
             to_predict = "Enrollment and event",
             enroll_fit = enroll_prior, enroll_pred = enroll_pred,
             event_fit = event_prior,
             dropout_fit = dropout_prior, event_pred = event_pred,
             subject_data = subject_data)
      } else {
        list(stage = "Design stage",
             to_predict = "Enrollment and event",
             enroll_fit = enroll_prior, enroll_pred = enroll_pred,
             event_fit = event_prior, event_pred = event_pred,
             subject_data = subject_data)
      }
    }
  } else { # analysis stage prediction
    if (tolower(to_predict) == "enrollment only") {
      if (showplot) print(enroll_pred$enroll_pred_plot)

      list(stage = "Real-time before enrollment completion",
           to_predict = "Enrollment only",
           observed = observed, enroll_fit = enroll_fit,
           enroll_pred = enroll_pred,
           subject_data = subject_data)
    } else if (tolower(to_predict) == "enrollment and event") {
      if (showplot) print(event_pred$event_pred_plot)

      if (tolower(dropout_model) != "none") {
        if (!is.null(covariates_event) &&
            !is.null(covariates_dropout)) {
          list(stage = "Real-time before enrollment completion",
               to_predict = "Enrollment and event",
               observed = observed, enroll_fit = enroll_fit,
               enroll_pred = enroll_pred,
               event_fit = event_fit,
               event_fit_with_covariates = event_fit_w_x,
               dropout_fit = dropout_fit,
               dropout_fit_with_covariates = dropout_fit_w_x,
               event_pred = event_pred,
               subject_data = subject_data)
        } else if (!is.null(covariates_event) &&
                   is.null(covariates_dropout)) {
          list(stage = "Real-time before enrollment completion",
               to_predict = "Enrollment and event",
               observed = observed, enroll_fit = enroll_fit,
               enroll_pred = enroll_pred,
               event_fit = event_fit,
               event_fit_with_covariates = event_fit_w_x,
               dropout_fit = dropout_fit,
               event_pred = event_pred,
               subject_data = subject_data)
        } else if (is.null(covariates_event) &&
                   !is.null(covariates_dropout)) {
          list(stage = "Real-time before enrollment completion",
               to_predict = "Enrollment and event",
               observed = observed, enroll_fit = enroll_fit,
               enroll_pred = enroll_pred,
               event_fit = event_fit,
               dropout_fit = dropout_fit,
               dropout_fit_with_covariates = dropout_fit_w_x,
               event_pred = event_pred,
               subject_data = subject_data)
        } else {
          list(stage = "Real-time before enrollment completion",
               to_predict = "Enrollment and event",
               observed = observed, enroll_fit = enroll_fit,
               enroll_pred = enroll_pred,
               event_fit = event_fit,
               dropout_fit = dropout_fit,
               event_pred = event_pred,
               subject_data = subject_data)
        }
      } else { # no dropout model
        if (!is.null(covariates_event)) {
          list(stage = "Real-time before enrollment completion",
               to_predict = "Enrollment and event",
               observed = observed, enroll_fit = enroll_fit,
               enroll_pred = enroll_pred,
               event_fit = event_fit,
               event_fit_with_covariates = event_fit_w_x,
               event_pred = event_pred,
               subject_data = subject_data)
        } else {
          list(stage = "Real-time before enrollment completion",
               to_predict = "Enrollment and event",
               observed = observed, enroll_fit = enroll_fit,
               enroll_pred = enroll_pred,
               event_fit = event_fit,
               event_pred = event_pred,
               subject_data = subject_data)
        }
      }
    } else if (tolower(to_predict) == "event only") {
      if (showplot) print(event_pred$event_pred_plot)

      if (tolower(dropout_model) != "none") {
        if (!is.null(covariates_event) &&
            !is.null(covariates_dropout)) {
          list(stage = "Real-time after enrollment completion",
               to_predict = "Event only",
               observed = observed,
               event_fit_with_covariates = event_fit_w_x,
               dropout_fit_with_covariates = dropout_fit_w_x,
               event_pred = event_pred,
               subject_data = subject_data)
        } else if (!is.null(covariates_event) &&
                   is.null(covariates_dropout)) {
          list(stage = "Real-time after enrollment completion",
               to_predict = "Event only",
               observed = observed,
               event_fit_with_covariates = event_fit_w_x,
               dropout_fit = dropout_fit,
               event_pred = event_pred,
               subject_data = subject_data)
        } else if (is.null(covariates_event) &&
                   !is.null(covariates_dropout)) {
          list(stage = "Real-time after enrollment completion",
               to_predict = "Event only",
               observed = observed,
               event_fit = event_fit,
               dropout_fit_with_covariates = dropout_fit_w_x,
               event_pred = event_pred,
               subject_data = subject_data)
        } else {
          list(stage = "Real-time after enrollment completion",
               to_predict = "Event only",
               observed = observed,
               event_fit = event_fit,
               dropout_fit = dropout_fit,
               event_pred = event_pred,
               subject_data = subject_data)
        }
      } else { # no dropout model
        if (!is.null(covariates_event)) {
          list(stage = "Real-time after enrollment completion",
               to_predict = "Event only",
               observed = observed,
               event_fit_with_covariates = event_fit_w_x,
               event_pred = event_pred,
               subject_data = subject_data)
        } else {
          list(stage = "Real-time after enrollment completion",
               to_predict = "Event only",
               observed = observed,
               event_fit = event_fit,
               event_pred = event_pred,
               subject_data = subject_data)
        }
      }
    }
  }
}
