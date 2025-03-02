#' @title Predict event
#' @description Utilizes pre-fitted time-to-event and time-to-dropout models
#'   to generate event and dropout times for ongoing subjects
#'   and new subjects. It also provides a
#'   prediction interval for the expected time to reach the target
#'   number of events.
#'
#' @param df The subject-level enrollment and event data, including
#'   \code{trialsdt}, \code{usubjid}, \code{randdt}, \code{cutoffdt},
#'   \code{time}, \code{event}, and \code{dropout}. The data should also
#'   include \code{treatment} coded as 1, 2, and so on, and
#'   \code{treatment_description} for by-treatment prediction. By default,
#'   it is set to \code{NULL} for event prediction at the design stage.
#' @param target_d The target number of events to reach in the study.
#' @param newSubjects The enrollment data for new subjects including
#'   \code{draw} and \code{arrivalTime}. The data should also include
#'   \code{treatment} for prediction by treatment. By default,
#'   it is set to \code{NULL},
#'   indicating the completion of subject enrollment.
#' @param event_fit The pre-fitted event model used to generate
#'   predictions.
#' @param dropout_fit The pre-fitted dropout model used to generate
#'   predictions. By default, it is set to \code{NULL},
#'   indicating no dropout.
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
#'   it is set to 500. If \code{newSubjects} is not \code{NULL},
#'   the number of draws in \code{newSubjects} should be \code{nreps}.
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
#'   show the prediction plot. By default, it is set to \code{TRUE}.
#' @param by_treatment A Boolean variable to control whether or not to
#'   predict event by treatment group. By default,
#'   it is set to \code{FALSE}.
#' @param covariates_event The names of baseline covariates from the input
#'   data frame to include in the event model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#' @param event_fit_with_covariates The pre-fitted event model with
#'   covariates used to generate event predictions for ongoing subjects.
#' @param covariates_dropout The names of baseline covariates from the input
#'   data frame to include in the dropout model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#' @param dropout_fit_with_covariates The pre-fitted dropout model with
#'   covariates used to generate dropout predictions for ongoing subjects.
#' @param fix_parameter Whether to fix parameters at the maximum
#'   likelihood estimates when generating new data for prediction.
#'   Defaults to FALSE, in which case, parameters will be drawn from
#'   their approximate posterior distribution.
#'
#' @details
#' To ensure successful event prediction at the design stage, it is
#' important to provide the \code{newSubjects} data set.
#'
#' To specify the event (dropout) model used during the design-stage event
#' prediction, the \code{event_fit} (\code{dropout_fit}) should be a list
#' with one element per treatment. For each treatment, the element
#' should include \code{model} to specify the event model
#' (exponential, weibull, log-logistic, log-normal, or piecewise
#' exponential), and \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix. For the piecewise
#' exponential event (dropout) model, the list should also include
#' \code{piecewiseSurvivalTime} (\code{piecewiseDropoutTime}) to indicate
#' the location of knots. It should be noted that the model averaging
#' and spline options are not appropriate for use during the design stage.
#'
#' Following the commencement of the trial, we obtain the event
#' model fit and the dropout model fit based on the observed data,
#' denoted as \code{event_fit} and \code{dropout_fit}, respectively.
#' These fitted models are subsequently utilized to generate event
#' and dropout times for both ongoing and new subjects in the trial.
#'
#' @return A list of prediction results which includes important
#' information such as the median, lower and upper percentiles for
#' the estimated day and date to reach the target number of events,
#' as well as simulated event data for both ongoing and new subjects.
#' The data for the prediction plot is also included
#' within this list.
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
#' @examples
#'
#' # Event prediction after enrollment completion
#'
#' event_fit <- fitEvent(df = interimData2,
#'                       event_model = "piecewise exponential",
#'                       piecewiseSurvivalTime = c(0, 140, 352))
#'
#' dropout_fit <- fitDropout(df = interimData2,
#'                           dropout_model = "exponential")
#'
#' event_pred <- predictEvent(df = interimData2, target_d = 200,
#'                            event_fit = event_fit$fit,
#'                            dropout_fit = dropout_fit$fit,
#'                            pilevel = 0.90, nreps = 100)
#'
#' @export
#'
predictEvent <- function(df = NULL, target_d, newSubjects = NULL,
                         event_fit = NULL, dropout_fit = NULL,
                         fixedFollowup = FALSE, followupTime = 365,
                         pilevel = 0.90, nyears = 4, nreps = 500,
                         showEnrollment = TRUE, showEvent = TRUE,
                         showDropout = FALSE, showOngoing = FALSE,
                         showsummary = TRUE, showplot = TRUE,
                         by_treatment = FALSE,
                         covariates_event = NULL,
                         event_fit_with_covariates = NULL,
                         covariates_dropout = NULL,
                         dropout_fit_with_covariates = NULL,
                         fix_parameter = FALSE) {

  if (!is.null(df)) erify::check_class(df, "data.frame")
  erify::check_n(target_d)
  if (!is.null(newSubjects)) erify::check_class(newSubjects, "data.frame")
  if (is.null(df) && is.null(newSubjects)) {
    stop("At least one of df and newSubjects must be specified")
  }

  erify::check_bool(by_treatment)

  if (is.null(df) && "treatment" %in% names(newSubjects))
    by_treatment = TRUE

  if (!is.null(df)) {
    setDT(df)
    setnames(df, tolower(names(df)))
  }


  # number of treatment groups and add treatment description
  if (by_treatment) {
    if (!is.null(df)) {
      if (!("treatment" %in% names(df))) {
        stop("df must contain treatment")
      }
      ngroups = df[, uniqueN(get("treatment"))]
    } else {
      ngroups = newSubjects[, uniqueN(get("treatment"))]
    }

    if (!is.null(df) && !is.null(newSubjects) &&
        length(table(df$treatment)) != length(table(newSubjects$treatment))) {
      stop("Number of treatments must match between df and newSubjects")
    }

    if (!is.null(df)) {
      if (!("treatment_description" %in% names(df))) {
        df[, `:=`(treatment_description =
                    paste("Treatment", get("treatment")))]
      }
    }

    if (!is.null(newSubjects)) {
      if (!("treatment_description" %in% names(newSubjects))) {
        newSubjects[, `:=`(treatment_description =
                             paste("Treatment", get("treatment")))]
      }
    }
  } else {  # treat as a special case of by-treatment calculation
    ngroups = 1
    if (!is.null(df)) {
      df[, `:=`(treatment = 1, treatment_description = "Overall")]
    }
    if (!is.null(newSubjects)) {
      newSubjects[, `:=`(treatment = 1, treatment_description = "Overall")]
    }
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }


  # check event_fit
  if (!is.null(event_fit)) {
    erify::check_class(event_fit, "list")

    if (!by_treatment) {  # convert event_fit to a list with 1 list element
      list1 = event_fit
      event_fit = list()
      event_fit[[1]] = list1
    }

    if (length(event_fit) != ngroups) {
      stop("event_fit must be a list with one element per treatment")
    }

    # check event_fit model
    if (!is.null(df)) {
      for (j in 1:ngroups) {
        erify::check_content(tolower(event_fit[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential",
                               "model averaging", "spline"))
      }
    } else {
      for (j in 1:ngroups) {
        erify::check_content(tolower(event_fit[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential"))
      }
    }

    # check event_fit parameters
    for (j in 1:ngroups) {
      model = tolower(event_fit[[j]]$model)
      p = length(event_fit[[j]]$theta)
      vtheta = event_fit[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with the length",
                   "of theta in event_fit"))
      }

      if ((model == "exponential" && p != 1) ||
          (model == "weibull" && p != 2) ||
          (model == "log-logistic" && p != 2) ||
          (model == "log-normal" && p != 2) ||
          (model == "piecewise exponential" &&
           p != length(event_fit[[j]]$piecewiseSurvivalTime)) ||
          (model == "model averaging" && p != 4) ||
          (model == "spline" && p != length(event_fit[[j]]$knots))) {
        stop(paste("Length of theta must be compatible with model",
                   "in event_fit"))
      }

      if (model == "piecewise exponential") {
        if (event_fit[[j]]$piecewiseSurvivalTime[1] != 0) {
          stop(paste("piecewiseSurvivalTime must start with 0",
                     "in event_fit"))
        }
        if (length(event_fit[[j]]$piecewiseSurvivalTime) > 1 &&
            any(diff(event_fit[[j]]$piecewiseSurvivalTime) <= 0)) {
          stop(paste("piecewiseSurvivalTime should be increasing",
                     "in event_fit"))
        }
      }
    }
  }


  # check dropout_fit
  if (!is.null(dropout_fit)) {
    erify::check_class(dropout_fit, "list")

    if (!by_treatment) { # convert dropout_fit to a list with 1 list element
      list1 = dropout_fit
      dropout_fit = list()
      dropout_fit[[1]] = list1
    }

    if (length(dropout_fit) != ngroups) {
      stop("dropout_fit must be a list with one element per treatment")
    }

    # check dropout_fit model
    if (!is.null(df)) {
      for (j in 1:ngroups) {
        erify::check_content(tolower(dropout_fit[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential",
                               "model averaging", "spline"))
      }
    } else {
      for (j in 1:ngroups) {
        erify::check_content(tolower(dropout_fit[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential"))
      }
    }

    # check dropout_fit parameters
    for (j in 1:ngroups) {
      model = tolower(dropout_fit[[j]]$model)
      p = length(dropout_fit[[j]]$theta)
      vtheta = dropout_fit[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with the length",
                   "of theta in dropout_fit"))
      }

      if ((model == "exponential" && p != 1) ||
          (model == "weibull" && p != 2) ||
          (model == "log-logistic" && p != 2) ||
          (model == "log-normal" && p != 2) ||
          (model == "piecewise exponential" &&
           p != length(dropout_fit[[j]]$piecewiseDropoutTime)) ||
          (model == "model averaging" && p != 4) ||
          (model == "spline" && p != length(dropout_fit[[j]]$knots))) {
        stop(paste("Length of theta must be compatible with model",
                   "in dropout_fit"))
      }

      if (model == "piecewise exponential") {
        if (dropout_fit[[j]]$piecewiseDropoutTime[1] != 0) {
          stop(paste("piecewiseDropoutTime must start with 0",
                     "in dropout_fit"))
        }
        if (length(dropout_fit[[j]]$piecewiseDropoutTime) > 1 &&
            any(diff(dropout_fit[[j]]$piecewiseDropoutTime) <= 0)) {
          stop(paste("piecewiseDropoutTime should be increasing",
                     "in dropout_fit"))
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
  erify::check_bool(fix_parameter)

  if (showplot && !(showEnrollment || showEvent ||
                    showDropout || showOngoing)) {
    stop("At least one parameter must be given for prediction plot")
  }


  # check event_fit_with_covariates
  event_fit_w_x <- event_fit_with_covariates
  if (!is.null(covariates_event)) {
    if (is.null(df)) {
      stop("df must be provided in the presence of covariates_event")
    }

    if (!all(covariates_event %in% colnames(df))) {
      stop("All covariates_event must exist in df")
    }

    xnames = paste(covariates_event, collapse = "+")
    formula = as.formula(paste("~", xnames))
    x_event = model.matrix(formula, df)  # design matrix with intercept
    q_event = ncol(x_event) - 1


    if (is.null(event_fit_w_x)) {
      stop("event_fit_with_covariates must be provided")
    }

    erify::check_class(event_fit_w_x, "list")

    if (!by_treatment) { # convert event_fit_w_x to a list with 1 element
      list1 = event_fit_w_x
      event_fit_w_x = list()
      event_fit_w_x[[1]] = list1
    }

    if (length(event_fit_w_x) != ngroups) {
      stop(paste("event_fit_with_covariates must be a list with",
                 "one element per treatment"))
    }

    # check event_fit_with_covariates model
    if (!is.null(df)) {
      for (j in 1:ngroups) {
        erify::check_content(tolower(event_fit_w_x[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential",
                               "model averaging", "spline"))
      }
    } else {
      for (j in 1:ngroups) {
        erify::check_content(tolower(event_fit_w_x[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential"))
      }
    }

    # check event_fit_with_covariates parameters
    for (j in 1:ngroups) {
      model = tolower(event_fit_w_x[[j]]$model)
      p = length(event_fit_w_x[[j]]$theta)
      vtheta = event_fit_w_x[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with the length",
                   "of theta in event_fit_with_covariates"))
      }

      if ((model == "exponential" && p != 1 + q_event) ||
          (model == "weibull" && p != 2 + q_event) ||
          (model == "log-logistic" && p != 2 + q_event) ||
          (model == "log-normal" && p != 2 + q_event) ||
          (model == "piecewise exponential" &&
           p != length(event_fit_w_x[[j]]$piecewiseSurvivalTime) +
           q_event) ||
          (model == "model averaging" && p != 2*(2 + q_event)) ||
          (model == "spline" && p != length(event_fit_w_x[[j]]$knots) +
           q_event)) {
        stop(paste("Length of theta must be compatible with model",
                   "in event_fit_with_covariates"))
      }

      if (model == "piecewise exponential") {
        if (event_fit_w_x[[j]]$piecewiseSurvivalTime[1] != 0) {
          stop(paste("piecewiseSurvivalTime must start with 0",
                     "in event_fit_with_covariates"))
        }
        if (length(event_fit_w_x[[j]]$piecewiseSurvivalTime) > 1 &&
            any(diff(event_fit_w_x[[j]]$piecewiseSurvivalTime) <= 0)) {
          stop(paste("piecewiseSurvivalTime should be increasing",
                     "in event_fit_with_covariates"))
        }
      }
    }
  }


  # check dropout_fit_with_covariates
  dropout_fit_w_x <- dropout_fit_with_covariates
  if (!is.null(covariates_dropout)) {
    if (is.null(df)) {
      stop("df must be provided in the presence of covariates_dropout")
    }

    if (!all(covariates_dropout %in% colnames(df))) {
      stop("All covariates_dropout must exist in df")
    }

    xnames = paste(covariates_dropout, collapse = "+")
    formula = as.formula(paste("~", xnames))
    x_dropout = model.matrix(formula, df)  # design matrix with intercept
    q_dropout = ncol(x_dropout) - 1


    if (is.null(dropout_fit_w_x)) {
      stop("dropout_fit_with_covariates must be provided")
    }

    erify::check_class(dropout_fit_w_x, "list")

    if (!by_treatment) { # convert dropout_fit_w_x to a list with 1 element
      list1 = dropout_fit_w_x
      dropout_fit_w_x = list()
      dropout_fit_w_x[[1]] = list1
    }

    if (length(dropout_fit_w_x) != ngroups) {
      stop(paste("dropout_fit_with_covariates must be a list with",
                 "one element per treatment"))
    }

    # check dropout_fit_with_covariates model
    if (!is.null(df)) {
      for (j in 1:ngroups) {
        erify::check_content(tolower(dropout_fit_w_x[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential",
                               "model averaging", "spline"))
      }
    } else {
      for (j in 1:ngroups) {
        erify::check_content(tolower(dropout_fit_w_x[[j]]$model),
                             c("exponential", "weibull", "log-logistic",
                               "log-normal", "piecewise exponential"))
      }
    }

    # check dropout_fit_with_covariates parameters
    for (j in 1:ngroups) {
      model = tolower(dropout_fit_w_x[[j]]$model)
      p = length(dropout_fit_w_x[[j]]$theta)
      vtheta = dropout_fit_w_x[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with the length",
                   "of theta in dropout_fit_with_covariates"))
      }

      if ((model == "exponential" && p != 1 + q_dropout) ||
          (model == "weibull" && p != 2 + q_dropout) ||
          (model == "log-logistic" && p != 2 + q_dropout) ||
          (model == "log-normal" && p != 2 + q_dropout) ||
          (model == "piecewise exponential" &&
           p != length(dropout_fit_w_x[[j]]$piecewiseDropoutTime) +
           q_dropout) ||
          (model == "model averaging" && p != 2*(2 + q_dropout)) ||
          (model == "spline" && p != length(dropout_fit_w_x[[j]]$knots) +
           q_dropout)) {
        stop(paste("Length of theta must be compatible with model",
                   "in dropout_fit_with_covariates"))
      }

      if (model == "piecewise exponential") {
        if (dropout_fit_w_x[[j]]$piecewiseDropoutTime[1] != 0) {
          stop(paste("piecewiseDropoutTime must start with 0",
                     "in dropout_fit_with_covariates"))
        }
        if (length(dropout_fit_w_x[[j]]$piecewiseDropoutTime) > 1 &&
            any(diff(dropout_fit_w_x[[j]]$piecewiseDropoutTime) <= 0)) {
          stop(paste("piecewiseDropoutTime should be increasing",
                     "in dropout_fit_with_covariates"))
        }
      }
    }
  }


  # check input data and extract ongoing subjects
  if (!is.null(df)) {
    df$trialsdt <- as.Date(df$trialsdt)
    df$randdt <- as.Date(df$randdt)
    df$cutoffdt <- as.Date(df$cutoffdt)

    trialsdt = df[1, get("trialsdt")]
    cutoffdt = df[1, get("cutoffdt")]
    t0 = as.numeric(cutoffdt - trialsdt + 1)

    if (df[, any(get("randdt") < get("trialsdt"))]) {
      stop("randdt must be greater than or equal to trialsdt")
    }

    if (df[, any(get("randdt") > get("cutoffdt"))]) {
      stop("randdt must be less than or equal to cutoffdt")
    }

    if (df[, any(get("time") < 1)]) {
      stop("time must be greater than or equal to 1")
    }

    if (df[, any(get("event") & get("dropout"))]) {
      stop("event and dropout cannot both be equal to 1 simultaneously")
    }

    if (df[, any(get("time") >
                 as.numeric(get("cutoffdt") - get("randdt") + 1))]) {
      stop("time must be less than or equal to cutoffdt - randdt + 1")
    }

    df[, `:=`(
      arrivalTime = as.numeric(get("randdt") - get("trialsdt") + 1),
      totalTime = as.numeric(get("randdt") - get("trialsdt")) + get("time"))]

    if (!("usubjid" %in% names(df))) {
      df[, `:=`(usubjid = paste0("A-", 100000 + .I))]
    }

    # subset to extract ongoing subjects
    iOngoing = df[, which(!get("event") & !get("dropout"))]
    ongoingSubjects = df[iOngoing]
    if (!is.null(covariates_event)) {
      x_eventOngoing <- x_event[iOngoing,]
    }
    if (!is.null(covariates_dropout)) {
      x_dropoutOngoing = x_dropout[iOngoing,]
    }

    usubjidOngoing <- ongoingSubjects$usubjid
    arrivalTimeOngoing <- ongoingSubjects$arrivalTime
    treatmentOngoing <- ongoingSubjects$treatment
    treatment_descriptionOngoing <- ongoingSubjects$treatment_description
    time0Ongoing <- ongoingSubjects$time
    tp = ongoingSubjects[, min(get("totalTime"))]
    cutofftpdt = as.Date(tp - 1, origin = trialsdt)
    n0 = df[, .N]
    d0 = df[, sum(get("event"))]
    c0 = df[, sum(get("dropout"))]
    r0 = ongoingSubjects[, .N]

    # subjects who have had the event or dropped out
    stoppedSubjects <- df[get("event") | get("dropout")]
  } else {
    t0 = 1
    tp = 1
    n0 = 0
    d0 = 0
    c0 = 0
    r0 = 0
  }

  t1 = t0 + 365*nyears
  d1 = target_d - d0  # number of new events
  erify::check_n(d1)


  # lower and upper percentages
  plower = (1 - pilevel)/2
  pupper = 1 - plower

  if (!is.null(newSubjects)) {
    # predicted enrollment end day
    new1 <- newSubjects[, .SD[.N], keyby = get("draw")]
    pred_day1 <- ceiling(quantile(new1$arrivalTime, c(0.5, plower, pupper)))

    # future time points at which to predict number of subjects
    t = sort(unique(c(seq(t0, t1, 30), t1, pred_day1)))
    t = t[t <= t1]
  }


  # enrollment prediction data
  df_copy <- data.table::copy(df)
  subjects_copy <- data.table::copy(newSubjects)
  if (!by_treatment) {
    if (!is.null(newSubjects)) {
      # predicted number of subjects enrolled after data cut
      subjects_copy[, `:=`(tmp_row = .I)]

      dfb1 <- CJ(
        t = t, tmp_row = subjects_copy$tmp_row, sorted = FALSE)[
          subjects_copy, on = "tmp_row"][
            , list(nenrolled = sum(get("arrivalTime") <= get("t")) + n0),
            keyby = c("t", "draw")][
              , list(n = quantile(get("nenrolled"), probs = 0.5),
                     pilevel = pilevel,
                     lower = quantile(get("nenrolled"), probs = plower),
                     upper = quantile(get("nenrolled"), probs = pupper),
                     mean = mean(get("nenrolled")),
                     var = var(get("nenrolled"))),
              keyby = "t"]
    }

    if (!is.null(df)) {
      # day 1
      df0 <- data.table(t = 1, n = 0, pilevel = pilevel,
                        lower = NA_real_, upper = NA_real_,
                        mean = 0, var = 0)

      # arrival time for subjects already enrolled before data cut
      dfa0 <- df[order(get("randdt")), list(
        t = as.numeric(get("randdt") - get("trialsdt") + 1),
        n = .I, pilevel = pilevel, lower = NA_real_, upper = NA_real_,
        mean = .I, var = 0)]

      new_row <- data.table(t = t0, n = n0, pilevel = pilevel,
                            lower = NA_real_, upper = NA_real_,
                            mean = n0, var = 0)

      dfa1 <- rbindlist(list(df0, dfa0, new_row),
                        use.names = TRUE)[, .SD[.N], by = "t"]
    }


    if (is.null(newSubjects)) { # existing subjects only
      # add predicted from data cut to specified years after data cut
      dfb1t0 <- dfa1[.N, .SD][, `:=`(
        pilevel = pilevel, lower = get("n"), upper = get("n"))]

      dfb1t1 <- data.table::copy(dfb1t0)[, `:=`(t = t1)]

      enroll_pred_df <- rbindlist(list(
        dfa1, dfb1t0, dfb1t1), use.names = TRUE)[
          , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
                parameter = "Enrollment")]
    } else if (is.null(df)) { # new subjects only
      enroll_pred_df <- dfb1[, `:=`(parameter = "Enrollment")]
    } else { # existing and new subjects
      enroll_pred_df <- rbindlist(list(
        dfa1, dfb1), use.names = TRUE)[
          , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
                 parameter = "Enrollment")]
    }

    enroll_pred_df <- enroll_pred_df[order(get("t"))]
  } else { # by treatment

    # summary of observed data by treatment
    if (!is.null(df)) {
      # add overall treatment
      df2 <- rbindlist(list(df, df_copy[, `:=`(
        treatment = 9999, treatment_description = "Overall")]),
        use.names = TRUE)

      sum_by_trt <- df2[, list(
        n0 = .I, d0 = sum(get("event")), c0 = sum(get("dropout")),
        r0 = sum(!(get("event") | get("dropout")))),
        keyby = c("treatment", "treatment_description")]
    }

    if (!is.null(newSubjects)) {
      # add overall treatment
      newSubjects2 <- rbindlist(list(
        newSubjects,
        subjects_copy[, `:=`(treatment = 9999,
                             treatment_description = "Overall")]),
        use.names = TRUE)

      if (is.null(df)) {
        sum_by_trt <- newSubjects2[
          , .SD[.N], keyby = c("treatment", "treatment_description")][
            , `:=`(n0 = 0, d0 = 0, c0 = 0, r0 = 0)][
              ,mget(c("treatment", "treatment_description",
                      "n0", "d0", "c0", "r0"))]
      }

      # predicted number of subjects enrolled by treatment after cutoff
      newSubjects2[, `:=`(tmp_row = .I)]

      dfb1 <- CJ(t = t, tmp_row = newSubjects2$tmp_row, sorted = FALSE)[
        newSubjects2, on = "tmp_row"][
          , list(nenrolled = sum(get("arrivalTime") <= get("t"))),
            keyby = c("treatment", "treatment_description", "t", "draw")][
              sum_by_trt, on = c("treatment", "treatment_description")][
                , `:=`(nenrolled = get("nenrolled") + get("n0"))][
                  , list(n = quantile(get("nenrolled"), probs = 0.5),
                         pilevel = pilevel,
                         lower = quantile(get("nenrolled"), probs = plower),
                         upper = quantile(get("nenrolled"), probs = pupper),
                         mean = mean(get("nenrolled")),
                         var = var(get("nenrolled"))),
                  keyby = c("treatment", "treatment_description", "t")]
    }


    if (!is.null(df)) {
      # day 1
      df0 <- sum_by_trt[, mget(c("treatment", "treatment_description"))][
        , `:=`(t = 1, n = 0, pilevel = pilevel, lower = NA_real_,
               upper = NA_real_, mean = 0, var = 0)]

      # arrival time for subjects already enrolled before data cut
      dfa1 <- df2[order(get("randdt")), `:=`(
        t = as.numeric(get("randdt") - get("trialsdt") + 1),
        n = .I), keyby = c("treatment", "treatment_description")][
          , `:=`(pilevel = pilevel, lower = NA_real_, upper = NA_real_,
                 mean = get("n"), var = 0)][
                   , mget(c("treatment", "treatment_description", "t", "n",
                            "pilevel", "lower", "upper", "mean", "var"))]

      sum_by_trt[, `:=`(t = t0, n = n0, pilevel = pilevel,
                        lower = NA_real_, upper = NA_real_,
                        mean = n0, var = 0)]

      dfa1 <- rbindlist(list(dfa1, df0, sum_by_trt), use.names = TRUE)[
        , .SD[.N], by = c("treatment", "treatment_description", "t")]
    }


    if (is.null(newSubjects)) { # existing subjects only
      # add predicted from data cut to specified years after data cut
      dfb1t0 <- dfa1[.N, .SD][, `:=`(
        pilevel = pilevel, lower = get("n"), upper = get("n"))]

      dfb1t1 <- data.table::copy(dfb1t0)[, `:=`(t = t1)]

      enroll_pred_df <- rbindlist(list(
        dfa1, dfb1t0, dfb1t1), use.names = TRUE)[
          , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
                 parameter = "Enrollment")]
    } else if (is.null(df)) { # new subjects only
      enroll_pred_df <- dfb1[, `:=`(parameter = "Enrollment")]
    } else { # existing and new subjects
      enroll_pred_df <- rbindlist(list(
        dfa1, dfb1), use.names = TRUE)[
          , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
                 parameter = "Enrollment")]
    }

    cols = c("treatment", "treatment_description", "t")
    enroll_pred_df <- enroll_pred_df[do.call("order", lapply(cols, as.name))]
  }


  # extract posterior draws of model parameters
  if (!is.null(event_fit)) {
    theta2 <- list()
    for (j in 1:ngroups) {
      if (!fix_parameter) {
        if (length(event_fit[[j]]$theta) == 1) {
          theta2[[j]] <- matrix(rnorm(
            nreps, mean = event_fit[[j]]$theta,
            sd = sqrt(event_fit[[j]]$vtheta)), ncol=1)
        } else {
          theta2[[j]] = mvtnorm::rmvnorm(
            nreps, mean = event_fit[[j]]$theta,
            sigma = event_fit[[j]]$vtheta)
        }
      } else {
        if (length(event_fit[[j]]$theta) == 1) {
          theta2[[j]] <- matrix(rep(event_fit[[j]]$theta, nreps), ncol=1)
        } else {
          theta2[[j]] = matrix(event_fit[[j]]$theta, nreps,
                               length(event_fit[[j]]$theta), byrow = TRUE)
        }
      }
    }
  }


  if (!is.null(dropout_fit)) {
    theta3 <- list()
    for (j in 1:ngroups) {
      if (!fix_parameter) {
        if (length(dropout_fit[[j]]$theta) == 1) {
          theta3[[j]] <- matrix(rnorm(
            nreps, mean = dropout_fit[[j]]$theta,
            sd = sqrt(dropout_fit[[j]]$vtheta)), ncol=1)
        } else {
          theta3[[j]] = mvtnorm::rmvnorm(
            nreps, mean = dropout_fit[[j]]$theta,
            sigma = dropout_fit[[j]]$vtheta)
        }
      } else {
        if (length(dropout_fit[[j]]$theta) == 1) {
          theta3[[j]] <- matrix(rep(dropout_fit[[j]]$theta, nreps), ncol=1)
        } else {
          theta3[[j]] = matrix(dropout_fit[[j]]$theta, nreps,
                               length(dropout_fit[[j]]$theta), byrow = TRUE)
        }
      }
    }
  }


  if (!is.null(event_fit_w_x)) {
    theta2_w_x <- list()
    for (j in 1:ngroups) {
      if (!fix_parameter) {
        if (length(event_fit_w_x[[j]]$theta) == 1) {
          theta2_w_x[[j]] <- matrix(rnorm(
            nreps, mean = event_fit_w_x[[j]]$theta,
            sd = sqrt(event_fit_w_x[[j]]$vtheta)), ncol=1)
        } else {
          theta2_w_x[[j]] = mvtnorm::rmvnorm(
            nreps, mean = event_fit_w_x[[j]]$theta,
            sigma = event_fit_w_x[[j]]$vtheta)
        }
      } else {
        if (length(event_fit_w_x[[j]]$theta) == 1) {
          theta2_w_x[[j]] <- matrix(rep(event_fit_w_x[[j]]$theta, nreps),
                                    ncol=1)
        } else {
          theta2_w_x[[j]] = matrix(event_fit_w_x[[j]]$theta, nreps,
                                   length(event_fit_w_x[[j]]$theta),
                                   byrow = TRUE)
        }
      }
    }
  }

  if (!is.null(dropout_fit_w_x)) {
    theta3_w_x <- list()
    for (j in 1:ngroups) {
      if (!fix_parameter) {
        if (length(dropout_fit_w_x[[j]]$theta) == 1) {
          theta3_w_x[[j]] <- matrix(rnorm(
            nreps, mean = dropout_fit_w_x[[j]]$theta,
            sd = sqrt(dropout_fit_w_x[[j]]$vtheta)), ncol=1)
        } else {
          theta3_w_x[[j]] = mvtnorm::rmvnorm(
            nreps, mean = dropout_fit_w_x[[j]]$theta,
            sigma = dropout_fit_w_x[[j]]$vtheta)
        }
      } else {
        if (length(dropout_fit_w_x[[j]]$theta) == 1) {
          theta3_w_x[[j]] <- matrix(rep(dropout_fit_w_x[[j]]$theta, nreps),
                                    ncol=1)
        } else {
          theta3_w_x[[j]] = matrix(dropout_fit_w_x[[j]]$theta, nreps,
                                   length(dropout_fit_w_x[[j]]$theta),
                                   byrow = TRUE)
        }
      }
    }
  }


  # generate the event and dropout times
  if (!is.null(newSubjects)) {
    m1 = nrow(newSubjects)
  } else {
    m1 = 0
  }

  n_rows = nreps*r0 + m1

  newEvents <- data.table(
    draw = rep(NA_real_, n_rows),
    usubjid = rep(NA_character_, n_rows),
    arrivalTime = rep(NA_real_, n_rows),
    treatment = rep(NA_real_, n_rows),
    treatment_description = rep(NA_character_, n_rows),
    time = rep(NA_real_, n_rows),
    event = rep(NA_real_, n_rows),
    dropout = rep(NA_real_, n_rows))

  offset = 0
  for (i in 1:nreps) {
    # number of new subjects in the simulated data set
    if (m1 > 0) {
      n1 = newSubjects[get("draw") == i, .N]
    } else {
      n1 = 0
    }

    m = r0 + n1

    # usubjid, arrival time, treatment, and time offset for new subjects
    if (n1 > 0) {
      newSubjects1 <- newSubjects[get("draw") == i]
      usubjidNew = newSubjects1$usubjid
      arrivalTimeNew = newSubjects1$arrivalTime
      treatmentNew = newSubjects1$treatment
      treatment_descriptionNew = newSubjects1$treatment_description
      time0New = rep(0, n1)
    }

    # concatenate ongoing and new subjects
    if (r0 == 0 && n1 > 0) {  # design stage
      usubjid = usubjidNew
      arrivalTime = arrivalTimeNew
      treatment = treatmentNew
      treatment_description = treatment_descriptionNew
      time0 = time0New
    } else if (r0 > 0 && n1 > 0) { # enrollment stage
      usubjid = c(usubjidOngoing, usubjidNew)
      arrivalTime = c(arrivalTimeOngoing, arrivalTimeNew)
      treatment = c(treatmentOngoing, treatmentNew)
      treatment_description = c(treatment_descriptionOngoing,
                                treatment_descriptionNew)
      time0 = c(time0Ongoing, time0New)
    } else if (r0 > 0 && n1 == 0) { # follow-up stage
      usubjid = usubjidOngoing
      arrivalTime = arrivalTimeOngoing
      treatment = treatmentOngoing
      treatment_description = treatment_descriptionOngoing
      time0 = time0Ongoing
    }


    # draw event time for new subjects
    if (n1 > 0) {
      survivalTimeNew = rep(NA_real_, n1)

      for (j in 1:ngroups) {
        cols = which(treatmentNew == j)
        ncols = length(cols)

        theta = theta2[[j]][i,]

        if (ncols > 0) {
          model = tolower(event_fit[[j]]$model)

          if (model == "exponential") {
            rate = exp(theta)
            survivalTimeNew[cols] = rexp(ncols, rate)
          } else if (model == "weibull") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            survivalTimeNew[cols] = rweibull(ncols, shape, scale)
          } else if (model == "log-logistic") {
            location = theta[1]
            scale = exp(theta[2])
            survivalTimeNew[cols] = exp(rlogis(ncols, location, scale))
          } else if (model == "log-normal") {
            meanlog = theta[1]
            sdlog = exp(theta[2])
            survivalTimeNew[cols] = rlnorm(ncols, meanlog, sdlog)
          } else if (model == "piecewise exponential") {
            J = length(theta)
            tcut = event_fit[[j]]$piecewiseSurvivalTime
            survivalTimeNew[cols] = qpwexp(
              runif(ncols), theta, J, tcut, lower.tail = FALSE)
          } else if (model == "model averaging") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            meanlog = theta[3]
            sdlog = exp(theta[4])
            w1 = event_fit[[j]]$w1

            # draw component indicator
            w = (runif(ncols) < w1)
            nw1 = sum(w)
            nw0 = ncols - nw1

            # draw from the corresponding component distribution
            if (nw1 > 0) {
              survivalTimeNew[cols][w==1] = rweibull(nw1, shape, scale)
            }

            if (nw0 > 0) {
              survivalTimeNew[cols][w==0] = rlnorm(nw0, meanlog, sdlog)
            }
          } else if (model == "spline") {
            knots = event_fit[[j]]$knots
            scale = event_fit[[j]]$scale

            survivalTimeNew[cols] = flexsurv::rsurvspline(
              ncols, theta, knots = knots, scale = scale)
          }
        }
      }
    }


    # draw event time for ongoing subjects without covariates
    if (r0 > 0 && is.null(event_fit_w_x)) {
      survivalTimeOngoing = rep(NA_real_, r0)

      for (j in 1:ngroups) {
        cols = which(treatmentOngoing == j)
        ncols = length(cols)

        u0 = time0[cols]
        theta = theta2[[j]][i,]

        if (ncols > 0) {
          model = tolower(event_fit[[j]]$model)

          if (model == "exponential") {
            rate = exp(theta)
            survivalTimeOngoing[cols] = rexp(ncols, rate) + u0
          } else if (model == "weibull") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            survivalTimeOngoing[cols] =
              (rexp(ncols)*scale^shape + u0^shape)^(1/shape)
          } else if (model == "log-logistic") {
            location = theta[1]
            scale = exp(theta[2])
            p = plogis(log(u0), location, scale, lower.tail = FALSE)
            survivalTimeOngoing[cols] =
              exp(qlogis(runif(ncols)*p, location, scale, lower.tail = FALSE))
          } else if (model == "log-normal") {
            meanlog = theta[1]
            sdlog = exp(theta[2])
            p = plnorm(u0, meanlog, sdlog, lower.tail = FALSE)
            survivalTimeOngoing[cols] =
              qlnorm(runif(ncols)*p, meanlog, sdlog, lower.tail = FALSE)
          } else if (model == "piecewise exponential") {
            J = length(theta)
            tcut = event_fit[[j]]$piecewiseSurvivalTime
            p = ppwexp(u0, theta, J, tcut, lower.tail = FALSE)
            survivalTimeOngoing[cols] = qpwexp(
              runif(ncols)*p, theta, J, tcut, lower.tail = FALSE)
          } else if (model == "model averaging") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            meanlog = theta[3]
            sdlog = exp(theta[4])
            w1 = event_fit[[j]]$w1

            # draw component indicator
            p1 = w1*pweibull(u0, shape, scale, lower.tail = FALSE)
            p2 = (1-w1)*plnorm(u0, meanlog, sdlog, lower.tail = FALSE)

            # p1/(p1+p2) is the posterior probability of w = 1 | T >= t0
            w = (runif(ncols) < p1/(p1+p2))
            nw1 = sum(w)
            nw0 = ncols - nw1

            # draw from the corresponding component distribution
            if (nw1 > 0) {
              survivalTimeOngoing[cols][w==1] =
                (rexp(nw1)*scale^shape + u0[w==1]^shape)^(1/shape)
            }

            if (nw0 > 0) {
              p = plnorm(u0[w==0], meanlog, sdlog, lower.tail = FALSE)
              survivalTimeOngoing[cols][w==0] =
                qlnorm(runif(nw0)*p, meanlog, sdlog, lower.tail = FALSE)
            }
          } else if (model == "spline") {
            knots = event_fit[[j]]$knots
            scale = event_fit[[j]]$scale

            p = flexsurv::psurvspline(
              u0, theta, knots = knots, scale = scale, lower.tail = FALSE)

            survivalTimeOngoing[cols] = flexsurv::qsurvspline(
              runif(ncols)*p, theta, knots = knots, scale = scale,
              lower.tail = FALSE)
          }
        }
      }
    }


    # draw event time for ongoing subjects with covariates
    if (r0 > 0 && !is.null(event_fit_w_x)) {
      survivalTimeOngoing = rep(NA_real_, r0)

      for (j in 1:ngroups) {
        cols = which(treatmentOngoing == j)
        ncols = length(cols)

        u0 = time0[cols]
        theta = theta2_w_x[[j]][i,]
        x1 = x_eventOngoing[cols,]

        if (ncols > 0) {
          model = tolower(event_fit_w_x[[j]]$model)

          if (model == "exponential") {
            rate = exp(as.numeric(x1 %*% theta))
            survivalTimeOngoing[cols] = rexp(ncols, rate) + u0
          } else if (model == "weibull") {
            shape = exp(-theta[q_event+2])
            scale = exp(as.numeric(x1 %*% theta[1:(q_event+1)]))
            survivalTimeOngoing[cols] =
              (rexp(ncols)*scale^shape + u0^shape)^(1/shape)
          } else if (model == "log-logistic") {
            location = as.numeric(x1 %*% theta[1:(q_event+1)])
            scale = exp(theta[q_event+2])
            p = plogis(log(u0), location, scale, lower.tail = FALSE)
            survivalTimeOngoing[cols] =
              exp(qlogis(runif(ncols)*p, location, scale, lower.tail = FALSE))
          } else if (model == "log-normal") {
            meanlog = as.numeric(x1 %*% theta[1:(q_event+1)])
            sdlog = exp(theta[q_event+2])
            p = plnorm(u0, meanlog, sdlog, lower.tail = FALSE)
            survivalTimeOngoing[cols] =
              qlnorm(runif(ncols)*p, meanlog, sdlog, lower.tail = FALSE)
          } else if (model == "piecewise exponential") {
            J = length(theta) - q_event    # number of intervals
            gamma = theta[1:J]
            tcut = event_fit[[j]]$piecewiseSurvivalTime
            xbeta = as.numeric(as.matrix(x1[,-1]) %*%
                                 theta[(J+1):(J+q_event)])
            p = ppwexp(u0, gamma, J, tcut, lower.tail = FALSE)
            survivalTimeOngoing[cols] = qpwexp(
              runif(ncols)^exp(-xbeta)*p, gamma, J, tcut, lower.tail = FALSE)
          } else if (model == "model averaging") {
            shape = exp(-theta[q_event+2])
            scale = exp(as.numeric(x1 %*% theta[1:(q_event+1)]))
            meanlog = as.numeric(x1 %*% theta[(q_event+3):(2*q_event+3)])
            sdlog = exp(theta[2*q_event+4])
            w1 = event_fit_w_x[[j]]$w1

            # draw component indicator
            p1 = w1*pweibull(u0, shape, scale, lower.tail = FALSE)
            p2 = (1-w1)*plnorm(u0, meanlog, sdlog, lower.tail = FALSE)

            # p1/(p1+p2) is the posterior probability of w = 1 | T >= t0
            w = (runif(ncols) < p1/(p1+p2))
            nw1 = sum(w)
            nw0 = ncols - nw1

            # draw from the corresponding component distribution
            if (nw1 > 0) {
              survivalTimeOngoing[cols][w==1] =
                (rexp(nw1)*scale[w==1]^shape + u0[w==1]^shape)^(1/shape)
            }

            if (nw0 > 0) {
              p = plnorm(u0[w==0], meanlog, sdlog, lower.tail = FALSE)
              survivalTimeOngoing[cols][w==0] =
                qlnorm(runif(nw0)*p, meanlog, sdlog, lower.tail = FALSE)
            }
          } else if (model == "spline") {
            k = length(theta) - q_event - 2
            gamma = theta[1:(k+2)]
            knots = event_fit_w_x[[j]]$knots
            scale = event_fit_w_x[[j]]$scale
            xbeta = as.numeric(as.matrix(x1[,-1]) %*%
                                 theta[(k+3):(k+q_event+2)])
            u1 = runif(ncols)

            survivalTimeOngoing[cols] = purrr::map_dbl(1:ncols, function(l) {
              p = flexsurv::psurvspline(
                u0[l], gamma, knots = knots, scale = scale,
                offset = xbeta[l], lower.tail = FALSE)
              flexsurv::qsurvspline(
                u1[l]*p, gamma, knots = knots, scale = scale,
                offset = xbeta[l], lower.tail = FALSE)
            })
          }
        }
      }
    }


    if (r0 == 0 && n1 > 0) {  # design stage
      survivalTime = survivalTimeNew
    } else if (r0 > 0 && n1 > 0) { # enrollment stage
      survivalTime = c(survivalTimeOngoing, survivalTimeNew)
    } else if (r0 > 0 && n1 == 0) { # follow-up stage
      survivalTime = survivalTimeOngoing
    }


    # draw dropout time for new subjects
    if (n1 > 0 && !is.null(dropout_fit)) {
      dropoutTimeNew = rep(NA_real_, n1)

      for (j in 1:ngroups) {
        cols = which(treatmentNew == j)
        ncols = length(cols)

        theta = theta3[[j]][i,]

        if (ncols > 0) {
          model = tolower(dropout_fit[[j]]$model)

          if (model == "exponential") {
            rate = exp(theta)
            dropoutTimeNew[cols] = rexp(ncols, rate)
          } else if (model == "weibull") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            dropoutTimeNew[cols] = rweibull(ncols, shape, scale)
          } else if (model == "log-logistic") {
            location = theta[1]
            scale = exp(theta[2])
            dropoutTimeNew[cols] = exp(rlogis(ncols, location, scale))
          } else if (model == "log-normal") {
            meanlog = theta[1]
            sdlog = exp(theta[2])
            dropoutTimeNew[cols] = rlnorm(ncols, meanlog, sdlog)
          } else if (model == "piecewise exponential") {
            J = length(theta)
            tcut = dropout_fit[[j]]$piecewiseDropoutTime
            dropoutTimeNew[cols] = qpwexp(
              runif(ncols), theta, J, tcut, lower.tail = FALSE)
          } else if (model == "model averaging") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            meanlog = theta[3]
            sdlog = exp(theta[4])
            w1 = dropout_fit[[j]]$w1

            # draw component indicator
            w = (runif(ncols) < w1)
            nw1 = sum(w)
            nw0 = ncols - nw1

            # draw from the corresponding component distribution
            if (nw1 > 0) {
              dropoutTimeNew[cols][w==1] = rweibull(nw1, shape, scale)
            }

            if (nw0 > 0) {
              dropoutTimeNew[cols][w==0] = rlnorm(nw0, meanlog, sdlog)
            }
          } else if (model == "spline") {
            knots = dropout_fit[[j]]$knots
            scale = dropout_fit[[j]]$scale

            dropoutTimeNew[cols] = flexsurv::rsurvspline(
              ncols, theta, knots = knots, scale = scale)
          }
        }
      }
    }


    # draw dropout time for ongoing subjects without covariates
    if (r0 > 0 && !is.null(dropout_fit) && is.null(dropout_fit_w_x)) {
      dropoutTimeOngoing = rep(NA_real_, r0)

      for (j in 1:ngroups) {
        cols = which(treatmentOngoing == j)
        ncols = length(cols)

        u0 = time0[cols]
        theta = theta3[[j]][i,]

        if (ncols > 0) {
          model = tolower(dropout_fit[[j]]$model)

          if (model == "exponential") {
            rate = exp(theta)
            dropoutTimeOngoing[cols] = rexp(ncols, rate) + u0
          } else if (model == "weibull") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            dropoutTimeOngoing[cols] =
              (rexp(ncols)*scale^shape + u0^shape)^(1/shape)
          } else if (model == "log-logistic") {
            location = theta[1]
            scale = exp(theta[2])
            p = plogis(log(u0), location, scale, lower.tail = FALSE)
            dropoutTimeOngoing[cols] =
              exp(qlogis(runif(ncols)*p, location, scale, lower.tail = FALSE))
          } else if (model == "log-normal") {
            meanlog = theta[1]
            sdlog = exp(theta[2])
            p = plnorm(u0, meanlog, sdlog, lower.tail = FALSE)
            dropoutTimeOngoing[cols] =
              qlnorm(runif(ncols)*p, meanlog, sdlog, lower.tail = FALSE)
          } else if (model == "piecewise exponential") {
            J = length(theta)
            tcut = dropout_fit[[j]]$piecewiseDropoutTime
            p = ppwexp(u0, theta, J, tcut, lower.tail = FALSE)
            dropoutTimeOngoing[cols] = qpwexp(
              runif(ncols)*p, theta, J, tcut, lower.tail = FALSE)
          } else if (model == "model averaging") {
            shape = exp(-theta[2])
            scale = exp(theta[1])
            meanlog = theta[3]
            sdlog = exp(theta[4])
            w1 = dropout_fit[[j]]$w1

            # draw component indicator
            p1 = w1*pweibull(u0, shape, scale, lower.tail = FALSE)
            p2 = (1-w1)*plnorm(u0, meanlog, sdlog, lower.tail = FALSE)

            # p1/(p1+p2) is the posterior probability of w = 1 | T >= t0
            w = (runif(ncols) < p1/(p1+p2))
            nw1 = sum(w)
            nw0 = ncols - nw1

            # draw from the corresponding component distribution
            if (nw1 > 0) {
              dropoutTimeOngoing[cols][w==1] =
                (rexp(nw1)*scale^shape + u0[w==1]^shape)^(1/shape)
            }

            if (nw0 > 0) {
              p = plnorm(u0[w==0], meanlog, sdlog, lower.tail = FALSE)
              dropoutTimeOngoing[cols][w==0] =
                qlnorm(runif(nw0)*p, meanlog, sdlog, lower.tail = FALSE)
            }
          } else if (model == "spline") {
            knots = dropout_fit[[j]]$knots
            scale = dropout_fit[[j]]$scale

            p = flexsurv::psurvspline(
              u0, theta, knots = knots, scale = scale, lower.tail = FALSE)

            dropoutTimeOngoing[cols] = flexsurv::qsurvspline(
              runif(ncols)*p, theta, knots = knots, scale = scale,
              lower.tail = FALSE)
          }
        }
      }
    }


    # draw dropout time for ongoing subjects with covariates
    if (r0 > 0 && !is.null(dropout_fit_w_x)) {
      dropoutTimeOngoing = rep(NA_real_, r0)

      for (j in 1:ngroups) {
        cols = which(treatmentOngoing == j)
        ncols = length(cols)

        u0 = time0[cols]
        theta = theta3_w_x[[j]][i,]
        x1 = x_dropoutOngoing[cols,]

        if (ncols > 0) {
          model = tolower(dropout_fit_w_x[[j]]$model)

          if (model == "exponential") {
            rate = exp(as.numeric(x1 %*% theta))
            dropoutTimeOngoing[cols] = rexp(ncols)/rate + u0
          } else if (model == "weibull") {
            shape = exp(-theta[q_dropout+2])
            scale = exp(as.numeric(x1 %*% theta[1:(q_dropout+1)]))
            dropoutTimeOngoing[cols] =
              (rexp(ncols)*scale^shape + u0^shape)^(1/shape)
          } else if (model == "log-logistic") {
            location = as.numeric(x1 %*% theta[1:(q_dropout+1)])
            scale = exp(theta[q_dropout+2])
            p = plogis(log(u0), location, scale, lower.tail = FALSE)
            dropoutTimeOngoing[cols] =
              exp(qlogis(runif(ncols)*p, location, scale, lower.tail = FALSE))
          } else if (model == "log-normal") {
            meanlog = as.numeric(x1 %*% theta[1:(q_dropout+1)])
            sdlog = exp(theta[q_dropout+2])
            p = plnorm(u0, meanlog, sdlog, lower.tail = FALSE)
            dropoutTimeOngoing[cols] =
              qlnorm(runif(ncols)*p, meanlog, sdlog, lower.tail = FALSE)
          } else if (model == "piecewise exponential") {
            J = length(theta) - q_dropout    # number of intervals
            gamma = theta[1:J]
            tcut = dropout_fit[[j]]$piecewiseDropoutTime
            xbeta = as.numeric(as.matrix(x1[,-1]) %*%
                                 theta[(J+1):(J+q_dropout)])
            p = ppwexp(u0, gamma, J, tcut, lower.tail = FALSE)
            dropoutTimeOngoing[cols] = qpwexp(
              runif(ncols)^exp(-xbeta)*p, gamma, J, tcut, lower.tail = FALSE)
          } else if (model == "model averaging") {
            shape = exp(-theta[q_dropout+2])
            scale = exp(as.numeric(x1 %*% theta[1:(q_dropout+1)]))
            meanlog = as.numeric(x1 %*% theta[(q_dropout+3):(2*q_dropout+3)])
            sdlog = exp(theta[2*q_dropout+4])
            w1 = dropout_fit_w_x[[j]]$w1

            # draw component indicator
            p1 = w1*pweibull(u0, shape, scale, lower.tail = FALSE)
            p2 = (1-w1)*plnorm(u0, meanlog, sdlog, lower.tail = FALSE)

            # p1/(p1+p2) is the posterior probability of w = 1 | T >= t0
            w = (runif(ncols) < p1/(p1+p2))
            nw1 = sum(w)
            nw0 = ncols - nw1

            # draw from the corresponding component distribution
            if (nw1 > 0) {
              dropoutTimeOngoing[cols][w==1] =
                (rexp(nw1)*scale[w==1]^shape + u0[w==1]^shape)^(1/shape)
            }

            if (nw0 > 0) {
              p = plnorm(u0[w==0], meanlog, sdlog, lower.tail = FALSE)
              dropoutTimeOngoing[cols][w==0] =
                qlnorm(runif(nw0)*p, meanlog, sdlog, lower.tail = FALSE)
            }
          } else if (model == "spline") {
            k = length(theta) - q_dropout - 2
            gamma = theta[1:(k+2)]
            knots = dropout_fit_w_x[[j]]$knots
            scale = dropout_fit_w_x[[j]]$scale
            xbeta = as.numeric(as.matrix(x1[,-1]) %*%
                                 theta[(k+3):(k+q_dropout+2)])
            u1 = runif(ncols)

            dropoutTimeOngoing[cols] = purrr::map_dbl(1:ncols, function(l) {
              p = flexsurv::psurvspline(
                u0[l], gamma, knots = knots, scale = scale,
                offset = xbeta[l], lower.tail = FALSE)
              flexsurv::qsurvspline(
                u1[l]*p, gamma, knots = knots, scale = scale,
                offset = xbeta[l], lower.tail = FALSE)
            })
          }
        }
      }
    }

    if (!is.null(dropout_fit) || !is.null(dropout_fit_w_x)) {
      if (r0 == 0 && n1 > 0) {  # design stage
        dropoutTime = dropoutTimeNew
      } else if (r0 > 0 && n1 > 0) { # enrollment stage
        dropoutTime = c(dropoutTimeOngoing, dropoutTimeNew)
      } else if (r0 > 0 && n1 == 0) { # follow-up stage
        dropoutTime = dropoutTimeOngoing
      }
    }


    # observed survival time and event indicator
    if (!fixedFollowup) {
      if (!is.null(dropout_fit) || !is.null(dropout_fit_w_x)) {
        time = pmin(survivalTime, dropoutTime)
        event = 1*(time == survivalTime)
        dropout = 1*(time == dropoutTime)
      } else {
        time = survivalTime
        event = 1
        dropout = 0
      }
    } else {
      if (!is.null(dropout_fit) || !is.null(dropout_fit_w_x)) {
        time = pmin(survivalTime, dropoutTime, followupTime)
        event = 1*(time == survivalTime)
        dropout = 1*(time == dropoutTime)
      } else {
        time = pmin(survivalTime, followupTime)
        event = 1*(time == survivalTime)
        dropout = 0
      }
    }

    # fill out the ith block of output data frame
    index = offset + (1:m)
    newEvents[index, "draw"] = i
    newEvents[index, "usubjid"] = usubjid
    newEvents[index, "arrivalTime"] = arrivalTime
    newEvents[index, "treatment"] = treatment
    newEvents[index, "treatment_description"] = treatment_description
    newEvents[index, "time"] = pmax(round(time), time0+1)
    newEvents[index, "event"] = event
    newEvents[index, "dropout"] = dropout
    offset = offset + m
  }


  # calculate total time since trial start
  newEvents <- newEvents[
    , `:=`(totalTime = get("arrivalTime") + get("time") - 1)]

  if (!is.null(df)) {
    # combined stopped, ongoing and new subjects
    cross_joined <- CJ(
      draw = 1:nreps, usubjid = stoppedSubjects$usubjid,
      sorted = FALSE)[stoppedSubjects, on = "usubjid",
                      mget(c("draw", "usubjid", "arrivalTime",
                             "treatment", "treatment_description",
                             "time", "event", "dropout", "totalTime"))]

    allSubjects <- rbindlist(list(cross_joined, newEvents),
                             use.names = TRUE)
  } else {
    allSubjects <- newEvents
  }

  # remove the dummy treatment from newEvents
  if (!by_treatment) {
    newEvents[, (c("treatment", "treatment_description")) := NULL]
  }

  # A general quantile method if there are data sets not reaching target_d
  # Find t such that sum(I{D_i(t) < target_d}, {i, 1, nreps}) / nreps = q.
  # This works because {D_i(t) < target_d} = {T_i(target_d) > t},
  # where D_i(t) is the cumulative number of events at time t, and
  # T_i(target_d) is the time to reach target_d for data set i.
  sdf <- function(t, target_d, d0, newEvents) {
    sumdata <- newEvents[, list(
      n = sum(get("totalTime") <= t & get("event"))), keyby = get("draw")]

    mean(sumdata$n < target_d - d0)
  }

  tmax = newEvents[get("event") == 1, max(get("totalTime"))]

  # obtain the quantiles
  if (sdf(tmax, target_d, d0, newEvents) == 0) { # target_d reached for all
    new1 <- newEvents[get("event") == 1][
      , setorder(.SD, "totalTime"), by = "draw"][
        , .SD[target_d - d0], by = "draw"]

    pred_day <- ceiling(quantile(new1$totalTime, c(0.5, plower, pupper)))
  } else {
    qs = 1 - c(0.5, plower, pupper)
    pred_day = rep(NA_real_, 3)
    for (j in 1:3) {
      # check if the quantile can be estimated from observed data
      if (sdf(tmax, target_d, d0, newEvents) <= qs[j]) {
        pred_day[j] = uniroot(function(x)
          sdf(x, target_d, d0, newEvents) - qs[j],
          c(tp, tmax), tol = 1)$root
        pred_day[j] = ceiling(pred_day[j])
      }
    }
    names(pred_day) <- names(quantile(1:100, c(0.5, plower, pupper)))
  }


  if (!is.null(df)) {
    pred_date <- as.Date(pred_day - 1, origin = trialsdt)

    str1 <- paste("Time from cutoff until", target_d, "events:",
                  pred_date[1] - cutoffdt + 1, "days")
    str2 <- paste("Median prediction date:", pred_date[1])
    str3 <- paste0("Prediction interval: ", pred_date[2], ", ", pred_date[3])
    s1 <- paste(str1, "\n", str2, "\n", str3, "\n")
  } else {
    str1 <- paste("Time from trial start until", target_d, "events")
    str2 <- paste("Median prediction day:", pred_day[1])
    str3 <- paste0("Prediction interval: ", pred_day[2], ", ", pred_day[3])
    s1 <- paste(str1, "\n", str2, "\n", str3, "\n")
  }


  # observed time points
  t2 = sort(unique(c(df$arrivalTime, df$totalTime)))

  # future time points at which to predict number of events
  t = unique(c(t2[t2 >= tp], seq(t0, t1, 30), t1))

  subjects_copy <- data.table::copy(allSubjects)
  if (!by_treatment) {
    # number of events, dropouts, and ongoing subjects after data cut
    subjects_copy[, `:=`(tmp_row = .I)]

    df1 = CJ(
      t = t, tmp_row = subjects_copy$tmp_row, sorted = FALSE)[
        subjects_copy, on = "tmp_row"][
          , list(nevents = sum(get("totalTime") <= get("t") &
                                 get("event")),
                 ndropouts = sum(get("totalTime") <= get("t") &
                                   get("dropout")),
                 nongoing = sum(get("arrivalTime") <= get("t") &
                                  get("totalTime") > get("t"))),
          keyby = c("t", "draw")]

    # predicted number of events after data cut
    dfb2 = df1[, list(n = quantile(get("nevents"), probs = 0.5),
                      pilevel = pilevel,
                      lower = quantile(get("nevents"), probs = plower),
                      upper = quantile(get("nevents"), probs = pupper),
                      mean = mean(get("nevents")),
                      var = var(get("nevents"))),
               keyby = c("t")]

    # predicted number of dropouts after data cut
    dfb3 = df1[, list(n = quantile(get("ndropouts"), probs = 0.5),
                      pilevel = pilevel,
                      lower = quantile(get("ndropouts"), probs = plower),
                      upper = quantile(get("ndropouts"), probs = pupper),
                      mean = mean(get("ndropouts")),
                      var = var(get("ndropouts"))),
               keyby = c("t")]

    # predicted number of subjects at risk after data cut
    dfb4 = df1[, list(n = quantile(get("nongoing"), probs = 0.5),
                      pilevel = pilevel,
                      lower = quantile(get("nongoing"), probs = plower),
                      upper = quantile(get("nongoing"), probs = pupper),
                      mean = mean(get("nongoing")),
                      var = var(get("nongoing"))),
               keyby = c("t")]

    if (!is.null(df)) {
      # day 1
      df0 <- data.table(t = 1, n = 0, pilevel = pilevel,
                        lower = NA_real_, upper = NA_real_,
                        mean = 0, var = 0)

      # observed number of events before data cut
      dt_copy <- df[get("totalTime") <= tp]
      setorder(dt_copy, "totalTime")

      dfa2 <- rbindlist(list(
        df0, dt_copy[, list(t = get("totalTime"),
                            n = cumsum(get("event")),
                            pilevel = pilevel,
                            lower = NA_real_,
                            upper = NA_real_,
                            mean = cumsum(get("event")),  # Reuses n
                            var = 0)]), use.names = TRUE)

      # observed number of dropouts before data cut
      dfa3 <- rbindlist(list(
        df0, dt_copy[, list(t = get("totalTime"),
                            n = cumsum(get("dropout")),
                            pilevel = pilevel,
                            lower = NA_real_,
                            upper = NA_real_,
                            mean = cumsum(get("dropout")),  # Reuses n
                            var = 0)]), use.names = TRUE)

      # observed number of ongoing subjects before data cutoff
      df_copy[, `:=`(tmp_row = .I)]

      dfa4 <- rbindlist(list(
        df0, CJ(t = t2, tmp_row = df_copy$tmp_row, sorted = FALSE)[
          df_copy, on = "tmp_row"][
            get("t") <= tp, list(
              n = sum(get("arrivalTime") <= get("t") &
                        (get("totalTime") > get("t") |
                           (!get("event") & !get("dropout")))),
              pilevel = pilevel,
              lower = NA_real_, upper = NA_real_,
              mean = sum(get("arrivalTime") <= get("t") &
                           (get("totalTime") > get("t") |
                              (!get("event") & !get("dropout")))),
              var = 0), by = "t"]), use.names = TRUE)

      # add time tp
      dfa2 <- rbindlist(list(
        dfa2, data.table::copy(dfa2)[.N][, `:=`(t = tp)]),
        use.names = TRUE)[, .SD[.N], by = "t"]

      dfa3 <- rbindlist(list(
        dfa3, data.table::copy(dfa3)[.N][, `:=`(t = tp)]),
        use.names = TRUE)[, .SD[.N], by = "t"]

      dfa4 <- rbindlist(list(
        dfa4, data.table::copy(dfa4)[.N][, `:=`(t = tp)]),
        use.names = TRUE)[, .SD[.N], by = "t"]


      # concatenate events before and after data cut
      event_pred_df <- rbindlist(list(dfa2, dfb2), use.names = TRUE)[
          , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
                 parameter = "Event")]

      # concatenate dropouts before and after data cut
      dropout_pred_df <- rbindlist(list(dfa3, dfb3), use.names = TRUE)[
          , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
                 parameter = "Dropout")]

      # concatenate ongoing subjects before and after data cut
      ongoing_pred_df <- rbindlist(list(dfa4, dfb4), use.names = TRUE)[
          , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
                 parameter = "Ongoing")]
    } else {
      event_pred_df <- dfb2[, `:=`(parameter = "Event")]
      dropout_pred_df <- dfb3[, `:=`(parameter = "Dropout")]
      ongoing_pred_df <- dfb4[, `:=`(parameter = "Ongoing")]
    }

    setorder(event_pred_df, "t")
    setorder(dropout_pred_df, "t")
    setorder(ongoing_pred_df, "t")


    # generate plot
    if (showEnrollment | showEvent | showDropout | showOngoing) {
      dt_list <- list()
      if (showEnrollment) dt_list <- c(dt_list, list(enroll_pred_df))
      if (showEvent) dt_list <- c(dt_list, list(event_pred_df))
      if (showDropout) dt_list <- c(dt_list, list(dropout_pred_df))
      if (showOngoing) dt_list <- c(dt_list, list(ongoing_pred_df))

      dfs <- rbindlist(dt_list, use.names = TRUE)[
        , `:=`(parameter = factor(get("parameter"), levels = c(
          "Enrollment", "Event", "Dropout", "Ongoing")))]

      if (!is.null(df)) { # after trial start
        dfa <- dfs[is.na(get("lower"))]
        dfb <- dfs[!is.na(get("lower"))]

        dfa_enrollment <- dfa[get("parameter") == "Enrollment"]
        dfb_enrollment <- dfb[get("parameter") == "Enrollment"]
        dfa_event <- dfa[get("parameter") == "Event"]
        dfb_event <- dfb[get("parameter") == "Event"]
        dfa_dropout <- dfa[get("parameter") == "Dropout"]
        dfb_dropout <- dfb[get("parameter") == "Dropout"]
        dfa_ongoing <- dfa[get("parameter") == "Ongoing"]
        dfb_ongoing <- dfb[get("parameter") == "Ongoing"]

        g1 <- plotly::plot_ly() %>%
          plotly::add_lines(
            data = dfa_enrollment, x = ~date, y = ~n,
            line = list(shape="hv", width=2),
            name = "observed enrollment") %>%
          plotly::add_lines(
            data = dfb_enrollment, x = ~date, y = ~n,
            line = list(width=2),
            name = "median prediction enrollment") %>%
          plotly::add_ribbons(
            data = dfb_enrollment, x = ~date, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval enrollment") %>%
          plotly::add_lines(
            data = dfa_event, x = ~date, y = ~n,
            line = list(shape="hv", width=2),
            name = "observed event") %>%
          plotly::add_lines(
            data = dfb_event, x = ~date, y = ~n,
            line = list(width=2),
            name = "median prediction event") %>%
          plotly::add_ribbons(
            data = dfb_event, x = ~date, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval event") %>%
          plotly::add_lines(
            data = dfa_dropout, x = ~date, y = ~n,
            line = list(shape="hv", width=2),
            name = "observed dropout") %>%
          plotly::add_lines(
            data = dfb_dropout, x = ~date, y = ~n,
            line = list(width=2),
            name = "median prediction dropout") %>%
          plotly::add_ribbons(
            data = dfb_dropout, x = ~date, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval dropout") %>%
          plotly::add_lines(
            data = dfa_ongoing, x = ~date, y = ~n,
            line = list(shape="hv", width=2),
            name = "observed ongoing") %>%
          plotly::add_lines(
            data = dfb_ongoing, x = ~date, y = ~n,
            line = list(width=2),
            name = "median prediction ongoing") %>%
          plotly::add_ribbons(
            data = dfb_ongoing, x = ~date, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval ongoing") %>%
          plotly::add_lines(
            x = rep(cutoffdt, 2), y = c(min(dfa$n), max(dfb$upper)),
            name = "cutoff", line = list(dash="dash"),
            showlegend = FALSE) %>%
          plotly::layout(
            annotations = list(
              x = cutoffdt, y = 0, text = 'cutoff',
              xanchor = "left", yanchor = "bottom",
              font = list(size=12), showarrow = FALSE),
            xaxis = list(title = "", zeroline = FALSE),
            yaxis = list(zeroline = FALSE))

        if (tp < t0) {
          g1 <- g1 %>%
            plotly::add_lines(
              x = rep(cutofftpdt, 2), y = c(min(dfa$n), max(dfb$upper)),
              name = "prediction start",
              line = list(dash="dash", color="grey"), showlegend = FALSE) %>%
            plotly::layout(
              annotations = list(
                x = cutofftpdt, y = 0, text = 'prediction start',
                xanchor = "left", yanchor = "bottom",
                font = list(size=12), showarrow = FALSE))
        }

        if (showEvent) {
          g1 <- g1 %>%
            plotly::add_lines(
              x = range(dfs$date), y = rep(target_d, 2),
              name = 'target events', showlegend = FALSE,
              line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
            plotly::layout(
              annotations = list(
                x = 0.95, xref = "paper", y = target_d,
                text = 'target events',
                xanchor = "right", yanchor = "bottom",
                font = list(size=12), showarrow = FALSE))
        }
      } else { # at design stage
        dfs_enrollment <- dfs[get("parameter") == "Enrollment"]
        dfs_event <- dfs[get("parameter") == "Event"]
        dfs_dropout <- dfs[get("parameter") == "Dropout"]
        dfs_ongoing <- dfs[get("parameter") == "Ongoing"]

        g1 <- plotly::plot_ly() %>%
          plotly::add_lines(
            data = dfs_enrollment, x = ~t, y = ~n,
            line = list(width=2),
            name = "median prediction enrollment") %>%
          plotly::add_ribbons(
            data = dfs_enrollment, x = ~t, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval enrollment") %>%
          plotly::add_lines(
            data = dfs_event, x = ~t, y = ~n,
            line = list(width=2),
            name = "median prediction event") %>%
          plotly::add_ribbons(
            data = dfs_event, x = ~t, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval event") %>%
          plotly::add_lines(
            data = dfs_dropout, x = ~t, y = ~n,
            line = list(width=2),
            name = "median prediction dropout") %>%
          plotly::add_ribbons(
            data = dfs_dropout, x = ~t, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval dropout") %>%
          plotly::add_lines(
            data = dfs_ongoing, x = ~t, y = ~n,
            line = list(width=2),
            name = "median prediction ongoing") %>%
          plotly::add_ribbons(
            data = dfs_ongoing, x = ~t, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", line = list(width=0),
            name = "prediction interval ongoing") %>%
          plotly::layout(
            xaxis = list(title = "Days since trial start", zeroline = FALSE),
            yaxis = list(zeroline = FALSE))

        if (showEvent) {
          g1 <- g1 %>%
            plotly::add_lines(
              x = range(dfs$t), y = rep(target_d, 2),
              name = 'target events', showlegend = FALSE,
              line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
            plotly::layout(
              annotations = list(
                x = 0.95, xref = "paper", y = target_d,
                text = 'target events',
                xanchor = "right", yanchor = "bottom",
                font = list(size=12), showarrow = FALSE))
        }
      }
    }
  } else {  # by treatment
    # add overall treatment
    allSubjects2 <- rbindlist(list(
      allSubjects, subjects_copy[, `:=`(treatment = 9999,
                                        treatment_description = "Overall")]),
      use.names = TRUE)

    # number of events, dropouts, and ongoing subjects after data cut
    allSubjects2[, `:=`(tmp_row = .I)]

    df1 = CJ(
      t = t, tmp_row = allSubjects2$tmp_row, sorted = FALSE)[
        allSubjects2, on = "tmp_row"][
          , list(nevents = sum(get("totalTime") <= get("t") &
                                 get("event")),
                 ndropouts = sum(get("totalTime") <= get("t") &
                                   get("dropout")),
                 nongoing = sum(get("arrivalTime") <= get("t") &
                                  get("totalTime") > get("t"))),
          keyby = c("treatment", "treatment_description", "t", "draw")]

    # predicted number of events after data cut
    cols <- c("treatment", "treatment_description", "t")

    dfb2 = df1[, list(n = quantile(get("nevents"), probs = 0.5),
                      pilevel = pilevel,
                      lower = quantile(get("nevents"), probs = plower),
                      upper = quantile(get("nevents"), probs = pupper),
                      mean = mean(get("nevents")),
                      var = var(get("nevents"))),
               keyby = cols]

    # predicted number of dropouts after data cut
    dfb3 = df1[, list(n = quantile(get("ndropouts"), probs = 0.5),
                      pilevel = pilevel,
                      lower = quantile(get("ndropouts"), probs = plower),
                      upper = quantile(get("ndropouts"), probs = pupper),
                      mean = mean(get("ndropouts")),
                      var = var(get("ndropouts"))),
               keyby = cols]

    # predicted number of subjects at risk after data cut
    dfb4 = df1[, list(n = quantile(get("nongoing"), probs = 0.5),
                      pilevel = pilevel,
                      lower = quantile(get("nongoing"), probs = plower),
                      upper = quantile(get("nongoing"), probs = pupper),
                      mean = mean(get("nongoing")),
                      var = var(get("nongoing"))),
               keyby = cols]


    if (!is.null(df)) {
      # day 1
      df0 <- sum_by_trt[, list(
        treatment, treatment_description, t = 1, n = 0, pilevel = pilevel,
        lower = NA_real_, upper = NA_real_, mean = 0, var = 0)]


      # observed number of events before data cut
      dt_copy <- df2[get("totalTime") <= tp]
      setorderv(dt_copy, c("treatment", "treatment_description", "totalTime"))

      dfa2 <- rbindlist(list(
        df0, dt_copy[, list(t = get("totalTime"),
                            n = cumsum(get("event")),
                            pilevel = pilevel,
                            lower = NA_real_,
                            upper = NA_real_,
                            mean = cumsum(get("event")),  # Reuses n
                            var = 0),
                     by = c("treatment", "treatment_description")]),
        use.names = TRUE)[do.call("order", lapply(cols, as.name))]

      # observed number of dropouts before data cut
      dfa3 <- rbindlist(list(
        df0, dt_copy[, list(t = get("totalTime"),
                            n = cumsum(get("dropout")),
                            pilevel = pilevel,
                            lower = NA_real_,
                            upper = NA_real_,
                            mean = cumsum(get("dropout")),  # Reuses n
                            var = 0),
                     by = c("treatment", "treatment_description")]),
        use.names = TRUE)[do.call("order", lapply(cols, as.name))]

      # observed number of ongoing subjects before data cutoff
      df2_copy <- data.table::copy(df2)[, `:=`(tmp_row = .I)]

      dfa4 <- rbindlist(list(
        df0, CJ(t = t2, tmp_row = df2_copy$tmp_row, sorted = FALSE)[
          df2_copy, on = "tmp_row"][
            get("t") <= tp, list(
              n = sum(get("arrivalTime") <= get("t") &
                        (get("totalTime") > get("t") |
                           (!get("event") & !get("dropout")))),
              pilevel = pilevel,
              lower = NA_real_, upper = NA_real_,
              mean = sum(get("arrivalTime") <= get("t") &
                           (get("totalTime") > get("t") |
                              (!get("event") & !get("dropout")))),
              var = 0), by = cols]),
        use.names = TRUE)[do.call("order", lapply(cols, as.name))]


      # add time tp
      dfa2 <- rbindlist(list(
        dfa2, data.table::copy(dfa2)[
          , .SD[.N], keyby = c("treatment", "treatment_description")][
            , `:=`(t = tp)]), use.names = TRUE)[
          , .SD[.N], keyby = cols]

      dfa3 <- rbindlist(list(
        dfa3, data.table::copy(dfa3)[
          , .SD[.N], keyby = c("treatment", "treatment_description")][
            , `:=`(t = tp)]), use.names = TRUE)[
              , .SD[.N], keyby = cols]

      dfa2 <- rbindlist(list(
        dfa2, data.table::copy(dfa2)[
          , .SD[.N], keyby = c("treatment", "treatment_description")][
            , `:=`(t = tp)]), use.names = TRUE)[
              , .SD[.N], keyby = cols]


      # concatenate events before and after data cut
      event_pred_df <- rbindlist(list(dfa2, dfb2), use.names = TRUE)[
        , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
               parameter = "Event")]

      # concatenate dropouts before and after data cut
      dropout_pred_df <- rbindlist(list(dfa3, dfb3), use.names = TRUE)[
        , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
               parameter = "Dropout")]

      # concatenate ongoing subjects before and after data cut
      ongoing_pred_df <- rbindlist(list(dfa4, dfb4), use.names = TRUE)[
        , `:=`(date = as.Date(get("t") - 1, origin = get("trialsdt")),
               parameter = "Ongoing")]
    } else {
      event_pred_df <- dfb2[, `:=`(parameter = "Event")]
      dropout_pred_df <- dfb3[, `:=`(parameter = "Dropout")]
      ongoing_pred_df <- dfb4[, `:=`(parameter = "Ongoing")]
    }

    setorderv(event_pred_df, cols)
    setorderv(dropout_pred_df, cols)
    setorderv(ongoing_pred_df, cols)

    # generate plot
    if (showEnrollment | showEvent | showDropout | showOngoing) {
      dt_list <- list()
      if (showEnrollment) dt_list <- c(dt_list, list(enroll_pred_df))
      if (showEvent) dt_list <- c(dt_list, list(event_pred_df))
      if (showDropout) dt_list <- c(dt_list, list(dropout_pred_df))
      if (showOngoing) dt_list <- c(dt_list, list(ongoing_pred_df))

      dfs <- rbindlist(dt_list, use.names = TRUE)[
        , `:=`(parameter = factor(get("parameter"), levels = c(
          "Enrollment", "Event", "Dropout", "Ongoing")))]

      if (!is.null(df)) { # after trial start
        dfa <- dfs[is.na(get("lower"))]
        dfb <- dfs[!is.na(get("lower"))]

        g1 <- list()
        for (i in c(9999, 1:ngroups)) {
          dfsi <- dfs[get("treatment") == i]
          dfai <- dfa[get("treatment") == i]
          dfbi <- dfb[get("treatment") == i]

          dfai_enrollment <- dfai[get("parameter") == "Enrollment"]
          dfbi_enrollment <- dfbi[get("parameter") == "Enrollment"]
          dfai_event <- dfai[get("parameter") == "Event"]
          dfbi_event <- dfbi[get("parameter") == "Event"]
          dfai_dropout <- dfai[get("parameter") == "Dropout"]
          dfbi_dropout <- dfbi[get("parameter") == "Dropout"]
          dfai_ongoing <- dfai[get("parameter") == "Ongoing"]
          dfbi_ongoing <- dfbi[get("parameter") == "Ongoing"]

          g1[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
            plotly::add_lines(
              data = dfai_enrollment, x = ~date, y = ~n,
              line = list(shape="hv", width=2),
              name = "observed enrollment") %>%
            plotly::add_lines(
              data = dfbi_enrollment, x = ~date, y = ~n,
              line = list(width=2),
              name = "median prediction enrollment") %>%
            plotly::add_ribbons(
              data = dfbi_enrollment, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval enrollment") %>%
            plotly::add_lines(
              data = dfai_event, x = ~date, y = ~n,
              line = list(shape="hv", width=2),
              name = "observed event") %>%
            plotly::add_lines(
              data = dfbi_event, x = ~date, y = ~n,
              line = list(width=2),
              name = "median prediction event") %>%
            plotly::add_ribbons(
              data = dfbi_event, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval event") %>%
            plotly::add_lines(
              data = dfai_dropout, x = ~date, y = ~n,
              line = list(shape="hv", width=2),
              name = "observed dropout") %>%
            plotly::add_lines(
              data = dfbi_dropout, x = ~date, y = ~n,
              line = list(width=2),
              name = "median prediction dropout") %>%
            plotly::add_ribbons(
              data = dfbi_dropout, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval dropout") %>%
            plotly::add_lines(
              data = dfai_ongoing, x = ~date, y = ~n,
              line = list(shape="hv", width=2),
              name = "observed ongoing") %>%
            plotly::add_lines(
              data = dfbi_ongoing, x = ~date, y = ~n,
              line = list(width=2),
              name = "median prediction ongoing") %>%
            plotly::add_ribbons(
              data = dfbi_ongoing, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval ongoing") %>%
            plotly::add_lines(
              x = rep(cutoffdt, 2), y = c(min(dfai$n), max(dfbi$upper)),
              name = "cutoff", line = list(dash="dash"),
              showlegend = FALSE) %>%
            plotly::layout(
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(zeroline = FALSE)) %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))


          if (tp < t0) {
            g1[[(i+1) %% 9999]] <- g1[[(i+1) %% 9999]] %>%
              plotly::add_lines(
                x = rep(cutofftpdt, 2), y = c(min(dfai$n), max(dfbi$upper)),
                name = "prediction start",
                line = list(dash="dash", color="grey"), showlegend = FALSE)
          }


          if (i == 9999) {
            g1[[1]] <- g1[[1]] %>%
              plotly::layout(
                annotations = list(
                  x = cutoffdt, y = 0, text = 'cutoff', xanchor = "left",
                  yanchor = "bottom", font = list(size=12),
                  showarrow = FALSE))

            if (tp < t0) {
              g1[[1]] <- g1[[1]] %>%
                plotly::layout(
                  annotations = list(
                    x = cutofftpdt, y = 0, text = 'prediction start',
                    xanchor = "left", yanchor = "bottom",
                    font = list(size=12), showarrow = FALSE))
            }

            if (showEvent) {
              g1[[1]]  <- g1[[1]] %>%
                plotly::add_lines(
                  x = range(dfsi$date), y = rep(target_d, 2),
                  name = 'target events', showlegend = FALSE,
                  line = list(dash="dot",
                              color="rgba(128, 128, 128, 0.5")) %>%
                plotly::layout(
                  annotations = list(
                    x = 0.95, xref = "paper", y = target_d,
                    text = 'target events', xanchor = "right",
                    yanchor = "bottom", font = list(size=12),
                    showarrow = FALSE))
            }
          }
        }
      } else { # prediction at design stage
        g1 <- list()
        for (i in c(9999, 1:ngroups)) {
          dfsi <- dfs[get("treatment") == i]

          dfsi_enrollment <- dfsi[get("parameter") == "Enrollment"]
          dfsi_event <- dfsi[get("parameter") == "Event"]
          dfsi_dropout <- dfsi[get("parameter") == "Dropout"]
          dfsi_ongoing <- dfsi[get("parameter") == "Ongoing"]

          g1[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
            plotly::add_lines(
              data = dfsi_enrollment, x = ~t, y = ~n,
              line = list(width=2),
              name = "median prediction enrollment") %>%
            plotly::add_ribbons(
              data = dfsi_enrollment, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval enrollment") %>%
            plotly::add_lines(
              data = dfsi_event, x = ~t, y = ~n,
              line = list(width=2),
              name = "median prediction event") %>%
            plotly::add_ribbons(
              data = dfsi_event, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval event") %>%
            plotly::add_lines(
              data = dfsi_dropout, x = ~t, y = ~n,
              line = list(width=2),
              name = "median prediction dropout") %>%
            plotly::add_ribbons(
              data = dfsi_dropout, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval dropout") %>%
            plotly::add_lines(
              data = dfsi_ongoing, x = ~t, y = ~n,
              line = list(width=2),
              name = "median prediction ongoing") %>%
            plotly::add_ribbons(
              data = dfsi_ongoing, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval ongoing") %>%
            plotly::layout(
              xaxis = list(title = "Days since trial start",
                           zeroline = FALSE),
              yaxis = list(zeroline = FALSE)) %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))


          if (i == 9999) {
            if (showEvent) {
              g1[[1]]  <- g1[[1]] %>%
                plotly::add_lines(
                  x = range(dfsi$t), y = rep(target_d, 2),
                  name = 'target events', showlegend = FALSE,
                  line = list(dash="dot",
                              color="rgba(128, 128, 128, 0.5")) %>%
                plotly::layout(
                  annotations = list(
                    x = 0.95, xref = "paper", y = target_d,
                    text = 'target events', xanchor = "right",
                    yanchor = "bottom", font = list(size=12),
                    showarrow = FALSE))
            }
          }
        }
      }
    }
  }


  if (showsummary) cat(s1)
  if (showplot) print(g1)

  if (!is.null(df)) {
    if (showEnrollment | showEvent | showDropout | showOngoing) {
      list(target_d = target_d,
           cutoffdt = cutoffdt, cutofftpdt = cutofftpdt,
           event_pred_day = pred_day, event_pred_date = pred_date,
           pilevel = pilevel, nyears = nyears, nreps = nreps,
           newEvents = newEvents,
           enroll_pred_df = enroll_pred_df,
           event_pred_df = event_pred_df,
           dropout_pred_df = dropout_pred_df,
           ongoing_pred_df = ongoing_pred_df,
           event_pred_summary = s1, event_pred_plot = g1)
    } else {
      list(target_d = target_d,
           cutoffdt = cutoffdt, cutofftpdt = cutofftpdt,
           event_pred_day = pred_day, event_pred_date = pred_date,
           pilevel = pilevel, nyears = nyears, nreps = nreps,
           newEvents = newEvents,
           enroll_pred_df = enroll_pred_df,
           event_pred_df = event_pred_df,
           dropout_pred_df = dropout_pred_df,
           ongoing_pred_df = ongoing_pred_df,
           event_pred_summary = s1)
    }
  } else {
    if (showEnrollment | showEvent | showDropout | showOngoing) {
      list(target_d = target_d,
           event_pred_day = pred_day,
           pilevel = pilevel, nyears = nyears, nreps = nreps,
           newEvents = newEvents,
           enroll_pred_df = enroll_pred_df,
           event_pred_df = event_pred_df,
           dropout_pred_df = dropout_pred_df,
           ongoing_pred_df = ongoing_pred_df,
           event_pred_summary = s1, event_pred_plot = g1)
    } else {
      list(target_d = target_d,
           event_pred_day = pred_day,
           pilevel = pilevel, nyears = nyears, nreps = nreps,
           newEvents = newEvents,
           enroll_pred_df = enroll_pred_df,
           event_pred_df = event_pred_df,
           dropout_pred_df = dropout_pred_df,
           ongoing_pred_df = ongoing_pred_df,
           event_pred_summary = s1)
    }
  }
}
