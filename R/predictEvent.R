#' @title Predict event
#' @description Utilizes pre-fitted time-to-event and time-to-dropout models
#'   to generate event and dropout times for ongoing subjects
#'   and new subjects. It also provides a
#'   prediction interval for the expected time to reach the target
#'   number of events.
#'
#' @param df The subject-level enrollment and event data,
#'   including \code{trialsdt}, \code{randdt}, \code{cutoffdt},
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
#'
#'
#' @details
#' To ensure successful event prediction at the design stage, it is
#' important to provide the \code{newSubjects} data set.
#'
#' To specify the event model used during the design-stage event
#' prediction, the \code{event_fit} be a list with one element
#' per treatment. For each treatment, the element should include \code{w}
#' to specify the weight of the treatment in a randomization block,
#' \code{model} to specify the event model
#' (exponential, weibull, log-normal, or piecewise exponential),
#' \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix.
#' For the piecewise exponential event model, the list
#' should also include \code{piecewiseSurvivalTime} to indicate
#' the location of knots. It should be noted that the model averaging
#' and spline options are not appropriate for use during the design stage.
#'
#' To specify the dropout model used during the design stage
#' event prediction, the \code{dropout_fit} should be a list
#' with one element per treatment. For each treatment, the element
#' should include \code{w} to specify the weight of the treatment
#' in a randomization block, \code{model} to specify the dropout model
#' (exponential, weibull, log-normal, or piecewise exponential),
#' \code{theta} and \code{vtheta} to indicate
#' the parameter values and the covariance matrix.
#' For the piecewise exponential dropout model, the list
#' should also include \code{piecewiseDropoutTime} to indicate
#' the location of knots.
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
#'                            event_fit = event_fit$event_fit,
#'                            dropout_fit = dropout_fit$dropout_fit,
#'                            pilevel = 0.90, nreps = 100)
#'
#' @export
#'
predictEvent <- function(df = NULL, target_d, newSubjects = NULL,
                         event_fit, dropout_fit = NULL,
                         fixedFollowup = FALSE, followupTime = 365,
                         pilevel = 0.90, nyears = 4, nreps = 500,
                         showEnrollment = TRUE, showEvent = TRUE,
                         showDropout = FALSE, showOngoing = FALSE,
                         showsummary = TRUE, showplot = TRUE,
                         by_treatment = FALSE) {

  if (!is.null(df)) erify::check_class(df, "data.frame")
  erify::check_n(target_d)
  if (!is.null(newSubjects)) erify::check_class(newSubjects, "data.frame")
  if (is.null(df) && is.null(newSubjects)) {
    stop("At least one of df and newSubjects must be specified.")
  }

  erify::check_bool(by_treatment)

  if (is.null(df) && "treatment" %in% names(newSubjects))
    by_treatment = TRUE

  # number of treatment groups
  if (by_treatment) {
    if (!is.null(df)) {
      ngroups = length(table(df$treatment))
    } else {
      ngroups = length(table(newSubjects$treatment))
    }

    if (!is.null(df) && !is.null(newSubjects) &&
        length(table(df$treatment)) != length(table(newSubjects$treatment))) {
      stop("Number of treatments must match between df and newSubjects.")
    }

    if (!is.null(df)) {
      if (!("treatment_description" %in% names(df))) {
        df <- df %>% dplyr::mutate(
            treatment_description = paste0("Treatment ", .data$treatment))
      }

      df$treatment_description = stats::reorder(
        as.factor(df$treatment_description), df$treatment)
    }

    if (!is.null(newSubjects)) {
      if (!("treatment_description" %in% names(newSubjects))) {
      newSubjects <- newSubjects %>% dplyr::mutate(
          treatment_description = paste0("Treatment ", .data$treatment))
      }

      newSubjects$treatment_description = stats::reorder(
        as.factor(newSubjects$treatment_description), newSubjects$treatment)
    }

  } else {  # treat as a special case of by-treatment calculation
    ngroups = 1
    if (!is.null(df)) {
      df$treatment = 1
      df$treatment_description = "Overall"
    }
    if (!is.null(newSubjects)) {
      newSubjects$treatment = 1
      newSubjects$treatment_description = "Overall"
    }
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }


  erify::check_class(event_fit, "list")

  if (!by_treatment) {  # convert event_fit to a list with 1 list element
    list1 = event_fit
    event_fit = list()
    event_fit[[1]] = list1
  }

  if (length(event_fit) != ngroups) {
    stop("event_fit must be a list with one element per treatment.")
  }

  # check event_fit model
  if (!is.null(df)) {
    for (j in 1:ngroups) {
      erify::check_content(tolower(event_fit[[j]]$model),
                           c("exponential", "weibull", "log-normal",
                             "piecewise exponential", "model averaging",
                             "spline"))
    }
  } else {
    for (j in 1:ngroups) {
      erify::check_content(tolower(event_fit[[j]]$model),
                           c("exponential", "weibull", "log-normal",
                             "piecewise exponential"))
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
        (model == "log-normal" && p != 2) ||
        (model == "piecewise exponential" &&
         p != length(event_fit[[j]]$piecewiseSurvivalTime))) {
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


  if (!is.null(dropout_fit)) {
    erify::check_class(dropout_fit, "list")

    if (!by_treatment) { # convert dropout_fit to a list with 1 list element
      list1 = dropout_fit
      dropout_fit = list()
      dropout_fit[[1]] = list1
    }

    if (length(dropout_fit) != ngroups) {
      stop("dropout_fit must be a list with one element per treatment.")
    }

    # check dropout_fit model
    for (j in 1:ngroups) {
      erify::check_content(tolower(dropout_fit[[j]]$model),
                           c("exponential", "weibull", "log-normal",
                             "piecewise exponential"))
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
          (model == "log-normal" && p != 2) ||
          (model == "piecewise exponential" &&
           p != length(dropout_fit[[j]]$piecewiseDropoutTime))) {
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

  if (!(showEnrollment | showEvent | showDropout | showOngoing)) {
    stop("At least one parameter must be given for prediction plot.")
  }


  # check input data and extract ongoing subjects
  if (!is.null(df)) {
    df <- dplyr::as_tibble(df)
    names(df) <- tolower(names(df))
    df$trialsdt <- as.Date(df$trialsdt)
    df$randdt <- as.Date(df$randdt)
    df$cutoffdt <- as.Date(df$cutoffdt)

    trialsdt = df$trialsdt[1]
    cutoffdt = df$cutoffdt[1]
    t0 = as.numeric(cutoffdt - trialsdt + 1)

    if (any(df$randdt < trialsdt)) {
      stop("randdt must be greater than or equal to trialsdt.")
    }

    if (any(df$randdt > cutoffdt)) {
      stop("randdt must be less than or equal to cutoffdt.")
    }

    if (any(df$time < 1)) {
      stop("time must be greater than or equal to 1.")
    }

    if (any(df$event == 1 & df$dropout == 1)) {
      stop("event and dropout cannot both be equal to 1 simultaneously.")
    }

    if (any(df$time > as.numeric(cutoffdt - df$randdt + 1))) {
      stop("time must be less than or equal to cutoffdt - randdt + 1.")
    }

    df <- df %>%
      dplyr::mutate(arrivalTime = as.numeric(.data$randdt - trialsdt + 1),
                    totalTime = .data$arrivalTime + .data$time - 1)

    # subset to extract ongoing subjects
    ongoingSubjects <- df %>%
      dplyr::filter(.data$event == 0 & .data$dropout == 0)

    arrivalTimeOngoing <- ongoingSubjects$arrivalTime
    treatmentOngoing <- ongoingSubjects$treatment
    treatment_descriptionOngoing <- ongoingSubjects$treatment_description
    time0Ongoing <- ongoingSubjects$time
    tp = min(ongoingSubjects$totalTime) # tp <= t0
    cutofftpdt <- as.Date(tp - 1, origin = trialsdt)
    n0 = nrow(df)
    d0 = sum(df$event)
    c0 = sum(df$dropout)
    r0 = nrow(ongoingSubjects)

    # subjects who have had the event or dropped out
    stoppedSubjects <- df %>%
      dplyr::filter(.data$event == 1 | .data$dropout == 1)
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
    new1 <- newSubjects %>%
      dplyr::group_by(.data$draw) %>%
      dplyr::slice(dplyr::n())

    pred_day1 <- ceiling(quantile(new1$arrivalTime, c(0.5, plower, pupper)))

    # future time points at which to predict number of subjects
    t = sort(unique(c(seq(t0, t1, 30), t1, pred_day1)))
    t = t[t <= t1]
  }


  # enrollment prediction data
  if (!by_treatment) {
    if (!is.null(newSubjects)) {
      # predicted number of subjects enrolled after data cut
      dfb1 <- dplyr::tibble(t = t) %>%
        dplyr::cross_join(newSubjects) %>%
        dplyr::group_by(.data$t, .data$draw) %>%
        dplyr::summarise(nenrolled = sum(.data$arrivalTime <= .data$t) + n0,
                         .groups = "drop_last") %>%
        dplyr::summarise(n = quantile(.data$nenrolled, probs = 0.5),
                         lower = quantile(.data$nenrolled, probs = plower),
                         upper = quantile(.data$nenrolled, probs = pupper),
                         mean = mean(.data$nenrolled),
                         var = var(.data$nenrolled)) %>%
        dplyr::ungroup()
    }


    if (!is.null(df)) {
      # day 1
      df0 <- dplyr::tibble(t = 1, n = 0, lower = NA, upper = NA,
                           mean = 0, var = 0)

      # arrival time for subjects already enrolled before data cut
      dfa1 <- df %>%
        dplyr::arrange(.data$randdt) %>%
        dplyr::mutate(t = as.numeric(.data$randdt - trialsdt + 1),
                      n = dplyr::row_number()) %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::bind_rows(dplyr::tibble(t = t0, n = n0, lower = NA,
                                       upper = NA, mean = n0, var= 0)) %>%
        dplyr::select(.data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::group_by(.data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()
    }


    if (is.null(newSubjects)) { # existing subjects only
      # add predicted from data cut to specified years after data cut
      dfb1t0 <- dfa1 %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(lower = .data$n, upper = .data$n)
      dfb1t1 <- dfb1t0 %>%
        dplyr::mutate(t = t1)

      enroll_pred_df <- dfa1 %>%
        dplyr::bind_rows(dfb1t0) %>%
        dplyr::bind_rows(dfb1t1) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Enrollment")
    } else if (is.null(df)) { # new subjects only
      enroll_pred_df <- dfb1 %>%
        dplyr::mutate(parameter = "Enrollment")
    } else { # existing and new subjects
      enroll_pred_df <- dfa1 %>%
        dplyr::bind_rows(dfb1) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Enrollment")
    }

    enroll_pred_df <- enroll_pred_df %>%
      dplyr::arrange(.data$t)
  } else { # by treatment

    # summary of observed data by treatment
    if (!is.null(df)) {
      # add overall treatment
      df2 <- df %>%
        dplyr::bind_rows(df %>% dplyr::mutate(
          treatment = 9999, treatment_description = "Overall"))

      sum_by_trt <- df2 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::summarise(n0 = dplyr::n(),
                         d0 = sum(.data$event),
                         c0 = sum(.data$dropout),
                         r0 = sum(!(.data$event | .data$dropout)),
                         .groups = "drop")
    }

    if (!is.null(newSubjects)) {
      # add overall treatment
      newSubjects2 <- newSubjects %>%
        dplyr::bind_rows(newSubjects %>% dplyr::mutate(
          treatment = 9999, treatment_description = "Overall"))

      if (is.null(df)) {
        sum_by_trt <- newSubjects2 %>%
          dplyr::group_by(.data$treatment, .data$treatment_description) %>%
          dplyr::slice(dplyr::n()) %>%
          dplyr::mutate(n0 = 0, d0 = 0, c0 = 0, r0 = 0) %>%
          dplyr::select(.data$treatment, .data$treatment_description,
                        .data$n0, .data$d0, .data$c0, .data$r0)
      }

      # predicted number of subjects enrolled by treatment after cutoff
      dfb1 <- dplyr::tibble(t = t) %>%
        dplyr::cross_join(newSubjects2) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$t, .data$draw) %>%
        dplyr::summarise(nenrolled = sum(.data$arrivalTime <= .data$t),
                         .groups = "drop_last") %>%
        dplyr::left_join(sum_by_trt,
                         by = c("treatment", "treatment_description")) %>%
        dplyr::mutate(nenrolled = .data$nenrolled + .data$n0) %>%
        dplyr::summarise(n = quantile(.data$nenrolled, probs = 0.5),
                         lower = quantile(.data$nenrolled, probs = plower),
                         upper = quantile(.data$nenrolled, probs = pupper),
                         mean = mean(.data$nenrolled),
                         var = var(.data$nenrolled),
                         .groups = "drop_last") %>%
        dplyr::ungroup()
    }


    if (!is.null(df)) {
      # day 1
      df0 <- sum_by_trt %>%
        dplyr::select(.data$treatment, .data$treatment_description) %>%
        dplyr::mutate(t = 1, n = 0, lower = NA, upper = NA, mean = 0, var = 0)

      # arrival time for subjects already enrolled before data cut
      dfa1 <- df2 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::arrange(.data$randdt) %>%
        dplyr::mutate(t = as.numeric(.data$randdt - trialsdt + 1),
                      n = dplyr::row_number()) %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::bind_rows(sum_by_trt %>%
                           dplyr::mutate(t = t0, n = n0, lower = NA,
                                         upper = NA, mean = n0, var = 0)) %>%
        dplyr::select(.data$treatment, .data$treatment_description,
                      .data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description)
    }


    if (is.null(newSubjects)) { # existing subjects only
      # add predicted from data cut to specified years after data cut
      dfb1t0 <- dfa1 %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(lower = .data$n, upper = .data$n)
      dfb1t1 <- dfb1t0 %>%
        dplyr::mutate(t = t1)

      enroll_pred_df <- dfa1 %>%
        dplyr::bind_rows(dfb1t0) %>%
        dplyr::bind_rows(dfb1t1) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Enrollment")
    } else if (is.null(df)) { # new subjects only
      enroll_pred_df <- dfb1 %>%
        dplyr::mutate(parameter = "Enrollment")
    } else { # existing and new subjects
      enroll_pred_df <- dfa1 %>%
        dplyr::bind_rows(dfb1) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Enrollment")
    }

    enroll_pred_df <- enroll_pred_df %>%
      dplyr::arrange(.data$treatment, .data$treatment_description, .data$t)
  }


  # extract posterior draws of model parameters
  theta2 <- list()
  for (j in 1:ngroups) {
    if (length(event_fit[[j]]$theta) == 1) {
      theta2[[j]] <- matrix(rnorm(nreps, mean = event_fit[[j]]$theta,
                                  sd = sqrt(event_fit[[j]]$vtheta)),
                            ncol=1)
    } else {
      theta2[[j]] = mvtnorm::rmvnorm(nreps, mean = event_fit[[j]]$theta,
                                     sigma = event_fit[[j]]$vtheta)
    }
  }


  if (!is.null(dropout_fit)) {
    theta3 <- list()
    for (j in 1:ngroups) {
      if (length(dropout_fit[[j]]$theta) == 1) {
        theta3[[j]] <- matrix(rnorm(nreps, mean = dropout_fit[[j]]$theta,
                                    sd = sqrt(dropout_fit[[j]]$vtheta)),
                              ncol=1)
      } else {
        theta3[[j]] = mvtnorm::rmvnorm(nreps, mean = dropout_fit[[j]]$theta,
                                       sigma = dropout_fit[[j]]$vtheta)
      }
    }
  }


  # generate the event and dropout times
  if (!is.null(newSubjects)) {
    m1 = nrow(newSubjects)
  } else {
    m1 = 0
  }

  newEvents = dplyr::as_tibble(matrix(
    nrow = nreps*r0 + m1, ncol = 7,
    dimnames = list(NULL, c("draw", "arrivalTime",
                            "treatment", "treatment_description",
                            "time", "event", "dropout"))))

  offset = 0
  for (i in 1:nreps) {
    # number of new subjects in the simulated data set
    if (m1 > 0) {
      n1 = nrow(newSubjects %>% dplyr::filter(.data$draw == i))
    } else {
      n1 = 0
    }

    m = r0 + n1

    # arrival time, treatment, and time offset for new subjects
    if (n1 > 0) {
      arrivalTimeNew = newSubjects$arrivalTime[newSubjects$draw == i]
      treatmentNew = newSubjects$treatment[newSubjects$draw == i]
      treatment_descriptionNew = newSubjects$treatment_description[
        newSubjects$draw == i]
      time0New = rep(0, n1)
    }

    # concatenate ongoing and new subjects
    if (r0 == 0 && n1 > 0) {  # design stage
      arrivalTime = arrivalTimeNew
      treatment = treatmentNew
      treatment_description = treatment_descriptionNew
      time0 = time0New
    } else if (r0 > 0 && n1 > 0) { # enrollment stage
      arrivalTime = c(arrivalTimeOngoing, arrivalTimeNew)
      treatment = c(treatmentOngoing, treatmentNew)
      treatment_description = c(treatment_descriptionOngoing,
                                treatment_descriptionNew)
      time0 = c(time0Ongoing, time0New)
    } else if (r0 > 0 && n1 == 0) { # follow-up stage
      arrivalTime = arrivalTimeOngoing
      treatment = treatmentOngoing
      treatment_description = treatment_descriptionOngoing
      time0 = time0Ongoing
    }


    # draw event time for ongoing and new subjects
    survivalTime = rep(NA, m)

    for (j in 1:ngroups) {
      cols = which(treatment == j)
      ncols = length(cols)

      if (ncols > 0) {
        model = tolower(event_fit[[j]]$model)

        if (model == "exponential") {
          rate = exp(theta2[[j]][i,])
          survivalTime[cols] = rexp(ncols, rate) + time0[cols]
        } else if (model == "weibull") {
          shape = exp(theta2[[j]][i,1])
          scale = exp(theta2[[j]][i,2])
          survivalTime[cols] = (rexp(ncols)*scale^shape +
                                  time0[cols]^shape)^(1/shape)
        } else if (model == "log-normal") {
          meanlog = theta2[[j]][i,1]
          sdlog = exp(theta2[[j]][i,2])
          survivalTime[cols] = exp(tmvtnsim::rtnorm(
            mean = rep(meanlog, ncols), sd = sdlog,
            lower = log(time0[cols]), upper = rep(Inf, ncols)))
        } else if (model == "piecewise exponential") {
          lambda = exp(theta2[[j]][i,]) # hazard rates in the intervals
          J = length(lambda) # number of intervals
          u = event_fit[[j]]$piecewiseSurvivalTime # left end points

          # partial sums of lambda*interval_width
          if (J > 1) {
            psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))
          } else {
            psum = 0
          }

          # find the interval containing time0
          j0 = findInterval(time0[cols], u)
          rhs = psum[j0] + lambda[j0]*(time0[cols] - u[j0]) + rexp(ncols)

          # find the interval containing time
          j1 = findInterval(rhs, psum)
          survivalTime[cols] = u[j1] + (rhs - psum[j1])/lambda[j1]
        } else if (model == "model averaging") {
          theta = theta2[[j]][i,]
          shape = exp(theta[1])
          scale = exp(theta[2])
          meanlog = theta[3]
          sdlog = exp(theta[4])
          w1 = event_fit[[j]]$w1

          # draw component indicator
          p1 = w1*pweibull(time0[cols], shape, scale, lower.tail = FALSE)
          p2 = (1-w1)*plnorm(time0[cols], meanlog, sdlog, lower.tail = FALSE)
          w = (runif(ncols) < p1/(p1+p2))
          nw1 = sum(w)
          nw0 = ncols - nw1

          # draw from the corresponding component distribution
          if (nw1 > 0) {
            survivalTime[cols][w==1] = (rexp(nw1)*scale^shape +
                                          time0[cols][w==1]^shape)^(1/shape)
          }

          if (nw0 > 0) {
            survivalTime[cols][w==0] = exp(tmvtnsim::rtnorm(
              mean = rep(meanlog, nw0), sd = sdlog,
              lower = log(time0[cols][w==0]), upper = rep(Inf, nw0)))
          }
        } else if (model == "spline") {
          gamma = theta2[[j]][i,]
          knots = event_fit[[j]]$knots
          scale = event_fit[[j]]$scale

          st0 = flexsurv::psurvspline(
            time0[cols], gamma, knots = knots, scale = scale,
            lower.tail = FALSE)

          survivalTime[cols] = flexsurv::qsurvspline(
            runif(ncols)*st0, gamma, knots = knots, scale = scale,
            lower.tail = FALSE)
        }
      }
    }

    # new subjects start with day 1 on arrival
    if (n1 > 0) survivalTime[(r0+1):m] = survivalTime[(r0+1):m] + 1


    # draw dropout time for ongoing and new subjects
    if (!is.null(dropout_fit)) {
      dropoutTime = rep(NA, m)
      for (j in 1:ngroups) {
        cols = which(treatment == j)
        ncols = length(cols)

        if (ncols > 0) {
          model = tolower(dropout_fit[[j]]$model)

          if (model == "exponential") {
            rate = exp(theta3[[j]][i,])
            dropoutTime[cols] = rexp(ncols, rate) + time0[cols]
          } else if (model == "weibull") {
            shape = exp(theta3[[j]][i,1])
            scale = exp(theta3[[j]][i,2])
            dropoutTime[cols] = (rexp(ncols)*scale^shape +
                                   time0[cols]^shape)^(1/shape)
          } else if (model == "log-normal") {
            meanlog = theta3[[j]][i,1]
            sdlog = exp(theta3[[j]][i,2])
            dropoutTime[cols] = exp(tmvtnsim::rtnorm(
              mean = rep(meanlog, ncols), sd = sdlog,
              lower = log(time0[cols]), upper = rep(Inf, ncols)))
          } else if (model == "piecewise exponential") {
            lambda = exp(theta3[[j]][i,]) # hazard rates in the intervals
            J = length(lambda) # number of intervals
            u = dropout_fit[[j]]$piecewiseDropoutTime # left end points

            # partial sums of lambda*interval_width
            if (J > 1) {
              psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))
            } else {
              psum = 0
            }

            # find the interval containing time0
            j0 = findInterval(time0[cols], u)
            rhs = psum[j0] + lambda[j0]*(time0[cols] - u[j0]) + rexp(ncols)

            # find the interval containing time
            j1 = findInterval(rhs, psum)
            dropoutTime[cols] = u[j1] + (rhs - psum[j1])/lambda[j1]
          }
        }
      }

      # new subjects start with day 1 on arrival
      if (n1 > 0) dropoutTime[(r0+1):m] = dropoutTime[(r0+1):m] + 1
    }


    # observed survival time and event indicator
    if (!fixedFollowup) {
      if (!is.null(dropout_fit)) {
        time = pmin(survivalTime, dropoutTime)
        event = 1*(time == survivalTime)
        dropout = 1*(time == dropoutTime)
      } else {
        time = survivalTime
        event = 1
        dropout = 0
      }
    } else {
      if (!is.null(dropout_fit)) {
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
    newEvents[index, "arrivalTime"] = arrivalTime
    newEvents[index, "treatment"] = treatment
    newEvents[index, "treatment_description"] = treatment_description
    newEvents[index, "time"] = time
    newEvents[index, "event"] = event
    newEvents[index, "dropout"] = dropout
    offset = offset + m
  }


  # calculate total time since trial start
  newEvents <- newEvents %>%
    dplyr::mutate(totalTime = .data$arrivalTime + .data$time - 1)

  if (!is.null(df)) {
    # combined stopped, ongoing and new subjects
    allSubjects <- dplyr::tibble(draw = 1:nreps) %>%
      dplyr::cross_join(stoppedSubjects) %>%
      dplyr::select(.data$draw, .data$arrivalTime,
                    .data$treatment, .data$treatment_description,
                    .data$time, .data$event, .data$dropout,
                    .data$totalTime) %>%
      dplyr::bind_rows(newEvents)
  } else {
    allSubjects <- newEvents
  }

  # remove the dummy treatment from newEvents
  if (!by_treatment) newEvents <- newEvents %>%
    dplyr::select(-c(treatment, treatment_description))

  # A general quantile method if there are data sets not reaching target_d
  # Find t such that sum(I{D_i(t) < target_d}, {i, 1, nreps}) / nreps = q.
  # This works because {D_i(t) < target_d} = {T_i(target_d) > t},
  # where D_i(t) is the cumulative number of events at time t, and
  # T_i(target_d) is the time to reach target_d for data set i.
  sdf <- function(t, target_d, d0, newEvents) {
    sumdata <- newEvents %>%
      dplyr::group_by(.data$draw) %>%
      dplyr::summarize(n = sum(.data$totalTime <= t & .data$event == 1) + d0)
    mean(sumdata$n < target_d)
  }

  # obtain the quantiles
  q = 1 - c(0.5, plower, pupper)
  pred_day = rep(NA, length(q))
  tmax = max(newEvents$totalTime[newEvents$event==1])
  for (j in 1:length(q)) {
    # check if the quantile can be estimated from observed data
    if (sdf(tmax, target_d, d0, newEvents) <= q[j]) {
      pred_day[j] = uniroot(function(x)
        sdf(x, target_d, d0, newEvents) - q[j],
        c(tp, tmax), tol = 1)$root
      pred_day[j] = ceiling(pred_day[j])
    }
  }
  names(pred_day) <- names(quantile(1:100, c(0.5, plower, pupper)))

  if (!is.null(df)) {
    pred_date <- as.Date(pred_day - 1, origin = trialsdt)

    str1 <- paste0("Time from cutoff until ", target_d, " events: ",
                   pred_date[1] - cutoffdt + 1, " days")
    str2 <- paste0("Median prediction date: ", pred_date[1])
    str3 <- paste0("Prediction interval: ", pred_date[2], ", ", pred_date[3])
    s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")
  } else {
    str1 <- paste0("Time from trial start until ", target_d, " events")
    str2 <- paste0("Median prediction day: ", pred_day[1])
    str3 <- paste0("Prediction interval: ", pred_day[2], ", ", pred_day[3])
    s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")
  }


  # observed time points
  t2 = sort(unique(c(df$arrivalTime, df$totalTime)))

  # future time points at which to predict number of events
  t = unique(c(t2[t2 >= tp], seq(t0, t1, 30), t1))

  if (!by_treatment) {
    # number of events, dropouts, and ongoing subjects after data cut
    df1 = dplyr::tibble(t = t) %>%
      dplyr::cross_join(allSubjects) %>%
      dplyr::group_by(.data$t, .data$draw) %>%
      dplyr::summarise(nevents = sum(.data$totalTime <= .data$t &
                                       .data$event == 1),
                       ndropouts = sum(.data$totalTime <= .data$t &
                                         .data$dropout == 1),
                       nongoing = sum(.data$arrivalTime <= .data$t &
                                        .data$totalTime > .data$t),
                       .groups = "drop_last")

    # predicted number of events after data cut
    dfb2 = df1 %>%
      dplyr::summarise(n = quantile(.data$nevents, probs = 0.5),
                       lower = quantile(.data$nevents, probs = plower),
                       upper = quantile(.data$nevents, probs = pupper),
                       mean = mean(.data$nevents),
                       var = var(.data$nevents))

    # predicted number of dropouts after data cut
    dfb3 = df1 %>%
      dplyr::summarise(n = quantile(.data$ndropouts, probs = 0.5),
                       lower = quantile(.data$ndropouts, probs = plower),
                       upper = quantile(.data$ndropouts, probs = pupper),
                       mean = mean(.data$ndropouts),
                       var = var(.data$ndropouts))

    # predicted number of subjects at risk after data cut
    dfb4 = df1 %>%
      dplyr::summarise(n = quantile(.data$nongoing, probs = 0.5),
                       lower = quantile(.data$nongoing, probs = plower),
                       upper = quantile(.data$nongoing, probs = pupper),
                       mean = mean(.data$nongoing),
                       var = var(.data$nongoing))

    if (!is.null(df)) {
      # day 1
      df0 <- dplyr::tibble(t = 1, n = 0, lower = NA, upper = NA,
                           mean = 0, var = 0)

      # observed number of events before data cut
      dfa2 <- df %>%
        dplyr::arrange(.data$totalTime) %>%
        dplyr::mutate(t = .data$totalTime, n = cumsum(.data$event)) %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::select(.data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::filter(.data$t <= tp) %>%
        dplyr::arrange(.data$t)

      # observed number of dropouts before data cut
      dfa3 <- df %>%
        dplyr::arrange(.data$totalTime) %>%
        dplyr::mutate(t = .data$totalTime, n = cumsum(.data$dropout)) %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::select(.data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::filter(.data$t <= tp) %>%
        dplyr::arrange(.data$t)


      # observed number of ongoing subjects before data cutoff
      dfa4 <- dplyr::tibble(t = t2) %>%
        dplyr::cross_join(df) %>%
        dplyr::group_by(.data$t) %>%
        dplyr::summarise(n = sum(.data$arrivalTime <= .data$t &
                                   (.data$totalTime > .data$t |
                                      (.data$event == 0 &
                                         .data$dropout == 0))),
                         .groups = "drop_last") %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::select(.data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::filter(.data$t <= tp) %>%
        dplyr::arrange(.data$t)


      # add time tp
      dfa2tp <- dfa2 %>%
        dplyr::ungroup() %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(t = tp)

      dfa2 <- dfa2 %>%
        dplyr::bind_rows(dfa2tp) %>%
        dplyr::group_by(.data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()


      dfa3tp <- dfa3 %>%
        dplyr::ungroup() %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(t = tp)

      dfa3 <- dfa3 %>%
        dplyr::bind_rows(dfa3tp) %>%
        dplyr::group_by(.data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()


      dfa4tp <- dfa4 %>%
        dplyr::ungroup() %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(t = tp)

      dfa4 <- dfa4 %>%
        dplyr::bind_rows(dfa4tp) %>%
        dplyr::group_by(.data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()


      # concatenate events before and after data cut
      event_pred_df <- dfa2 %>%
        dplyr::bind_rows(dfb2) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Event")

      # concatenate dropouts before and after data cut
      dropout_pred_df <- dfa3 %>%
        dplyr::bind_rows(dfb3) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Dropout")

      # concatenate ongoing subjects before and after data cut
      ongoing_pred_df <- dfa4 %>%
        dplyr::bind_rows(dfb4) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Ongoing")
    } else {
      event_pred_df <- dfb2 %>%
        dplyr::mutate(parameter = "Event")

      dropout_pred_df <- dfb3 %>%
        dplyr::mutate(parameter = "Dropout")

      ongoing_pred_df <- dfb4 %>%
        dplyr::mutate(parameter = "Ongoing")
    }

    event_pred_df <- event_pred_df %>%
      dplyr::arrange(.data$t)

    dropout_pred_df <- dropout_pred_df %>%
      dplyr::arrange(.data$t)

    ongoing_pred_df <- ongoing_pred_df %>%
      dplyr::arrange(.data$t)


    dfs <- dplyr::tibble()
    if (showEnrollment) dfs <- dfs %>% dplyr::bind_rows(enroll_pred_df)
    if (showEvent) dfs <- dfs %>% dplyr::bind_rows(event_pred_df)
    if (showDropout) dfs <- dfs %>% dplyr::bind_rows(dropout_pred_df)
    if (showOngoing) dfs <- dfs %>% dplyr::bind_rows(ongoing_pred_df)

    dfs$parameter <- factor(dfs$parameter, levels = c(
      "Enrollment", "Event", "Dropout", "Ongoing"))

    if (!is.null(df)) {
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      g1 <- plotly::plot_ly() %>%
        plotly::add_ribbons(
          data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", fillcolor = ~parameter,
          line = list(width=0)) %>%
        plotly::add_lines(
          data = dfb, x = ~date, y = ~n, color = ~parameter,
          line = list(width=2)) %>%
        plotly::add_lines(
          data = dfa, x = ~date, y = ~n, color = ~parameter,
          line = list(shape="hv", width=2)) %>%
        plotly::add_lines(
          x = rep(cutoffdt, 2), y = range(dfs$n), name = "cutoff",
          line = list(dash="dash"), showlegend = FALSE) %>%
        plotly::layout(
          annotations = list(
            x = cutoffdt, y = 0, text = 'cutoff', xanchor = "left",
            yanchor = "bottom", font = list(size=12), showarrow = FALSE),
          xaxis = list(title = "", zeroline = FALSE),
          yaxis = list(zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = 'h'))

      if (tp < t0) {
        g1 <- g1 %>%
          plotly::add_lines(
            x = rep(cutofftpdt, 2), y = range(dfs$n),
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
              text = 'target events', xanchor = "right", yanchor = "bottom",
              font = list(size=12), showarrow = FALSE))
      }

    } else {
      g1 <- plotly::plot_ly() %>%
        plotly::add_ribbons(
          data = dfs, x = ~t, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", fillcolor = ~parameter,
          line = list(width=0)) %>%
        plotly::add_lines(
          data = dfs, x = ~t, y = ~n, color = ~parameter,
          line = list(width=2)) %>%
        plotly::layout(
          xaxis = list(title = "Days since trial start", zeroline = FALSE),
          yaxis = list(zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = 'h'))

      if (showEvent) {
        g1 <- g1 %>%
          plotly::add_lines(
            x = range(dfs$t), y = rep(target_d, 2),
            name = 'target events', showlegend = FALSE,
            line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
          plotly::layout(
            annotations = list(
              x = 0.95, xref = "paper", y = target_d,
              text = 'target events', xanchor = "right", yanchor = "bottom",
              font = list(size=12), showarrow = FALSE))
      }
    }
  } else {  # by treatment
    # add overall treatment
    allSubjects2 <- allSubjects %>%
      dplyr::bind_rows(allSubjects %>% dplyr::mutate(
        treatment = 9999, treatment_description = "Overall"))

    # number of events, dropouts, and ongoing subjects after data cut
    df1 = dplyr::tibble(t = t) %>%
      dplyr::cross_join(allSubjects2) %>%
      dplyr::group_by(.data$treatment, .data$treatment_description,
                      .data$t, .data$draw) %>%
      dplyr::summarise(nevents = sum(.data$totalTime <= .data$t &
                                       .data$event == 1),
                       ndropouts = sum(.data$totalTime <= .data$t &
                                         .data$dropout == 1),
                       nongoing = sum(.data$arrivalTime <= .data$t &
                                        .data$totalTime > .data$t),
                       .groups = "drop_last")

    # predicted number of events after data cut
    dfb2 = df1 %>%
      dplyr::summarise(n = quantile(.data$nevents, probs = 0.5),
                       lower = quantile(.data$nevents, probs = plower),
                       upper = quantile(.data$nevents, probs = pupper),
                       mean = mean(.data$nevents),
                       var = var(.data$nevents),
                       .groups = "drop_last")

    # predicted number of dropouts after data cut
    dfb3 = df1 %>%
      dplyr::summarise(n = quantile(.data$ndropouts, probs = 0.5),
                       lower = quantile(.data$ndropouts, probs = plower),
                       upper = quantile(.data$ndropouts, probs = pupper),
                       mean = mean(.data$ndropouts),
                       var = var(.data$ndropouts),
                       .groups = "drop_last")

    # predicted number of subjects at risk after data cut
    dfb4 = df1 %>%
      dplyr::summarise(n = quantile(.data$nongoing, probs = 0.5),
                       lower = quantile(.data$nongoing, probs = plower),
                       upper = quantile(.data$nongoing, probs = pupper),
                       mean = mean(.data$nongoing),
                       var = var(.data$nongoing),
                       .groups = "drop_last")


    if (!is.null(df)) {
      # day 1
      df0 <- sum_by_trt %>%
        dplyr::select(.data$treatment, .data$treatment_description) %>%
        dplyr::mutate(t = 1, n = 0, lower = NA, upper = NA, mean = 0, var = 0)


      # observed number of events before data cut
      dfa2 <- df2 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::arrange(.data$totalTime) %>%
        dplyr::mutate(t = .data$totalTime, n = cumsum(.data$event)) %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::select(.data$treatment, .data$treatment_description,
                      .data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::filter(.data$t <= tp) %>%
        dplyr::arrange(.data$treatment, .data$treatment_description, .data$t)

      # observed number of dropouts before data cut
      dfa3 <- df2 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::arrange(.data$totalTime) %>%
        dplyr::mutate(t = .data$totalTime, n = cumsum(.data$dropout)) %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::select(.data$treatment, .data$treatment_description,
                      .data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::filter(.data$t <= tp) %>%
        dplyr::arrange(.data$treatment, .data$treatment_description, .data$t)


      # observed number of ongoing subjects before data cutoff
      dfa4 <- dplyr::tibble(t = t2) %>%
        dplyr::cross_join(df2) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$t) %>%
        dplyr::summarise(
          n = sum(.data$arrivalTime <= .data$t &
                    (.data$totalTime > .data$t |
                       (.data$event == 0 & .data$dropout == 0))),
          .groups = "drop_last") %>%
        dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
        dplyr::select(.data$treatment, .data$treatment_description,
                      .data$t, .data$n, .data$lower, .data$upper,
                      .data$mean, .data$var) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::filter(.data$t <= tp) %>%
        dplyr::arrange(.data$treatment, .data$treatment_description, .data$t)


      # add time tp
      dfa2tp <- dfa2 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(t = tp)

      dfa2 <- dfa2 %>%
        dplyr::bind_rows(dfa2tp) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()


      dfa3tp <- dfa3 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(t = tp)

      dfa3 <- dfa3 %>%
        dplyr::bind_rows(dfa3tp) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()


      dfa4tp <- dfa4 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(t = tp)

      dfa4 <- dfa4 %>%
        dplyr::bind_rows(dfa4tp) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()


      # concatenate events before and after data cut
      event_pred_df <- dfa2 %>%
        dplyr::bind_rows(dfb2) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Event")

      # concatenate dropouts before and after data cut
      dropout_pred_df <- dfa3 %>%
        dplyr::bind_rows(dfb3) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Dropout")

      # concatenate ongoing subjects before and after data cut
      ongoing_pred_df <- dfa4 %>%
        dplyr::bind_rows(dfb4) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Ongoing")
    } else {
      event_pred_df <- dfb2 %>%
        dplyr::mutate(parameter = "Event")

      dropout_pred_df <- dfb3 %>%
        dplyr::mutate(parameter = "Dropout")

      ongoing_pred_df <- dfb4 %>%
        dplyr::mutate(parameter = "Ongoing")
    }

    event_pred_df <- event_pred_df %>%
      dplyr::arrange(.data$treatment, .data$treatment_description, .data$t)

    dropout_pred_df <- dropout_pred_df %>%
      dplyr::arrange(.data$treatment, .data$treatment_description, .data$t)

    ongoing_pred_df <- ongoing_pred_df %>%
      dplyr::arrange(.data$treatment, .data$treatment_description, .data$t)

    dfs <- dplyr::tibble()
    if (showEnrollment) dfs <- dfs %>% dplyr::bind_rows(enroll_pred_df)
    if (showEvent) dfs <- dfs %>% dplyr::bind_rows(event_pred_df)
    if (showDropout) dfs <- dfs %>% dplyr::bind_rows(dropout_pred_df)
    if (showOngoing) dfs <- dfs %>% dplyr::bind_rows(ongoing_pred_df)

    dfs$parameter <- factor(dfs$parameter, levels = c(
      "Enrollment", "Event", "Dropout", "Ongoing"))

    if (!is.null(df)) {
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      g <- list()
      for (i in c(9999, 1:ngroups)) {
        dfsi <- dfs %>%
          dplyr::filter(.data$treatment == i)
        dfbi <- dfb %>%
          dplyr::filter(.data$treatment == i)
        dfai <- dfa %>%
          dplyr::filter(.data$treatment == i)

        g[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
          plotly::add_ribbons(
            data = dfbi, x = ~date, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", fillcolor = ~parameter,
            line = list(width=0)) %>%
          plotly::add_lines(
            data = dfbi, x = ~date, y = ~n, color = ~parameter,
            line = list(width=2)) %>%
          plotly::add_lines(
            data = dfai, x = ~date, y = ~n, color = ~parameter,
            line = list(shape="hv", width=2)) %>%
          plotly::add_lines(
            x = rep(cutoffdt, 2), y = range(dfsi$n), name = "cutoff",
            line = list(dash="dash"), showlegend = FALSE) %>%
          plotly::layout(
            xaxis = list(title = "", zeroline = FALSE),
            yaxis = list(zeroline = FALSE),
            legend = list(x = 0, y = 1.05, yanchor = "bottom",
                          orientation = 'h')) %>%
          plotly::layout(
            annotations = list(
              x = 0.5, y = 1,
              text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
              xanchor = "center", yanchor = "bottom",
              showarrow = FALSE, xref='paper', yref='paper'))


        if (tp < t0) {
          g[[(i+1) %% 9999]] <- g[[(i+1) %% 9999]] %>%
            plotly::add_lines(
              x = rep(cutofftpdt, 2), y = range(dfsi$n),
              name = "prediction start",
              line = list(dash="dash", color="grey"), showlegend = FALSE)
        }


        if (i == 9999) {
          g[[1]] <- g[[1]] %>%
            plotly::layout(
              annotations = list(
                x = cutoffdt, y = 0, text = 'cutoff', xanchor = "left",
                yanchor = "bottom", font = list(size=12),
                showarrow = FALSE))

          if (tp < t0) {
            g[[1]] <- g[[1]] %>%
              plotly::layout(
                annotations = list(
                  x = cutofftpdt, y = 0, text = 'prediction start',
                  xanchor = "left", yanchor = "bottom",
                  font = list(size=12), showarrow = FALSE))
          }

          if (showEvent) {
            g[[1]]  <- g[[1]] %>%
              plotly::add_lines(
                x = range(dfsi$date), y = rep(target_d, 2),
                name = 'target events', showlegend = FALSE,
                line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
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
      g <- list()
      for (i in c(9999, 1:ngroups)) {
        dfsi <- dfs %>%
          dplyr::filter(.data$treatment == i)

        g[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
          plotly::add_ribbons(
            data = dfsi, x = ~t, ymin = ~lower, ymax = ~upper,
            fill = "tonexty", fillcolor = ~parameter,
            line = list(width=0)) %>%
          plotly::add_lines(
            data = dfsi, x = ~t, y = ~n, color = ~parameter,
            line = list(width=2)) %>%
          plotly::layout(
            xaxis = list(title = "Days since trial start", zeroline = FALSE),
            yaxis = list(zeroline = FALSE),
            legend = list(x = 0, y = 1.05, yanchor = "bottom",
                          orientation = 'h')) %>%
          plotly::layout(
            annotations = list(
              x = 0.5, y = 1,
              text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
              xanchor = "center", yanchor = "bottom",
              showarrow = FALSE, xref='paper', yref='paper'))


        if (i == 9999) {
          if (showEvent) {
            g[[1]]  <- g[[1]] %>%
              plotly::add_lines(
                x = range(dfsi$t), y = rep(target_d, 2),
                name = 'target events', showlegend = FALSE,
                line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
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

    g1 <- plotly::subplot(g, nrows = ngroups + 1, margin = 0.05)
  }


  if (showsummary) cat(s1)
  if (showplot) print(g1)

  if (!is.null(df)) {
    list(target_d = target_d,
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
         event_pred_day = pred_day,
         pilevel = pilevel, nyears = nyears, nreps = nreps,
         newEvents = newEvents,
         enroll_pred_df = enroll_pred_df,
         event_pred_df = event_pred_df,
         dropout_pred_df = dropout_pred_df,
         ongoing_pred_df = ongoing_pred_df,
         event_pred_summary = s1, event_pred_plot = g1)
  }
}
