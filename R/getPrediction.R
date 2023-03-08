#' @title Enrollment and event prediction
#' @description Performs enrollment and event prediction based on
#'   observed data and specified enrollment and event models.
#'
#' @param df Observed subject-level enrollment and event data.
#' @param target_n Target number of subjects to enroll.
#' @param target_d Target number of events to reach.
#' @param to_predict What to predict: enrollment only, event only,
#'   enrollment and event. Defaults to enrollment and event.
#' @param enroll_model Enrollment model: Poisson, time-decay, B-spline.
#'   Defaults to B-spine, in which case, need to specify number of knots
#'   and day lags to compute the average enrollment rate to carry forward.
#' @param nknots Number of inner knots for B-spline enrollment model.
#'   Defaults to 1.
#' @param lags Day lags for averaging enrollment rates from B-spline.
#'   Defaults to 30.
#' @param event_model Event model: exponential, Weibull, log-normal,
#'   piecewise exponential, and model averaging (of Weibull and
#'   log-normal). If piecewise exponential, need to specify
#'   number of pieces. Defaults to model averaging.
#' @param npieces Number of pieces for piecewise exponential event model.
#'   Defaults to 3.
#' @param dropout_model Dropout model: none, exponential, Weibull, log-normal.
#'   Defaults to Weibull.
#' @param fixedFollowup Whether a fixed follow-up design is used.
#'   Defaults to \code{FALSE} for variable follow-up.
#' @param followupTime Follow-up time for a fixed follow-up design.
#'   Defaults to 365 days.
#' @param pilevel Prediction interval level. Defaults to 0.90.
#' @param nreps Number of replications for simulation. Defaults to 500.
#'
#' @return A list of observed data model fits, simulated enrollment data
#'   for new subjects, and simulated event data for ongoing and new
#'   subjects.
#'
#' @examples
#'
#' ret <- getPrediction(df = observedData, target_n = 400, target_d = 120,
#'                      to_predict = "enrollment and event",
#'                      enroll_model = "b-spline", nknots = 1, lags = 30,
#'                      event_model = "piecewise exponential", npieces = 3,
#'                      dropout_model = "exponential",
#'                      fixedFollowup = FALSE, followupTime = 365,
#'                      pilevel = 0.90, nreps = 500)
#'
#' @export
#'
getPrediction <- function(
    df, target_n, target_d,
    to_predict = "enrollment and event",
    enroll_model = "B-spline", nknots = 1, lags = 30,
    event_model = "model averaging", npieces = 3,
    dropout_model = "weibull",
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nreps = 500) {

  erify::check_class(df, "data.frame")
  erify::check_n(target_n)
  erify::check_n(target_d)
  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))
  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay", "b-spline"))
  erify::check_n(nknots)
  erify::check_n(lags, zero=TRUE)
  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-normal",
                         "piecewise exponential", "model averaging"))
  erify::check_n(npieces)
  erify::check_content(tolower(dropout_model),
                       c("exponential", "weibull", "log-normal"))
  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(nreps)

  # make a copy of the input data set
  dfin <- df

  # trial start date and cutoff date for the data set
  trialsdt = min(dfin$RANDDT)
  cutoffdt = dfin$CUTOFFDT[1]

  # summarize observed data
  observed <- summarizeObserved(dfin, to_predict, dropout_model)

  # fit and predict enrollment
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    erify::check_n(target_n - observed$n0,
                   supplement = "Enrollment target reached.")

    fitEnr <- fitEnrollment(observed$adsl, enroll_model, nknots)
    predEnr <- predictEnrollment(target_n, observed$adsl, fitEnr,
                                 lags, pilevel, nreps)
  }

  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    erify::check_n(target_d - observed$d0,
                   supplement = "Event target reached.")

    fitEvt <- fitEvent(observed$adtte, observed$kmdfEvt, event_model, npieces)

    if (tolower(dropout_model) != "none") {
      fitDrp <- fitDropout(observed$adtte, observed$kmdfDrp, dropout_model)
      if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
        predEvt <- predictEvent(target_d, observed$adtte, predEnr$newSubjects,
                                fitEvt, fitDrp, fixedFollowup, followupTime,
                                pilevel, lags, nreps)
      } else {
        predEvt <- predictEvent(target_d, observed$adtte, NULL,
                                fitEvt, fitDrp, fixedFollowup, followupTime,
                                pilevel, lags, nreps)
      }
    } else {
      if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
        predEvt <- predictEvent(target_d, observed$adtte, predEnr$newSubjects,
                                fitEvt, NULL, fixedFollowup, followupTime,
                                pilevel, lags, nreps)
      } else {
        predEvt <- predictEvent(target_d, observed$adtte, NULL,
                                fitEvt, NULL, fixedFollowup, followupTime,
                                pilevel, lags, nreps)
      }
    }
  }

  if (tolower(to_predict) == "enrollment only") {
    list(observed = observed, fitEnr = fitEnr, predEnr = predEnr)
  } else if (tolower(to_predict) == "event only") {
    if (tolower(dropout_model) != "none") {
      list(observed = observed, fitEvt = fitEvt, fitDrp = fitDrp,
           predEvt = predEvt)
    } else {
      list(observed = observed, fitEvt = fitEvt, predEvt = predEvt)
    }
  } else if (tolower(to_predict) == "enrollment and event") {
    if (tolower(dropout_model) != "none") {
      list(observed = observed, fitEnr = fitEnr, predEnr = predEnr,
           fitEvt = fitEvt, fitDrp = fitDrp, predEvt = predEvt)
    } else {
      list(observed = observed, fitEnr = fitEnr, predEnr = predEnr,
           fitEvt = fitEvt, predEvt = predEvt)
    }
  }

}
