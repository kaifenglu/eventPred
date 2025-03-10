#' @title Summarize observed data
#' @description Provides an overview of the observed data, including
#' the trial start date, data cutoff date, enrollment duration,
#' number of subjects enrolled, number of events and dropouts,
#' number of subjects at risk, cumulative enrollment and event data,
#' daily enrollment rates, and Kaplan-Meier plots for time to event
#' and time to dropout.
#'
#' @param df The subject-level data, including \code{trialsdt},
#'   \code{usubjid}, \code{randdt}, and \code{cutoffdt} for enrollment
#'   prediction, as well as \code{time}, \code{event} and \code{dropout}
#'   for event prediction, and \code{treatment} coded as 1, 2,
#'   and so on, and \code{treatment_description} for prediction
#'   by treatment group.
#' @param to_predict Specifies what to predict: "enrollment only",
#'   "event only", or "enrollment and event". By default, it is set to
#'   "event only".
#' @param showplot A Boolean variable to control whether or not to
#'   show the observed data plots. By default, it is set to \code{TRUE}.
#' @param by_treatment A Boolean variable to control whether or not to
#'   summarize observed data by treatment group. By default,
#'   it is set to \code{FALSE}.
#'
#' @return A list that includes a range of summary statistics,
#' data sets, and plots depending on the value of \code{to_predict}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' observed1 <- summarizeObserved(
#'   df = interimData1,
#'   to_predict = "enrollment and event")
#'
#' observed2 <- summarizeObserved(
#'   df = interimData2,
#'   to_predict = "event only")
#'
#' @export
#'
summarizeObserved <- function(df, to_predict = "event only",
                              showplot = TRUE, by_treatment = FALSE) {

  erify::check_class(df, "data.frame")
  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))
  erify::check_bool(showplot)
  erify::check_bool(by_treatment)

  data.table::setDT(df)

  df$trialsdt <- as.Date(df$trialsdt)
  df$randdt <- as.Date(df$randdt)
  df$cutoffdt <- as.Date(df$cutoffdt)

  trialsdt = df[1, get("trialsdt")]
  cutoffdt = df[1, get("cutoffdt")]
  t0 = as.numeric(cutoffdt - trialsdt + 1)
  n0 = nrow(df)  # current number of subjects enrolled

  if (df[, any(get("randdt") < get("trialsdt"))]) {
    stop("randdt must be greater than or equal to trialsdt")
  }

  if (df[, any(get("randdt") > get("cutoffdt"))]) {
    stop("randdt must be less than or equal to cutoffdt")
  }

  if (grepl("event", to_predict, ignore.case = TRUE)) {
    d0 = df[, sum(get("event"))]  # current number of events
    c0 = df[, sum(get("dropout"))] # current number of dropouts
    r0 = df[, sum(!get("event") & !get("dropout"))] # number at risk

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

    ongoingSubjects <- df[!get("event") & !get("dropout")]

    # number of ongoing subjects with the last known date before cutoff
    rp = df[, sum(get("time") <
                    as.numeric(get("cutoffdt") - get("randdt") + 1)
                  & !get("event") & !get("dropout"))]

    # minimum calendar time for event prediction
    tp = ongoingSubjects[, min(
      as.numeric(get("randdt") - get("trialsdt")) + get("time"))]

    cutofftpdt = as.Date(tp - 1, origin = trialsdt)
  }

  if (by_treatment) {
    ngroups = df[, data.table::uniqueN(get("treatment"))]

    if (!("treatment_description" %in% names(df))) {
      df[, `:=`(treatment_description =
                  paste("Treatment", get("treatment")))]
    }

    # order treatment description based on treatment
    df[, `:=`(treatment_description = stats::reorder(as.factor(
      get("treatment_description")), get("treatment")))]

    treatment_mapping <- df[, mget(c("treatment", "treatment_description"))][
      , .SD[.N], by = "treatment"]
  } else {
    ngroups = 1
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }


  # enrollment and event data
  if (!by_treatment) {
    adsl <- df[order(get("randdt"))][
      , `:=`(n = .I, parameter = "Enrollment", date = get("randdt"))]

    # columns to keep
    cols = c("n", "parameter", "date")

    # remove duplicate
    adslu <- adsl[, .SD[.N], by = "randdt"][, mget(cols)]

    # dummy subject to initialize time axis at trial start
    adsl0 <- data.table(n = 0, parameter = "Enrollment", date = trialsdt)

    # extend enrollment information to cutoff date
    adsl1 <- adsl[.N][, `:=`(date = get("cutoffdt"))][, mget(cols)]

    if (grepl("event", to_predict, ignore.case = TRUE)) {
      # time to event data
      adtte <- data.table::copy(df)[
        , `:=`(adt = as.Date(get("time") - 1, origin = get("randdt")))][
        order(get("adt")), `:=`(n = cumsum(get("event")),
                                parameter = "Event", date = get("adt"))]

      # remove duplicate
      adtteu <- adtte[, .SD[.N], by = "adt"][, mget(cols)]

      # dummy subject to initialize time axis at trial start
      adtte0 <- data.table(n = 0, parameter = "Event", date = trialsdt)

      # combine enrollment and time to event data
      ad <- data.table::rbindlist(list(
        adsl0, adslu, adsl1, adtte0, adtteu), use.names = TRUE)
    } else {
      ad <- data.table::rbindlist(list(
        adsl0, adslu, adsl1), use.names = TRUE)
    }
  } else { # by treatment
    trtcols = c("treatment", "treatment_description")

    adsl <- df[, .SD[order(get("randdt"))], by = trtcols][
      , `:=`(n = seq_len(.N), parameter = "Enrollment", date = get("randdt")),
      by = trtcols]

    # columns to keep
    cols = c("treatment", "treatment_description", "n", "parameter", "date")

    # remove duplicate
    adslu <- adsl[, .SD[.N], by = c(
      "treatment", "treatment_description", "randdt")][, mget(cols)]

    # dummy subject to initialize time axis at trial start
    adsl0 <- data.table(
      treatment = 1:ngroups, n = 0, parameter = "Enrollment",
      date = trialsdt)[treatment_mapping, on = "treatment"]

    # extend enrollment information to cutoff date
    adsl1 <- adsl[, .SD[.N],
                  by = c("treatment", "treatment_description")][
                    , `:=`(date = get("cutoffdt"))][, mget(cols)]

    if (grepl("event", to_predict, ignore.case = TRUE)) {
      # time to event data
      trtcols = c("treatment", "treatment_description")

      adtte <- data.table::copy(df)[
        , `:=`(adt = as.Date(get("time") - 1, origin = get("randdt")))][
          , .SD[order(get("adt"))], by = trtcols][
            , `:=`(n = cumsum(get("event")),
                   parameter = "Event", date = get("adt")),
            by = trtcols]

      # remove duplicate
      adtteu <- adtte[, .SD[.N], by = c(
        "treatment", "treatment_description", "adt")][, mget(cols)]

      # dummy subject to initialize time axis at trial start
      adtte0 <- data.table(
        treatment = 1:ngroups, n = 0, parameter = "Event",
        date = trialsdt)[treatment_mapping, on = "treatment"]

      # combine enrollment and time to event data
      ad <- data.table::rbindlist(list(
        adsl0, adslu, adsl1, adtte0, adtteu), use.names = TRUE)
    } else {
      ad <- data.table::rbindlist(list(
        adsl0, adslu, adsl1), use.names = TRUE)
    }
  }


  # plot cumulative enrollment and event data
  if (!by_treatment) {
    if (ad[, data.table::uniqueN(get("parameter")) > 1]) {
      cumAccrual <- plotly::plot_ly(
        ad, x=~date, y=~n, color=~parameter, colors=c("blue", "red")) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = ""),
          yaxis = list(zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = 'h'))
    } else {
      cumAccrual <- plotly::plot_ly(ad, x=~date, y=~n) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = ""),
          yaxis = list(zeroline = FALSE),
          title = list(text = "Cumulative enrollment"))
    }

    if (showplot) print(cumAccrual)
  } else { # by treatment
    if (ad[, data.table::uniqueN(get("parameter")) > 1]) {
      cumAccrual <- plotly::plot_ly(
        ad, x=~date, y=~n, color=~parameter, colors=c("blue", "red"),
        linetype=~treatment_description) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = ""),
          yaxis = list(zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = 'h'))
    } else {
      cumAccrual <- plotly::plot_ly(
        ad, x=~date, y=~n, linetype=~treatment_description) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = ""),
          yaxis = list(zeroline = FALSE),
          legend = list(x = 0, y = 1, yanchor = "middle",
                        orientation = 'h'),
          title = list(text = "Cumulative enrollment"))
    }

    if (showplot) print(cumAccrual)
  }


  # daily enrollment plot with loess smoothing
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    t = adsl[, as.numeric(get("randdt") - get("trialsdt") + 1)]
    days = seq(1, t0)
    n = as.numeric(table(factor(t, levels = days)))

    enroll <- data.table(day = days, n = n,
                         date = as.Date(days - 1, origin = trialsdt))

    fit <- loess.smooth(enroll$date, enroll$n,
                        span = 1/3, degree = 1, family = "gaussian")

    dailyAccrual <- plotly::plot_ly(
      enroll, x=~date, y=~n, name="observed", type='scatter',
      mode='markers') %>%
      plotly::add_lines(x = fit$x, y = fit$y, name="loess") %>%
      plotly::layout(
        xaxis = list(title = ""),
        yaxis = list(zeroline = FALSE),
        title = list(text = "Daily enrollment")) %>%
      plotly::hide_legend()

    if (showplot) print(dailyAccrual)
  }


  # Kaplan-Meier plot
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!by_treatment) {
      kmfitEvent <- survival::survfit(survival::Surv(time, event) ~ 1,
                                      data = adtte)

      kmdfEvent <- data.table(time = kmfitEvent$time,
                              surv = kmfitEvent$surv)
      # add day 1
      if (kmdfEvent[, min(get("time")) > 1]) {
        kmdfEvent <- data.table::rbindlist(list(
          data.table(time = 1, surv = 1), kmdfEvent), use.names = TRUE)
      }

      kmEvent <- plotly::plot_ly(kmdfEvent, x=~time, y=~surv) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = "Days since randomization", zeroline = FALSE),
          yaxis = list(title = "Survival probability", zeroline = FALSE),
          title = list(text = "Kaplan-Meier plot for time to event"))

      if (showplot) print(kmEvent)

      # time to dropout
      kmfitDropout <- survival::survfit(survival::Surv(time, dropout) ~ 1,
                                        data = adtte)

      kmdfDropout <- data.table(time = kmfitDropout$time,
                                surv = kmfitDropout$surv)
      # add day 1
      if (kmdfDropout[, min(get("time")) > 1]) {
        kmdfDropout <- data.table::rbindlist(list(
          data.table(time = 1, surv = 1), kmdfDropout), use.names = TRUE)
      }

      kmDropout <- plotly::plot_ly(kmdfDropout, x=~time, y=~surv) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = "Days since randomization", zeroline = FALSE),
          yaxis = list(title = "Survival probability", zeroline = FALSE),
          title = list(text = "Kaplan-Meier plot for time to dropout"))

      if (showplot) print(kmDropout)
    } else { # by treatment
      kmfitEvent <- survival::survfit(survival::Surv(time, event) ~
                                        treatment, data = adtte)

      treatment <- as.numeric(substring(names(kmfitEvent$strata), 11))

      treatment_description <- treatment_mapping[
        data.table(treatment = treatment), on = "treatment",
        get("treatment_description")]

      kmdfEvent <- data.table::rbindlist(list(
        data.table(treatment = treatment,
                   treatment_description = treatment_description,
                   time = 1, surv = 1),
        data.table(treatment = rep(treatment, kmfitEvent$strata),
                   treatment_description =
                     rep(treatment_description, kmfitEvent$strata),
                   time = kmfitEvent$time,
                   surv = kmfitEvent$surv)),
        use.names = TRUE)[
          , .SD[.N], by = c("treatment", "treatment_description", "time")]

      kmEvent <- plotly::plot_ly(
        kmdfEvent, x=~time, y=~surv, linetype=~treatment_description) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = "Days since randomization", zeroline = FALSE),
          yaxis = list(title = "Survival probability", zeroline = FALSE),
          legend = list(x = 0, y = 1,  yanchor = "middle",
                        orientation = 'h'),
          title = list(text = "Kaplan-Meier plot for time to event"))

      if (showplot) print(kmEvent)

      # time to dropout
      kmfitDropout <- survival::survfit(survival::Surv(time, dropout) ~
                                          treatment, data = adtte)

      treatment <- as.numeric(substring(names(kmfitDropout$strata), 11))

      treatment_description <- treatment_mapping[
        data.table(treatment = treatment), on = "treatment",
        get("treatment_description")]

      kmdfDropout <- data.table::rbindlist(list(
        data.table(treatment = treatment,
                   treatment_description = treatment_description,
                   time = 1, surv = 1),
        data.table(treatment = rep(treatment, kmfitDropout$strata),
                   treatment_description =
                     rep(treatment_description, kmfitDropout$strata),
                   time = kmfitDropout$time,
                   surv = kmfitDropout$surv)),
        use.names = TRUE)[
          , .SD[.N], by = c("treatment", "treatment_description", "time")]

      kmDropout <- plotly::plot_ly(
        kmdfDropout, x=~time, y=~surv, linetype=~treatment_description) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = "Days since randomization", zeroline = FALSE),
          yaxis = list(title = "Survival probability", zeroline = FALSE),
          legend = list(x = 0, y = 1, yanchor = "middle",
                        orientation = 'h'),
          title = list(text = "Kaplan-Meier plot for time to dropout"))

      if (showplot) print(kmDropout)
    }
  }


  # output
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
      # enrollment and event
      list(trialsdt = trialsdt, cutoffdt = cutoffdt, t0 = t0,
           n0 = n0, d0 = d0, c0 = c0, r0 = r0, rp = rp,
           tp = tp, cutofftpdt = cutofftpdt,
           adsl = adsl, adtte = adtte,
           cum_accrual_df = ad,
           daily_accrual_df = enroll,
           event_km_df = kmdfEvent,
           dropout_km_df = kmdfDropout,
           cum_accrual_plot = cumAccrual,
           daily_accrual_plot = dailyAccrual,
           event_km_plot = kmEvent,
           dropout_km_plot = kmDropout)
    } else { # event only
      list(trialsdt = trialsdt, cutoffdt = cutoffdt, t0 = t0,
           n0 = n0, d0 = d0, c0 = c0, r0 = r0, rp = rp,
           tp = tp, cutofftpdt = cutofftpdt,
           adsl = adsl, adtte = adtte,
           cum_accrual_df = ad,
           event_km_df = kmdfEvent,
           dropout_km_df = kmdfDropout,
           cum_accrual_plot = cumAccrual,
           event_km_plot = kmEvent,
           dropout_km_plot = kmDropout)
    }
  } else { # enrollment only
    list(trialsdt = trialsdt, cutoffdt = cutoffdt, t0 = t0,
         n0 = n0, adsl = adsl,
         cum_accrual_df = ad,
         daily_accrual_df = enroll,
         cum_accrual_plot = cumAccrual,
         daily_accrual_plot = dailyAccrual)
  }
}
