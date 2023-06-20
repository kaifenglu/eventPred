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
#' @examples
#'
#' observed1 <- summarizeObserved(df = interimData1,
#'                                to_predict = "enrollment and event")
#'
#' observed2 <- summarizeObserved(df = interimData2,
#'                                to_predict = "event only")
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

  df <- dplyr::as_tibble(df)
  names(df) <- tolower(names(df))
  df$trialsdt <- as.Date(df$trialsdt)
  df$randdt <- as.Date(df$randdt)
  df$cutoffdt <- as.Date(df$cutoffdt)

  trialsdt = df$trialsdt[1]
  cutoffdt = df$cutoffdt[1]
  t0 = as.numeric(cutoffdt - trialsdt + 1)
  n0 = nrow(df)  # current number of subjects enrolled

  if (any(df$randdt < trialsdt)) {
    stop("randdt must be greater than or equal to trialsdt.")
  }

  if (any(df$randdt > cutoffdt)) {
    stop("randdt must be less than or equal to cutoffdt.")
  }

  if (grepl("event", to_predict, ignore.case = TRUE)) {
    d0 = sum(df$event)  # current number of events
    c0 = sum(df$dropout) # current number of dropouts
    r0 = sum(!(df$event | df$dropout)) # number of subjects at risk

    # number of ongoing subjects with the last known date before cutoff
    rp = sum((df$time < as.numeric(cutoffdt - df$randdt + 1)) &
               !(df$event | df$dropout))

    if (any(df$time < 1)) {
      stop("time must be greater than or equal to 1.")
    }

    if (any(df$event == 1 & df$dropout == 1)) {
      stop("event and dropout cannot both be equal to 1 simultaneously.")
    }

    if (any(df$time > as.numeric(cutoffdt - df$randdt + 1))) {
      stop("time must be less than or equal to cutoffdt - randdt + 1.")
    }

    ongoingSubjects <- df %>%
      dplyr::filter(.data$event == 0 & .data$dropout == 0)

    # minimum calendar time for event prediction
    tp = min(as.numeric(ongoingSubjects$randdt - trialsdt + 1) +
               ongoingSubjects$time - 1)
    cutofftpdt = as.Date(tp - 1, origin = trialsdt)
  }

  if (by_treatment) {
    ngroups = length(table(df$treatment))
    if (!("treatment_description" %in% names(df))) {
      df <- df %>% dplyr::mutate(
        treatment_description = paste0("Treatment ", .data$treatment))
    }

    # order treatment description based on treatment
    df$treatment_description = stats::reorder(
      as.factor(df$treatment_description), df$treatment)

    treatment_mapping <- df %>%
      dplyr::select(.data$treatment, .data$treatment_description) %>%
      dplyr::arrange(.data$treatment) %>%
      dplyr::group_by(.data$treatment) %>%
      dplyr::slice(dplyr::n())
  } else {
    ngroups = 1
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }

  # enrollment and event data
  if (!by_treatment) {
    adsl <- df %>%
      dplyr::arrange(.data$randdt) %>%
      dplyr::mutate(n = dplyr::row_number(),
                    parameter = "Enrollment",
                    date = .data$randdt)

    # remove duplicate
    adslu <- adsl %>%
      dplyr::group_by(.data$randdt) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$n, .data$parameter, .data$date)

    # dummy subject to initialize time axis at trial start
    adsl0 <- dplyr::tibble(n = 0, parameter = "Enrollment", date = trialsdt)

    # extend enrollment information to cutoff date
    adsl1 <- adsl %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::mutate(date = cutoffdt) %>%
      dplyr::select(.data$n, .data$parameter, .data$date)


    if (grepl("event", to_predict, ignore.case = TRUE)) {
      # time to event data
      adtte <- df %>%
        dplyr::mutate(adt = as.Date(.data$time - 1, origin = .data$randdt)) %>%
        dplyr::arrange(.data$adt) %>%
        dplyr::mutate(n = cumsum(.data$event),
                      parameter = "Event",
                      date = .data$adt)

      # remove duplicate
      adtteu <- adtte %>%
        dplyr::group_by(.data$adt) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::select(.data$n, .data$parameter, .data$date)

      # dummy subject to initialize time axis at trial start
      adtte0 <- dplyr::tibble(n = 0, parameter = "Event", date = trialsdt)

      # combine enrollment and time to event data
      ad <- adsl0 %>%
        dplyr::bind_rows(adslu) %>%
        dplyr::bind_rows(adsl1) %>%
        dplyr::bind_rows(adtte0) %>%
        dplyr::bind_rows(adtteu)
    } else {
      ad <- adsl0 %>%
        dplyr::bind_rows(adslu) %>%
        dplyr::bind_rows(adsl1)
    }
  } else { # by treatment
    adsl <- df %>%
      dplyr::group_by(.data$treatment, .data$treatment_description) %>%
      dplyr::arrange(.data$randdt) %>%
      dplyr::mutate(n = dplyr::row_number(),
                    parameter = "Enrollment",
                    date = .data$randdt) %>%
      dplyr::ungroup()

    # remove duplicate
    adslu <- adsl %>%
      dplyr::group_by(.data$treatment, .data$treatment_description,
                      .data$randdt) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$treatment, .data$treatment_description,
                    .data$n, .data$parameter, .data$date)

    # dummy subject to initialize time axis at trial start
    adsl0 <- dplyr::tibble(treatment = 1:ngroups,
                           n = 0,
                           parameter = "Enrollment",
                           date = trialsdt) %>%
      dplyr::left_join(treatment_mapping, by = "treatment")

    # extend enrollment information to cutoff date
    adsl1 <- adsl %>%
      dplyr::group_by(.data$treatment, .data$treatment_description) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::mutate(date = cutoffdt) %>%
      dplyr::select(.data$treatment, .data$treatment_description,
                    .data$n, .data$parameter, .data$date) %>%
      dplyr::ungroup()


    if (grepl("event", to_predict, ignore.case = TRUE)) {
      # time to event data
      adtte <- df %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::mutate(adt = as.Date(.data$time - 1, origin = .data$randdt)) %>%
        dplyr::arrange(.data$adt) %>%
        dplyr::mutate(n = cumsum(.data$event),
                      parameter = "Event",
                      date = .data$adt) %>%
        dplyr::ungroup()

      # remove duplicate
      adtteu <- adtte %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$adt) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::select(.data$treatment, .data$treatment_description,
                      .data$n, .data$parameter, .data$date)

      # dummy subject to initialize time axis at trial start
      adtte0 <- dplyr::tibble(treatment = 1:ngroups,
                              n = 0,
                              parameter = "Event",
                              date = trialsdt) %>%
        dplyr::left_join(treatment_mapping, by = "treatment")

      # combine enrollment and time to event data
      ad <- adsl0 %>%
        dplyr::bind_rows(adslu) %>%
        dplyr::bind_rows(adsl1) %>%
        dplyr::bind_rows(adtte0) %>%
        dplyr::bind_rows(adtteu)
    } else {
      ad <- adsl0 %>%
        dplyr::bind_rows(adslu) %>%
        dplyr::bind_rows(adsl1)
    }
  }

  # plot cumulative enrollment and event data
  if (!by_treatment) {
    if (length(unique(ad$parameter)) > 1) {
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
    if (length(unique(ad$parameter)) > 1) {
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
    t = as.numeric(adsl$randdt - trialsdt + 1)
    days = seq(1, t0)
    n = as.numeric(table(factor(t, levels = days)))

    enroll <- dplyr::tibble(day = days, n = n) %>%
      dplyr::mutate(date = as.Date(.data$day - 1, origin = trialsdt))

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

      kmdfEvent <- dplyr::tibble(time = kmfitEvent$time,
                                 surv = kmfitEvent$surv)
      # add day 1
      if (min(kmdfEvent$time) > 1) {
        kmdfEvent <- dplyr::tibble(time = 1, surv = 1) %>%
          dplyr::bind_rows(kmdfEvent)
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

      kmdfDropout <- dplyr::tibble(time = kmfitDropout$time,
                                   surv = kmfitDropout$surv)
      if (min(kmdfDropout$time) > 1) {
        kmdfDropout <- dplyr::tibble(time = 1, surv = 1) %>%
          dplyr::bind_rows(kmdfDropout)
      }

      kmDropout <- plotly::plot_ly(kmdfDropout, x=~time, y=~surv) %>%
        plotly::add_lines(line = list(shape = "hv")) %>%
        plotly::layout(
          xaxis = list(title = "Days since randomization", zeroline = FALSE),
          yaxis = list(title = "Survival probability", zeroline = FALSE),
          title = list(text = "Kaplan-Meier plot for time to dropout"))

      if (showplot) print(kmDropout)
    } else { # by treatment
      kmfitEvent <- survival::survfit(survival::Surv(time, event) ~ treatment,
                                      data = adtte)
      treatment <- as.numeric(substring(attr(kmfitEvent$strata, "names"), 11))
      treatment_description <-
        (treatment_mapping %>% dplyr::right_join(dplyr::tibble(
          treatment = treatment), by = "treatment"))$treatment_description

      kmdfEvent <- dplyr::tibble(
        treatment = treatment, treatment_description = treatment_description,
        time = 1, surv = 1) %>%
        dplyr::bind_rows(dplyr::tibble(
          treatment = rep(treatment, kmfitEvent$strata),
          treatment_description = rep(treatment_description,
                                      kmfitEvent$strata),
          time = kmfitEvent$time,
          surv = kmfitEvent$surv)) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$time) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()

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

      treatment <- as.numeric(substring(attr(kmfitDropout$strata, "names"),
                                        11))
      treatment_description <-
        (treatment_mapping %>% dplyr::right_join(dplyr::tibble(
          treatment = treatment), by = "treatment"))$treatment_description

      kmdfDropout <- dplyr::tibble(
        treatment = treatment, treatment_description = treatment_description,
        time = 1, surv = 1) %>%
        dplyr::bind_rows(dplyr::tibble(
          treatment = rep(treatment, kmfitDropout$strata),
          treatment_description = rep(treatment_description,
                                      kmfitDropout$strata),
          time = kmfitDropout$time,
          surv = kmfitDropout$surv)) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$time) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()

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
           adsl = adsl,
           adtte = adtte, event_km_df = kmdfEvent,
           dropout_km_df = kmdfDropout,
           cum_accrual_plot = cumAccrual,
           event_km_plot = kmEvent,
           dropout_km_plot = kmDropout)
    }
  } else { # enrollment only
    list(trialsdt = trialsdt, cutoffdt = cutoffdt, t0 = t0,
         n0 = n0, adsl = adsl,
         cum_accrual_plot = cumAccrual,
         daily_accrual_plot = dailyAccrual)
  }
}
