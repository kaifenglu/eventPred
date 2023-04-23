#' @title Summarize observed data
#' @description Provides an overview of the observed data, including
#' the trial start date, data cutoff date, number of subjects
#' enrolled, enrollment duration, number of events and dropouts,
#' number of subjects at risk, cumulative enrollment and event data,
#' daily enrollment rates, and Kaplan-Meier plots for time to event
#' and time to dropout.
#'
#' @param df The subject-level data, including \code{trialsdt},
#'   \code{randdt}, \code{cutoffdt} for enrollment prediction,
#'   as well as \code{time}, \code{event} and \code{dropout}
#'   for event prediction.
#' @param to_predict Specifies what to predict: "enrollment only",
#'   "event only", or "enrollment and event". By default, it is set to
#'   "enrollment and event".
#' @param showplot A Boolean variable to control whether or not to
#'   show the observed data plots. By default, it is set to \code{TRUE}.
#'
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
summarizeObserved <- function(df, to_predict = "enrollment and event",
                              showplot = TRUE) {
  erify::check_class(df, "data.frame")
  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))

  df <- dplyr::as_tibble(df)
  names(df) <- tolower(names(df))
  df$trialsdt <- as.Date(df$trialsdt)
  df$randdt <- as.Date(df$randdt)
  df$cutoffdt <- as.Date(df$cutoffdt)

  trialsdt = df$trialsdt[1]
  cutoffdt = df$cutoffdt[1]
  t0 = as.numeric(cutoffdt - trialsdt + 1)
  n0 = nrow(df)  # current number of subjects enrolled

  if (grepl("event", to_predict, ignore.case = TRUE)) {
    d0 = sum(df$event)  # current number of events
    c0 = sum(df$dropout) # current number of dropouts
    r0 = sum(!(df$event | df$dropout)) # number of subjects at risk
  }

  # enrollment data
  adsl <- df %>%
    dplyr::arrange(.data$randdt) %>%
    dplyr::mutate(adt = as.Date(.data$time - 1, origin = .data$randdt),
                  n = dplyr::row_number(),
                  parameter = "Enrollment",
                  date = .data$randdt)

  adslu <- adsl %>%
    dplyr::group_by(.data$randdt) %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::ungroup()

  # extend enrollment information to cutoff date
  adsl1 <- adsl %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::mutate(date = cutoffdt)



  if (grepl("event", to_predict, ignore.case = TRUE)) {
    # time to event data
    adtte <- df %>%
      dplyr::mutate(adt = as.Date(.data$time - 1, origin = .data$randdt)) %>%
      dplyr::arrange(.data$adt) %>%
      dplyr::mutate(n = cumsum(.data$event),
                    parameter = "Event",
                    date = .data$adt)

    adtteu <- adtte %>%
      dplyr::group_by(.data$adt) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::ungroup()

    # dummy subject to initialize time to event axis at trial start
    adtte0 <- adtte %>%
      dplyr::slice(1) %>%
      dplyr::mutate(randdt = trialsdt, adt = trialsdt, time = 1,
                    event = 0, dropout = 0,
                    n = 0, parameter = "Event", date = trialsdt)

    # combine enrollment and time to event data
    ad <- adslu %>%
      dplyr::bind_rows(adsl1) %>%
      dplyr::bind_rows(adtte0) %>%
      dplyr::bind_rows(adtteu)
  } else {
    ad <- adslu %>%
      dplyr::bind_rows(adsl1)
  }


  # plot cumulative enrollment and event data
  if (length(unique(ad$parameter)) > 1) {
    cumAccrual <- plotly::plot_ly(
      ad, x=~date, y=~n, color=~parameter, colors=c("blue", "red")) %>%
      plotly::add_lines(line = list(shape = "hv")) %>%
      plotly::layout(
        xaxis = list(title = ""),
        yaxis = list(zeroline = FALSE),
        legend = list(x = 0, y = 1.05, yanchor = "bottom", orientation = 'h'))
  } else {
    cumAccrual <- plotly::plot_ly(ad, x=~date, y=~n) %>%
      plotly::add_lines(line = list(shape = "hv")) %>%
      plotly::layout(
        xaxis = list(title = ""),
        yaxis = list(zeroline = FALSE),
        title = list(text = "Cumulative enrollment"))
  }

  if (showplot) print(cumAccrual)

  # daily enrollment plot with loess smoothing
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    t = as.numeric(adsl$randdt - trialsdt + 1)
    days = seq(1, t0)
    n = as.numeric(table(factor(t, levels = days)))

    enroll <- dplyr::tibble(day = days, n = n) %>%
      dplyr::mutate(date = as.Date(.data$day - 1, origin = trialsdt))

    fit <- loess.smooth(enroll$date, enroll$n,
                        span = 1/3, degree = 1,
                        family = "gaussian")

    dailyAccrual <- plotly::plot_ly(enroll, x=~date, y=~n, name="observed",
                                    type='scatter', mode='markers') %>%
      plotly::add_lines(x = fit$x, y = fit$y, name="loess") %>%
      plotly::layout(xaxis = list(title = ""),
                     yaxis = list(zeroline = FALSE),
                     title = list(text = "Daily enrollment")) %>%
      plotly::hide_legend()
    if (showplot) print(dailyAccrual)
  }


  # Kaplan-Meier plot
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    kmfitEvent <- survival::survfit(survival::Surv(time, event) ~ 1,
                                    data = adtte)

    kmdfEvent <- dplyr:: tibble(time = 0, surv = 1) %>%
      dplyr::bind_rows(dplyr::tibble(time = kmfitEvent$time,
                                     surv = kmfitEvent$surv))

    kmEvent <- plotly::plot_ly(kmdfEvent, x=~time, y=~surv) %>%
      plotly::add_lines(line = list(shape = "hv")) %>%
      plotly::layout(xaxis = list(title = "Days since randomization",
                                  zeroline = FALSE),
                     yaxis = list(title = "Survival probability",
                                  zeroline = FALSE),
                     title = list(
                       text = "Kaplan-Meier plot for time to event"))
    if (showplot) print(kmEvent)


    # time to dropout
    kmfitDropout <- survival::survfit(survival::Surv(time, dropout) ~ 1,
                                      data = adtte)
    kmdfDropout <- dplyr::tibble(time = 0, surv = 1) %>%
      dplyr::bind_rows(dplyr::tibble(time = kmfitDropout$time,
                                     surv = kmfitDropout$surv))

    kmDropout <- plotly::plot_ly(kmdfDropout, x=~time, y=~surv) %>%
      plotly::add_lines(line = list(shape = "hv")) %>%
      plotly::layout(xaxis = list(title = "Days since randomization",
                                  zeroline = FALSE),
                     yaxis = list(title = "Survival probability",
                                  zeroline = FALSE),
                     title = list(
                       text = "Kaplan-Meier plot for time to dropout"))
    if (showplot) print(kmDropout)

  }


  # output
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
      # enrollment and event
      list(trialsdt = trialsdt, cutoffdt = cutoffdt, t0 = t0,
           n0 = n0, d0 = d0, c0 = c0, r0 = r0, adsl = adsl,
           adtte = adtte, event_km_df = kmdfEvent,
           dropout_km_df = kmdfDropout,
           cum_accrual_plot = cumAccrual,
           daily_accrual_plot = dailyAccrual,
           event_km_plot = kmEvent,
           dropout_km_plot = kmDropout)
    } else { # event only
      list(trialsdt = trialsdt, cutoffdt = cutoffdt, t0 = t0,
           n0 = n0, d0 = d0, c0 = c0, r0 = r0, adsl = adsl,
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
