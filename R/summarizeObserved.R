#' @title Summarize observed data
#' @description Provides an overview of the observed data, including
#' the trial start date, data cutoff date, number of subjects
#' enrolled, enrollment duration, number of events and dropouts,
#' number of subjects at risk, cumulative enrollment and event data,
#' daily enrollment rates, and Kaplan-Meier plots for time to event
#' and time to dropout.
#'
#' @param df The subject-level data, including \code{randdt} and
#'   \code{cutoffdt} for enrollment prediction,
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
  df$randdt <- as.Date(df$randdt)
  df$cutoffdt <- as.Date(df$cutoffdt)

  trialsdt = min(df$randdt)
  cutoffdt = df$cutoffdt[1]
  n0 = nrow(df)  # current number of subjects enrolled
  t0 = as.numeric(cutoffdt - trialsdt + 1)

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
                  date = .data$randdt) %>%
    dplyr::mutate(year = format(.data$date, format = "%Y"))

  adsl1 <- adsl %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::mutate(date = cutoffdt,
                  year = format(.data$date, format = "%Y"))

  if (grepl("event", to_predict, ignore.case = TRUE)) {
    # time to event data
    adtte <- df %>%
      dplyr::mutate(adt = as.Date(.data$time - 1, origin = .data$randdt)) %>%
      dplyr::arrange(.data$adt) %>%
      dplyr::mutate(n = cumsum(.data$event),
                    parameter = "Event",
                    date = .data$adt) %>%
      dplyr::mutate(year = format(.data$date, format = "%Y"))

    # dummy subject to initialize time to event axis at trial start
    adtte0 <- df %>% dplyr::slice(1) %>%
      dplyr::mutate(randdt = trialsdt, adt = trialsdt, time = 1,
                    event = 0, dropout = 0,
                    n = 0, parameter = "Event", date = trialsdt) %>%
      dplyr::mutate(year = format(.data$date, format = "%Y"))

    # combine enrollment and time to event data
    ad <- adsl %>%
      dplyr::bind_rows(adsl1) %>%
      dplyr::bind_rows(adtte0) %>%
      dplyr::bind_rows(adtte)
  } else {
    ad <- adsl %>%
      dplyr::bind_rows(adsl1)
  }

  df2 <- dplyr::tibble(
    parameter = c("Enrollment", "Event", "Dropout", "Ongoing"),
    colorvalues = c("blue", "red", "#D55E00", "#009E73"))

  # only show legends for parameters appearing in data set
  df3 <- df2 %>%
    dplyr::filter(.data$parameter %in% unique(ad$parameter))

  # use number of months between first and last dates to determine ticks
  n_months = lubridate::interval(min(ad$date), max(ad$date)) %/% months(1)
  bw = fbw(n_months)

  # plot cumulative enrollment and event data
  if (length(unique(ad$parameter)) > 1) {
    g1 <- ggplot2::ggplot() +
      ggplot2::geom_step(data = ad, ggplot2::aes(
        x = .data$date, y = .data$n,
        group = .data$parameter, color = .data$parameter)) +
      ggplot2::scale_x_date(name = NULL,
                            labels = scales::date_format("%b"),
                            breaks = scales::breaks_width(bw),
                            minor_breaks = NULL,
                            expand = c(0.01, 0.01)) +
      ggplot2::scale_color_manual(name = NULL,
                                  values = df3$colorvalues) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  } else {
    g1 <- ggplot2::ggplot() +
      ggplot2::geom_step(data = ad, ggplot2::aes(
        x = .data$date, y = .data$n,
        group = .data$parameter)) +
      ggplot2::scale_x_date(name = NULL,
                            labels = scales::date_format("%b"),
                            breaks = scales::breaks_width(bw),
                            minor_breaks = NULL,
                            expand = c(0.01, 0.01)) +
      ggplot2::labs(y = "Subjects", title = "Cumulative enrollment") +
      ggplot2::theme_bw()
  }

  # generate the year labels
  g2 <- flabel(ad, trialsdt)

  # stack them together
  cumAccrual <- g1 + g2 + patchwork::plot_layout(nrow = 2,
                                                 heights = c(15, 1))
  if (showplot) print(cumAccrual)

  # daily enrollment plot with loess smoothing
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    t = as.numeric(adsl$randdt - trialsdt + 1)
    days = seq(1, t0)
    n = as.numeric(table(factor(t, levels = days)))

    enroll <- dplyr::tibble(day = days, n = n) %>%
      dplyr::mutate(date = as.Date(.data$day - 1, origin = trialsdt),
                    year = format(.data$date, format = "%Y"))

    n_months = lubridate::interval(min(enroll$date),
                                   max(enroll$date)) %/% months(1)
    bw = fbw(n_months)

    g1 <- ggplot2::ggplot(data = enroll,
                          ggplot2::aes(x = .data$date, y = .data$n)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = FALSE) +
      ggplot2::scale_x_date(name = NULL,
                            labels = scales::date_format("%b"),
                            breaks = scales::breaks_width(bw),
                            minor_breaks = NULL,
                            expand = c(0.01, 0.01)) +
      ggplot2::labs(y = "Subjects", title = "Daily enrollment") +
      ggplot2::theme_bw()

    g2 <- flabel(enroll, trialsdt)

    dailyAccrual <- g1 + g2 +
      patchwork::plot_layout(nrow = 2, heights = c(15, 1))
    if (showplot) print(dailyAccrual)
  }


  # Kaplan-Meier plot
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    kmfitEvent <- survival::survfit(survival::Surv(time, event) ~ 1,
                                    data = adtte)

    kmdfEvent <- dplyr:: tibble(time = 0, surv = 1) %>%
      dplyr::bind_rows(dplyr::tibble(time = kmfitEvent$time,
                                     surv = kmfitEvent$surv))

    kmEvent <- ggplot2::ggplot() +
      ggplot2::geom_step(data = kmdfEvent, ggplot2::aes(x = .data$time,
                                                        y = .data$surv)) +
      ggplot2::labs(x = "Days since randomization",
                    y = "Survival probability",
                    title = "Kaplan-Meier plot for time to event") +
      ggplot2::theme_bw()
    if (showplot) print(kmEvent)

    # time to dropout
    kmfitDropout <- survival::survfit(survival::Surv(time, dropout) ~ 1,
                                      data = adtte)
    kmdfDropout <- dplyr::tibble(time = 0, surv = 1) %>%
      dplyr::bind_rows(dplyr::tibble(time = kmfitDropout$time,
                                     surv = kmfitDropout$surv))

    kmDropout <- ggplot2::ggplot() +
      ggplot2::geom_step(data = kmdfDropout, ggplot2::aes(x = .data$time,
                                                          y = .data$surv)) +
      ggplot2::labs(x = "Days since randomization",
                    y = "Survival probability",
                    title = "Kaplan-Meier plot for time to dropout") +
      ggplot2::theme_bw()
    if (showplot) print(kmDropout)
  }


  # output
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           n0 = n0, t0 = t0, d0 = d0, c0 = c0, r0 = r0, adsl = adsl,
           adtte = adtte, event_km_df = kmdfEvent,
           dropout_km_df = kmdfDropout,
           cum_accrual_plot = cumAccrual,
           daily_accrual_plot = dailyAccrual,
           event_km_plot = kmEvent,
           dropout_km_plot = kmDropout)
    } else {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           n0 = n0, t0 = t0, d0 = d0, c0 = c0, r0 = r0, adsl = adsl,
           adtte = adtte, event_km_df = kmdfEvent,
           dropout_km_df = kmdfDropout,
           cum_accrual_plot = cumAccrual,
           event_km_plot = kmEvent,
           dropout_km_plot = kmDropout)
    }
  } else {
    if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           n0 = n0, t0 = t0, adsl = adsl,
           cum_accrual_plot = cumAccrual,
           daily_accrual_plot = dailyAccrual)
    } else {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           n0 = n0, t0 = t0, adsl = adsl,
           cum_accrual_plot = cumAccrual)
    }
  }
}
