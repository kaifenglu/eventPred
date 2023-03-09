#' @title Summarize observed data
#' @description Provides an overview of the observed data, including
#' the trial start date, data cutoff date, number of subjects
#' enrolled, enrollment duration, number of events and dropouts,
#' number of subjects at risk, cumulative enrollment and event data,
#' daily enrollment rate, and Kaplan-Meier plots for time to event
#' and time to dropout.
#'
#' @param df The subject-level data, including \code{randdt} and
#'   \code{cutoffdt} for enrollment prediction,
#'   as well as \code{time}, \code{event} and \code{dropout}
#'   for event prediction.
#' @param to_predict Specifies what to predict: "enrollment only",
#'   "event only", or "enrollment and event". By default, it is set to
#'   "enrollment and event".
#' @param dropout_model The dropout model with options including "none",
#'   "exponential", "Weibull", and "log-normal". By default, it is
#'   set to "Weibull".
#' @return A list that includes a range of summary statistics and
#' data sets depending on the value of \code{to_predict}.
#'
#' @examples
#'
#' observed <- summarizeObserved(df = observedData,
#'                               to_predict = "enrollment and event")
#'
#' @export
#'
summarizeObserved <- function(df, to_predict = "enrollment and event",
                              dropout_model = "weibull") {
  erify::check_class(df, "data.frame")
  erify::check_content(tolower(to_predict),
                       c("enrollment only", "event only",
                         "enrollment and event"))
  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-normal"))

  names(df) <- tolower(names(df))
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
    arrange(.data$randdt) %>%
    mutate(n = row_number(),
           parameter = "subjects",
           date = .data$randdt) %>%
    mutate(year = format(.data$date, format = "%Y"))


  if (grepl("event", to_predict, ignore.case = TRUE)) {
    # time to event data
    adtte <- df %>%
      mutate(adt = as.Date(.data$time - 1, origin = .data$randdt)) %>%
      arrange(.data$adt) %>%
      mutate(n = cumsum(.data$event),
             parameter = "events",
             date = .data$adt) %>%
      mutate(year = format(.data$date, format = "%Y"))

    # dummy subject to initialize time to event axis at trial start
    adtte0 <- first(df) %>%
      mutate(randdt = trialsdt, adt = trialsdt,
             event = 0, dropout = 0,
             n = 0, parameter = "events", date = trialsdt) %>%
      mutate(year = format(.data$date, format = "%Y"))

    # combine enrollment and time to event data
    ad <- adsl %>%
      bind_rows(adtte0) %>%
      bind_rows(adtte)

    ylab = "Subjects / Events"
    title = "Observed cumulative subjects and events over time"
  } else {
    ad <- adsl
    ylab = "Subjects"
    title = "Observed cumulative subjects over time"
  }


  # use number of months between first and last dates to determine ticks
  n_months = lubridate::interval(min(ad$date), max(ad$date)) %/% months(1)
  bw = fbw(n_months)

  # plot cumulative enrollment and event data
  g1 <- ggplot() +
    geom_step(data = ad, aes(x = .data$date, y = .data$n,
                             group = .data$parameter)) +
    scale_x_date(name = NULL,
                 labels = scales::date_format("%b"),
                 breaks = scales::breaks_width(bw),
                 minor_breaks = NULL,
                 expand = c(0.01, 0.01)) +
    labs(y = ylab, title = title) +
    theme_bw()

  # generate the year labels
  g2 <- flabel(ad, trialsdt)

  # stack them together
  cumAccrual <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
  print(cumAccrual)



  # daily enrollment plot with loess smoothing
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    adsl <- adsl %>%
      mutate(time = as.numeric(.data$randdt - trialsdt + 1))

    days = seq(1, t0)
    n = sapply(days, function(i) sum(adsl$time == i))
    enroll <- tibble(day = days, n = n)

    dailyAccrual <- ggplot(data = enroll, aes(x = .data$day, y = .data$n)) +
      geom_point() +
      geom_smooth(formula = y ~ x, method = loess, se = FALSE) +
      labs(x = "Days since trial start",
           y = "Subjects enrolled daily",
           title = "Daily subject enrollment") +
      theme_bw()
    print(dailyAccrual)
  }


  # Kaplan-Meier plot
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    adtte <- adtte %>%
      mutate(time = as.numeric(.data$adt - .data$randdt + 1))

    kmfitEvt <- survival::survfit(survival::Surv(time, event) ~ 1,
                                  data = adtte)
    kmdfEvt <- tibble(time = kmfitEvt$time, surv = kmfitEvt$surv)
    kmdfEvt <- tibble(time = 0, surv = 1) %>%
      bind_rows(kmdfEvt)

    kmEvent <- ggplot() +
      geom_step(data = kmdfEvt, aes(x = .data$time, y = .data$surv)) +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Kaplan-Meier plot for time to event") +
      theme_bw()
    print(kmEvent)

    # time to dropout
    if (tolower(dropout_model) != "none") {
      kmfitDrp <- survival::survfit(survival::Surv(time, dropout) ~ 1,
                                    data = adtte)
      kmdfDrp <- tibble(time = kmfitDrp$time, surv = kmfitDrp$surv)
      kmdfDrp <- tibble(time = 0, surv = 1) %>%
        bind_rows(kmdfDrp)

      kmDropout <- ggplot() +
        geom_step(data = kmdfDrp, aes(x = .data$time, y = .data$surv)) +
        labs(x = "Days since randomization",
             y = "Survival probability",
             title = "Kaplan-Meier plot for time to dropout") +
        theme_bw()
      print(kmDropout)
    }
  }


  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (tolower(dropout_model) != "none") {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           n0 = n0, t0 = t0, d0 = d0, c0 = c0, r0 = r0, adsl = adsl,
           adtte = adtte, kmdfEvt = kmdfEvt, kmdfDrp = kmdfDrp)
    } else {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           n0 = n0, t0 = t0, d0 = d0, r0 = r0, adsl = adsl,
           adtte = adtte, kmdfEvt = kmdfEvt)
    }
  } else {
    list(trialsdt = trialsdt, cutoffdt = cutoffdt,
         n0 = n0, t0 = t0, adsl = adsl)
  }
}
