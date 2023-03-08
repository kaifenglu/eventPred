#' @title Predict event
#' @description Use the fitted time-to-event and time-to-dropout model to
#'   generate event and dropout times for ongoing subjects (conditional on
#'   their current study exposure) and new subjects, and provide the
#'   prediction interval for time to reach the target number of events.
#'
#' @param target_d Target number of events to reach.
#' @param df Observed subject-level enrollment and event data.
#' @param newSubjects Enrollment data for new subjects. Defaults to
#'   \code{NULL} implying the completion of subject enrollment.
#' @param fitEvt Event model fit to observed data.
#' @param fitDrp Dropout model fit to observed data. Defaults to
#'   \code{NULL} implying the absence of dropout.
#' @param fixedFollowup Whether a fixed follow-up design is used.
#'   Defaults to \code{FALSE} for variable follow-up.
#' @param followupTime Follow-up time for fixed follow-up design. Defaults
#'   to 365 days.
#' @param pilevel Prediction interval level. Defaults to 0.90.
#' @param lags Day lags for averaging enrollment rates from B-spline.
#'   Defaults to 30.
#' @param nreps Number of replications for simulation. Defaults to 500.
#'
#' @return A list of prediction results consisting of the median time to
#'   reach event target, date of median, lower and upper percentiles
#'   to reach event target, and simulated event data for ongoing and new
#'   subjects.
#'
#' @examples
#'
#' observed <- summarizeObserved(df = observedData,
#'                               to_predict = "enrollment and event")
#'
#' fitEnr <- fitEnrollment(df = observed$adsl, enroll_model = "b-spline",
#'                         nknots = 1)
#'
#' predEnr <- predictEnrollment(target_n = 400, df = observed$adsl,
#'                              fit = fitEnr, lags = 30, pilevel = 0.90,
#'                              nreps = 500)
#'
#' fitEvt <- fitEvent(df = observed$adtte, kmdf = observed$kmdfEvt,
#'                    event_model = "piecewise exponential", npieces = 3)
#'
#' fitDrp <- fitDropout(df = observed$adtte, kmdf = observed$kmdfDrp,
#'                      dropout_model = "exponential")
#'
#' predEvt <- predictEvent(target_d = 120, df = observed$adtte,
#'                         newSubjects = predEnr$newSubjects,
#'                         fitEvt = fitEvt, fitDrp = fitDrp,
#'                         fixedFollowup = FALSE, followupTime = 365,
#'                         pilevel = 0.90, lags = 30, nreps = 500)
#'
#' @export
#'
predictEvent <- function(target_d, df, newSubjects = NULL,
                         fitEvt, fitDrp = NULL,
                         fixedFollowup = FALSE, followupTime = 365,
                         pilevel = 0.90, lags = 30, nreps = 500) {

  erify::check_n(target_d)
  erify::check_class(df, "data.frame")
  if (!is.null(newSubjects)) erify::check_class(newSubjects, "data.frame")
  erify::check_class(fitEvt, "list")
  if (!is.null(fitDrp)) erify::check_class(fitDrp, "list")
  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(lags, zero=TRUE)
  erify::check_n(nreps)


  trialsdt = min(df$RANDDT)
  cutoffdt = df$CUTOFFDT[1]
  d0 = sum(df$event)
  t0 = as.numeric(cutoffdt - trialsdt + 1)
  d1 = target_d - d0  # number of new events

  erify::check_n(d1)

  ongoingSubjects <- df %>%
    filter(.data$event == 0 & .data$dropout == 0) %>%
    mutate(arrivalTime = as.numeric(.data$RANDDT - trialsdt + 1))

  arrivalTimeOngoing <- ongoingSubjects$arrivalTime
  time0 = t0 - arrivalTimeOngoing
  r0 = nrow(ongoingSubjects)

  if (!is.null(newSubjects)) {
    n1 = nrow(newSubjects)/nreps
  } else {
    n1 = 0
  }

  # extract posterior draws of model parameters
  if (tolower(fitEvt$model) == "exponential") {
    theta2 = rnorm(nreps, mean=fitEvt$theta, sd=sqrt(fitEvt$vtheta))
  } else {
    theta2 = mvtnorm::rmvnorm(nreps, mean=fitEvt$theta, sigma=fitEvt$vtheta)
  }

  if (!is.null(fitDrp)) {
    if (tolower(fitDrp$model) == "exponential") {
      theta3 = rnorm(nreps, mean=fitDrp$theta, sd=sqrt(fitDrp$vtheta))
    } else {
      theta3 = mvtnorm::rmvnorm(nreps, mean=fitDrp$theta, sigma=fitDrp$vtheta)
    }
  }


  # output data frame
  m = r0 + n1 # number of ongoing and new subjects
  newEvents = data.frame(matrix(nrow = nreps*m, ncol = 5))
  colnames(newEvents) = c("draw", "arrivalTime", "time", "event", "dropout")

  for (i in 1:nreps) {
    # concatenate arrival time for ongoing and new subjects
    if (n1 > 0) {
      arrivalTimeNew = newSubjects$arrivalTime[newSubjects$draw == i]
      arrivalTime = c(arrivalTimeOngoing, arrivalTimeNew)
    } else {
      arrivalTime = arrivalTimeOngoing
    }

    # draw event time for ongoing and new subjects
    if (tolower(fitEvt$model) == "exponential") {
      survivalTimeOngoing = rexp(r0, rate=exp(theta2[i])) + time0

      if (n1 > 0) {
        survivalTimeNew = rexp(n1, rate=exp(theta2[i]))
        survivalTime = c(survivalTimeOngoing, survivalTimeNew)
      } else {
        survivalTime = survivalTimeOngoing
      }
    } else if (tolower(fitEvt$model) == "weibull") {
      shape = exp(theta2[i,1])
      scale = exp(theta2[i,2])
      survivalTimeOngoing = (rexp(r0)*scale^shape + time0^shape)^(1/shape)

      if (n1 > 0) {
        survivalTimeNew = rweibull(n1, shape, scale)
        survivalTime = c(survivalTimeOngoing, survivalTimeNew)
      } else {
        survivalTime = survivalTimeOngoing
      }
    } else if (tolower(fitEvt$model) == "log-normal") {
      meanlog = theta2[i,1]
      sdlog = exp(theta2[i,2])

      # first draw truncated normal on the log scale
      y = tmvtnsim::rtnorm(mean=rep(meanlog, r0), sd=sdlog,
                           lower=log(time0), upper=rep(Inf,r0))
      survivalTimeOngoing = exp(y)

      if (n1 > 0) {
        survivalTimeNew = rlnorm(n1, meanlog, sdlog)
        survivalTime = c(survivalTimeOngoing, survivalTimeNew)
      } else {
        survivalTime = survivalTimeOngoing
      }
    } else if (tolower(fitEvt$model) == "piecewise exponential") {
      lambda = exp(theta2[i,]) # hazard rates in the intervals
      J = length(lambda) # number of intervals
      u = c(0, fitEvt$knots) # left end points of the intervals

      # partial sums of lambda*interval_width
      psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))

      # find the interval containing time0
      j0 = findInterval(time0, u)
      rhs = psum[j0] + lambda[j0]*(time0 - u[j0]) - log(runif(r0))

      # find the interval containing time
      j1 = findInterval(rhs, psum)
      survivalTimeOngoing = u[j1] + (rhs - psum[j1])/lambda[j1]

      if (n1 > 0) {
        rhs = -log(runif(n1))
        j1 = findInterval(rhs, psum)
        survivalTimeNew = u[j1] + (rhs - psum[j1])/lambda[j1]
        survivalTime = c(survivalTimeOngoing, survivalTimeNew)
      } else {
        survivalTime = survivalTimeOngoing
      }
    } else if (tolower(fitEvt$model) == "model averaging") {
      # distribution function for model averaging of Weibull and log-normal
      pmodavg <- function(t, theta, w1, lower.tail = TRUE) {
        shape = exp(theta[1])
        scale = exp(theta[2])
        meanlog = theta[3]
        sdlog = exp(theta[4])

        p1 = pweibull(pmax(0,t), shape, scale)
        p2 = plnorm(pmax(0,t), meanlog, sdlog)
        p = w1*p1 + (1-w1)*p2

        if (!lower.tail) p = 1 - p
        p
      }

      # inverse transform method to generate event time given surviving time0
      fOngoing_avg <- function(time, theta, w1, time0, u) {
        pmodavg(time, theta, w1, lower.tail = FALSE) -
          u*pmodavg(time0, theta, w1, lower.tail = FALSE)
      }

      u = runif(r0)
      interval = cbind(time0, time0 + 365*10)
      survivalTimeOngoing = rstpm2::vuniroot(fOngoing_avg, interval,
                                             theta = theta2[i,],
                                             fitEvt$w1, time0, u,
                                             extendInt = "yes")$root
      if (n1 > 0) {
        theta = theta2[i,]
        shape = exp(theta[1])
        scale = exp(theta[2])
        meanlog = theta[3]
        sdlog = exp(theta[4])

        # draw component indicator
        w = rbinom(n1, 1, fitEvt$w1)

        # draw from the corresponding component distribution
        survivalTimeNew = rep(NA, n1)
        survivalTimeNew[w==1] = rweibull(sum(w), shape, scale)
        survivalTimeNew[w==0] = rlnorm(sum(1-w), meanlog, sdlog)

        survivalTime = c(survivalTimeOngoing, survivalTimeNew)
      } else {
        survivalTime = survivalTimeOngoing
      }
    }


    # draw dropout time for ongoing and new subjects
    if (!is.null(fitDrp)) {
      if (tolower(fitDrp$model) == "exponential") {
        dropoutTimeOngoing = rexp(r0, rate=exp(theta3[i])) + time0

        if (n1 > 0) {
          dropoutTimeNew = rexp(n1, rate=exp(theta3[i]))
          dropoutTime = c(dropoutTimeOngoing, dropoutTimeNew)
        } else {
          dropoutTime = dropoutTimeOngoing
        }
      } else if (tolower(fitDrp$model) == "weibull") {
        shape = exp(theta3[i,1])
        scale = exp(theta3[i,2])
        dropoutTimeOngoing = (rexp(r0)*scale^shape + time0^shape)^(1/shape)

        if (n1 > 0) {
          dropoutTimeNew = rweibull(n1, shape, scale)
          dropoutTime = c(dropoutTimeOngoing, dropoutTimeNew)
        } else {
          dropoutTime = dropoutTimeOngoing
        }
      } else if (tolower(fitDrp$model) == "log-normal") {
        meanlog = theta3[i,1]
        sdlog = exp(theta3[i,2])

        # first draw truncated normal on the log scale
        y = tmvtnsim::rtnorm(mean = rep(meanlog, r0), sd = sdlog,
                             lower = log(time0), upper = rep(Inf, r0))
        dropoutTimeOngoing = exp(y)

        if (n1 > 0) {
          dropoutTimeNew = rlnorm(n1, meanlog, sdlog)
          dropoutTime = c(dropoutTimeOngoing, dropoutTimeNew)
        } else {
          dropoutTime = dropoutTimeOngoing
        }
      }

    }


    # observed survival time and event indicator
    if (!fixedFollowup) {
      if (!is.null(fitDrp)) {
        time = pmin(survivalTime, dropoutTime)
        event = 1*(time == survivalTime)
        dropout = 1*(time == dropoutTime)
      } else {
        time = survivalTime
        event = 1
        dropout = 0
      }
    } else {
      if (!is.null(fitDrp)) {
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
    index = (i-1)*m + (1:m)
    newEvents[index, "draw"] = i
    newEvents[index, "arrivalTime"] = arrivalTime
    newEvents[index, "time"] = time
    newEvents[index, "event"] = event
    newEvents[index, "dropout"] = dropout
  }

  # calculate total time since trial start
  newEvents <- newEvents %>%
    mutate(totalTime = .data$arrivalTime + .data$time)

  t1 = min(max(newEvents$totalTime[newEvents$event == 1]), t0 + 365*4)

  # lower and upper percentages
  plower = (1 - pilevel)/2
  pupper = 1 - plower

  # A general quantile method if there are data sets not reaching target_d
  # Find t such that Sum[I(D_i(t) < target_d), {i, 1, nreps}] / nreps = q
  # This works because {D_i(t) < target_d} = {T_i(target_d) > t},
  # where D_i(t) is the cumulative number of events at time t, and
  # T_i(target_d) is the time to reach target_d for data set i
  # empirical survival function for T(target_d) based on the above method
  f <- function(t, target_d, d0, newEvents) {
    sumdata <- newEvents %>%
      group_by(.data$draw) %>%
      summarize(d = sum(.data$totalTime <= t & .data$event == 1) + d0)
    mean(sumdata$d < target_d)
  }

  # obtain the quantiles
  q = 1 - c(0.5, plower, pupper)
  pred2dy = rep(NA, length(q))
  for (j in 1:length(q)) {
    # check if the quantile can be estimated from observed data
    if (f(t1, target_d, d0, newEvents) <= q[j]) {
      pred2dy[j] = uniroot(function(x) f(x, target_d, d0, newEvents) - q[j],
                           c(t0, t1))$root
      pred2dy[j] = ceiling(pred2dy[j])
    }
  }
  pred2dt = as.Date(pred2dy - 1, origin = trialsdt)


  if (!is.na(pred2dy[1])) {
    r2month = round(pred2dy[1]/30.4375, digits = 1)
  } else {
    r2month = NA
  }


  # set up time points for plotting
  t = c(seq(t0, t1, 30), t1)


  # observed number of events before data cut
  dfa <- df %>%
    filter(.data$event == 1 | .data$dropout == 1) %>%
    mutate(arrivalTime = as.numeric(.data$RANDDT - trialsdt + 1)) %>%
    mutate(totalTime = .data$arrivalTime + .data$time) %>%
    arrange(.data$totalTime) %>%
    mutate(time = .data$totalTime,
           d = cumsum(.data$event),
           lower = NA,
           upper = NA) %>%
    select(.data$time, .data$d, .data$lower, .data$upper)

  # predicted number of events after data cut
  dfb = data.frame(matrix(nrow=length(t), ncol=4))
  colnames(dfb) = c('time', 'd', 'lower', 'upper')
  for (i in 1:length(t)) {
    # number of events after data cut in each simulated data set
    sumdata <- newEvents %>%
      group_by(.data$draw) %>%
      summarize(d = sum(.data$totalTime <= t[i] & .data$event == 1) + d0)

    # summary across simulated data sets
    dfb[i, 'time'] = t[i]
    dfb[i, c('d', 'lower', 'upper')] =
      quantile(sumdata$d, probs = c(0.5, plower, pupper))
  }


  # add time zero and concatenate subjects enrolled before and after data cut
  df2 <- tibble(time = 0, d = 0, lower = NA, upper = NA) %>%
    bind_rows(dfa) %>%
    bind_rows(dfb) %>%
    mutate(date = as.Date(.data$time - 1, origin = trialsdt)) %>%
    mutate(year = format(.data$date, format = "%Y"))


  # separate data into observed and predicted
  df2a <- df2 %>% filter(is.na(.data$lower))
  df2b <- df2 %>% filter(!is.na(.data$lower))


  n_months = lubridate::interval(min(df2$date), max(df2$date)) %/% months(1)
  bw = fbw(n_months)

  g2 <- flabel(df2, trialsdt)

  # plot the enrollment and time to event data with month as x-axis label
  # generate the plot
  g1 <- ggplot() +
    geom_ribbon(data=df2b,
                aes(x=.data$date, ymin=.data$lower, ymax=.data$upper),
                alpha=0.5, fill="lightblue") +
    geom_step(data=df2a, aes(x=.data$date, y=.data$d), color="black") +
    geom_line(data=df2b, aes(x=.data$date, y=.data$d), color="blue") +
    geom_vline(xintercept = cutoffdt, linetype = 2) +
    geom_hline(yintercept = target_d, linetype = 2) +
    scale_x_date(name = NULL,
                 labels = scales::date_format("%b"),
                 breaks = scales::breaks_width(bw),
                 minor_breaks = NULL,
                 expand = c(0.01, 0.01)) +
    labs(y = "Events", title = "Predicted events") +
    theme_bw()


  # stack them together
  p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
  print(p1)

  list(remEvtMonth = r2month, predEvtDT = pred2dt, pilevel = pilevel,
       newEvents = newEvents)
}
