#' @title Predict event
#' @description Utilizes pre-fitted time-to-event and time-to-dropout models
#'   to generate event and dropout times for ongoing subjects
#'   and new subjects. It also provides a
#'   prediction interval for the expected time to reach the target
#'   number of events.
#'
#' @param df The subject-level enrollment and event data,
#'   including \code{randdt}, \code{cutoffdt},
#'   \code{time}, \code{event}, and \code{dropout}. By default, it
#'   is set to \code{NULL} for event prediction at the design stage.
#' @param target_d The target number of events to reach in the study.
#' @param newSubjects The enrollment data for new subjects. By default,
#'   it is set to \code{NULL}, indicating the completion of
#'   subject enrollment.
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
#' @param nreps The number of replications for simulation. By default,
#'   it is set to 500.
#' @param showplot A Boolean variable to control whether or not
#'   the prediction plot is displayed. By default, it is set to
#'   \code{TRUE}.
#'
#' @details
#' To ensure successful event prediction at the design stage, it is
#' important to provide the \code{newSubjects} data set.
#'
#' To specify the event model used during the design-stage event
#' prediction, the \code{event_fit} list should include the event model,
#' \code{model}, the number of treatment groups,
#' \code{ngroups}, the randomization probabilities for each group,
#' \code{prob}, the model parameters, \code{theta},
#' and the covariance matrix, \code{vtheta}, both of which
#' have \code{ngroups} blocks with the \code{j}-th block
#' specifying the prior distribution of model
#' parameters for the \code{j}-th treatment group. For the
#' piecewise exponential event model, \code{knots} should also
#' be included to indicate the location of inner knots.
#' It should be noted that the model averaging option is not
#' appropriate for use during the design stage.
#'
#' To specify the dropout model used during the design stage
#' event prediction, the \code{dropout_fit} list should include
#' the dropout model, \code{model}, the number of treatment groups,
#' \code{ngroups}, the randomization probabilities for each
#' group, \code{prob}, the model parameters, \code{theta}, and
#' the covariance matrix ,\code{vtheta}, both of which have
#' \code{ngroups} blocks with the \code{j}-th block specifying
#' the prior distribution of model parameters for the
#' \code{j}-th treatment group.
#'
#' @return A list of prediction results which includes important
#' information such as the median, lower and upper percentiles for
#' the estimated time and date to reach the target number of events,
#' as well as simulated event data for both ongoing and new subjects.
#' Additionally, the data for the prediction plot is also included
#' within this list.
#'
#' @examples
#'
#' # Example 1: Event prediction at analysis stage after enrollment ends
#'
#' event_fit <- fitEvent(df = observedData,
#'                       event_model = "piecewise exponential", npieces = 3)
#'
#' dropout_fit <- fitDropout(df = observedData,
#'                           dropout_model = "exponential")
#'
#' event_pred <- predictEvent(df = observedData, target_d = 200,
#'                            event_fit = event_fit,
#'                            dropout_fit = dropout_fit,
#'                            pilevel = 0.90, nreps = 500)
#'
#'
#' # Example 2: Event prediction at analysis stage before enrollment ends
#'
#' observed <- summarizeObserved(df = observedData,
#'                               to_predict = "enrollment and event")
#'
#' enroll_fit <- fitEnrollment(df = observed$adsl, enroll_model = "b-spline",
#'                             nknots = 1)
#'
#' enroll_pred <- predictEnrollment(df = observed$adsl, target_n = 480,
#'                                  enroll_fit = enroll_fit,
#'                                  lags = 30, pilevel = 0.90, nreps = 500)
#'
#' event_fit <- fitEvent(df = observed$adtte,
#'                       event_model = "piecewise exponential", npieces = 3)
#'
#' dropout_fit <- fitDropout(df = observed$adtte,
#'                           dropout_model = "exponential")
#'
#' event_pred <- predictEvent(df = observed$adtte, target_d = 200,
#'                            newSubjects = enroll_pred$newSubjects,
#'                            event_fit = event_fit,
#'                            dropout_fit = dropout_fit,
#'                            pilevel = 0.90, nreps = 500)
#'
#' # Example 3: Event prediction at design stage
#'
#' parameter_enroll_model <- list(
#'   model = "piecewise poisson",
#'   accrualTime = seq(0, 8)*30.4375,
#'   accrualIntensity = 26/9*seq(1, 9)/30.4375)
#'
#' parameter_event_model <- list(
#'   model = "piecewise exponential",
#'   ngroups = 2,
#'   prob = c(0.5, 0.5),
#'   theta = log(c(0.0533, 0.0309, 0.0533, 0.0533)/30.4375),
#'   vtheta = diag(4)*1e-8,
#'   knots = 6*30.4375)
#'
#' parameter_dropout_model <- list(
#'   model = "exponential",
#'   ngroups = 2,
#'   prob = c(0.5, 0.5),
#'   theta = log(rep(-log(1-0.05)/12, 2)/30.4375),
#'   vtheta = diag(2)*1e-8)
#'
#' enroll_pred <- predictEnrollment(
#'   enroll_fit = parameter_enroll_model,
#'   target_n = 480,
#'   pilevel = 0.90, nreps = 500)
#'
#' event_pred <- predictEvent(
#'   target_d = 200,
#'   newSubjects = enroll_pred$newSubjects,
#'   event_fit = parameter_event_model,
#'   dropout_fit = parameter_dropout_model,
#'   pilevel = 0.90, nreps = 500)
#'
#' @export
#'
predictEvent <- function(df = NULL, target_d, newSubjects = NULL,
                         event_fit, dropout_fit = NULL,
                         fixedFollowup = FALSE, followupTime = 365,
                         pilevel = 0.90, nreps = 500, showplot = TRUE) {

  if (!is.null(df)) erify::check_class(df, "data.frame")
  erify::check_n(target_d)
  if (!is.null(newSubjects)) erify::check_class(newSubjects, "data.frame")
  if (is.null(df) & is.null(newSubjects)) {
    stop("At least one of df and newSubjects must be specified.")
  }
  erify::check_class(event_fit, "list")
  if (!is.null(dropout_fit)) erify::check_class(dropout_fit, "list")
  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(nreps)
  erify::check_bool(showplot)

  if (!is.null(df)) {
    names(df) <- tolower(names(df))
    trialsdt = min(df$randdt)
    cutoffdt = df$cutoffdt[1]
    d0 = sum(df$event)
    t0 = as.numeric(cutoffdt - trialsdt + 1)
  } else {
    d0 = 0
    t0 = 1
  }
  d1 = target_d - d0  # number of new events

  erify::check_n(d1)

  if (!is.null(df)) {
    ongoingSubjects <- df %>%
      filter(.data$event == 0 & .data$dropout == 0) %>%
      mutate(arrivalTime = as.numeric(.data$randdt - trialsdt + 1))

    arrivalTimeOngoing <- ongoingSubjects$arrivalTime
    time0 = t0 - arrivalTimeOngoing
    r0 = nrow(ongoingSubjects)
  } else {
    time0 = t0
    r0 = 0
  }

  if (!is.null(newSubjects)) {
    n1 = nrow(newSubjects)/nreps
  } else {
    n1 = 0
  }



  # extract posterior draws of model parameters
  if (!is.null(df)) { # analysis stage event prediction
    if (tolower(event_fit$model) == "exponential") {
      theta2 = rnorm(nreps, mean = event_fit$theta,
                     sd = sqrt(event_fit$vtheta))
    } else {
      theta2 = mvtnorm::rmvnorm(nreps, mean = event_fit$theta,
                                sigma = event_fit$vtheta)
    }

    if (!is.null(dropout_fit)) {
      if (tolower(dropout_fit$model) == "exponential") {
        theta3 = rnorm(nreps, mean=dropout_fit$theta,
                       sd = sqrt(dropout_fit$vtheta))
      } else {
        theta3 = mvtnorm::rmvnorm(nreps, mean = dropout_fit$theta,
                                  sigma = dropout_fit$vtheta)
      }
    }
  } else { # design stage event prediction
    if (tolower(event_fit$model) == "exponential") {
      k = event_fit$ngroups # number of treatment groups
      theta2 <- list() # create an empty list
      # posterior draws of model parameters for each treatment group
      for (j in 1:k) {
        theta2[[j]] <- rnorm(nreps, mean = event_fit$theta[j],
                             sd = sqrt(event_fit$vtheta[j]))
      }
    } else {
      k = event_fit$ngroups # number of treatment groups
      p = length(event_fit$theta)/k # number of parameters for each group
      theta2 <- list() # create an empty list
      # posterior draws of model parameters for each treatment group
      for (j in 1:k) {
        jp = (j-1)*p + (1:p)
        theta2[[j]] = mvtnorm::rmvnorm(nreps, mean = event_fit$theta[jp],
                                       sigma = event_fit$vtheta[jp,jp])
      }
    }

    if (!is.null(dropout_fit)) {
      if (tolower(dropout_fit$model) == "exponential") {
        k = dropout_fit$ngroups # number of treatment groups
        theta3 <- list() # create an empty list
        # posterior draws of model parameters for each treatment group
        for (j in 1:k) {
          theta3[[j]] <- rnorm(nreps, mean = dropout_fit$theta[j],
                               sd = sqrt(dropout_fit$vtheta[j]))
        }
      } else {
        k = dropout_fit$ngroups # number of treatment groups
        p = length(dropout_fit$theta)/k # number of parameters for each group
        theta3 <- list() # create an empty list
        # posterior draws of model parameters for each treatment group
        for (j in 1:k) {
          jp = (j-1)*p + (1:p)
          theta3[[j]] = mvtnorm::rmvnorm(nreps, mean = dropout_fit$theta[jp],
                                         sigma = dropout_fit$vtheta[jp,jp])
        }
      }
    }
  }


  # define distribution function for model averaging
  if (tolower(event_fit$model) == "model averaging") {
    if (is.null(df)) {
      stop("Model averaging is only used at the analysis stage.")
    }

    # distribution function for model averaging of Weibull and log-normal
    pmodavg <- function(t, theta, w1, lower.tail = TRUE, log.p = FALSE) {
      shape = exp(theta[1])
      scale = exp(theta[2])
      meanlog = theta[3]
      sdlog = exp(theta[4])

      p1 = pweibull(pmax(0,t), shape, scale)
      p2 = plnorm(pmax(0,t), meanlog, sdlog)
      p = w1*p1 + (1-w1)*p2

      if (!lower.tail) p = 1 - p
      if (log.p) p = log(p)
      p
    }

    # inverse transform method to generate event time given surviving time0
    fOngoing_avg <- function(time, theta, w1, time0, u) {
      pmodavg(time, theta, w1, lower.tail = FALSE) -
        u*pmodavg(time0, theta, w1, lower.tail = FALSE)
    }
  }


  # generate the event and dropout times
  if (!is.null(df)) { # analysis stage event prediction
    # output data frame
    m = r0 + n1 # number of ongoing and new subjects
    newEvents = data.frame(matrix(nrow = nreps*m, ncol = 5))
    colnames(newEvents) = c("draw", "arrivalTime", "time", "event",
                            "dropout")

    for (i in 1:nreps) {
      # concatenate arrival time for ongoing and new subjects
      if (n1 > 0) {
        arrivalTimeNew = newSubjects$arrivalTime[newSubjects$draw == i]
        arrivalTime = c(arrivalTimeOngoing, arrivalTimeNew)
      } else {
        arrivalTime = arrivalTimeOngoing
      }

      # draw event time for ongoing and new subjects
      if (tolower(event_fit$model) == "exponential") {
        survivalTimeOngoing = rexp(r0, rate = exp(theta2[i])) + time0

        if (n1 > 0) {
          survivalTimeNew = rexp(n1, rate = exp(theta2[i]))
          survivalTime = c(survivalTimeOngoing, survivalTimeNew)
        } else {
          survivalTime = survivalTimeOngoing
        }
      } else if (tolower(event_fit$model) == "weibull") {
        shape = exp(theta2[i,1])
        scale = exp(theta2[i,2])
        survivalTimeOngoing = (rexp(r0)*scale^shape + time0^shape)^(1/shape)

        if (n1 > 0) {
          survivalTimeNew = rweibull(n1, shape, scale)
          survivalTime = c(survivalTimeOngoing, survivalTimeNew)
        } else {
          survivalTime = survivalTimeOngoing
        }
      } else if (tolower(event_fit$model) == "log-normal") {
        meanlog = theta2[i,1]
        sdlog = exp(theta2[i,2])

        # first draw truncated normal on the log scale
        y = tmvtnsim::rtnorm(mean = rep(meanlog, r0), sd = sdlog,
                             lower = log(time0), upper = rep(Inf, r0))
        survivalTimeOngoing = exp(y)

        if (n1 > 0) {
          survivalTimeNew = rlnorm(n1, meanlog, sdlog)
          survivalTime = c(survivalTimeOngoing, survivalTimeNew)
        } else {
          survivalTime = survivalTimeOngoing
        }
      } else if (tolower(event_fit$model) == "piecewise exponential") {
        lambda = exp(theta2[i,]) # hazard rates in the intervals
        J = length(lambda) # number of intervals
        u = c(0, event_fit$knots) # left end points of the intervals

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
      } else if (tolower(event_fit$model) == "model averaging") {
        u = runif(r0)
        interval = cbind(time0, time0 + 365*10)
        survivalTimeOngoing = rstpm2::vuniroot(
          fOngoing_avg, interval, theta = theta2[i,],
          w1 = event_fit$w1, time0, u, extendInt = "yes")$root

        if (n1 > 0) {
          theta = theta2[i,]
          shape = exp(theta[1])
          scale = exp(theta[2])
          meanlog = theta[3]
          sdlog = exp(theta[4])

          # draw component indicator
          w = rbinom(n1, size = 1, prob = event_fit$w1)

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
      if (!is.null(dropout_fit)) {
        if (tolower(dropout_fit$model) == "exponential") {
          dropoutTimeOngoing = rexp(r0, rate = exp(theta3[i])) + time0

          if (n1 > 0) {
            dropoutTimeNew = rexp(n1, rate = exp(theta3[i]))
            dropoutTime = c(dropoutTimeOngoing, dropoutTimeNew)
          } else {
            dropoutTime = dropoutTimeOngoing
          }
        } else if (tolower(dropout_fit$model) == "weibull") {
          shape = exp(theta3[i,1])
          scale = exp(theta3[i,2])
          dropoutTimeOngoing = (rexp(r0)*scale^shape + time0^shape)^(1/shape)

          if (n1 > 0) {
            dropoutTimeNew = rweibull(n1, shape, scale)
            dropoutTime = c(dropoutTimeOngoing, dropoutTimeNew)
          } else {
            dropoutTime = dropoutTimeOngoing
          }
        } else if (tolower(dropout_fit$model) == "log-normal") {
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
      index = (i-1)*m + (1:m)
      newEvents[index, "draw"] = i
      newEvents[index, "arrivalTime"] = arrivalTime
      newEvents[index, "time"] = time
      newEvents[index, "event"] = event
      newEvents[index, "dropout"] = dropout
    }

  } else { # design stage event prediction
    # output data frame
    m = r0 + n1 # number of ongoing and new subjects
    newEvents = data.frame(matrix(nrow = nreps*m, ncol = 6))
    colnames(newEvents) = c("draw", "arrivalTime", "treatment", "time",
                            "event", "dropout")

    k = event_fit$ngroups # number of treatment groups
    prob = event_fit$prob # randomization probabilities

    for (i in 1:nreps) {
      arrivalTime = newSubjects$arrivalTime[newSubjects$draw == i]

      # draw treatment for each subject
      w = rmultinom(n1, size = 1, prob = prob)
      treatment = as.numeric((1:k) %*% w)

      # draw event time for new subjects
      survivalTime = rep(NA, n1)
      if (tolower(event_fit$model) == "exponential") {
        for (j in 1:k) {
          cols = which(w[j,] == 1)
          if (length(cols) > 0) {
            rate = exp(theta2[[j]][i])
            survivalTime[cols] = rexp(length(cols), rate)
          }
        }
      } else if (tolower(event_fit$model) == "weibull") {
        for (j in 1:k) {
          cols = which(w[j,] == 1)
          if (length(cols) > 0) {
            shape = exp(theta2[[j]][i,1])
            scale = exp(theta2[[j]][i,2])
            survivalTime[cols] = rweibull(length(cols), shape, scale)
          }
        }
      } else if (tolower(event_fit$model) == "log-normal") {
        for (j in 1:k) {
          cols = which(w[j,] == 1)
          if (length(cols) > 0) {
            meanlog = theta2[[j]][i,1]
            sdlog = exp(theta2[[j]][i,2])
            survivalTime[cols] = rlnorm(length(cols), meanlog, sdlog)
          }
        }
      } else if (tolower(event_fit$model) == "piecewise exponential") {
        u = c(0, event_fit$knots) # left end points of the intervals
        J = length(u) # number of intervals
        for (j in 1:k) {
          cols = which(w[j,] == 1)
          if (length(cols) > 0) {
            lambda = exp(theta2[[j]][i,]) # hazard rates in the intervals
            # partial sums of lambda*interval_width
            psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))
            rhs = -log(runif(length(cols))) # cumulative hazard
            j1 = findInterval(rhs, psum)
            survivalTime[cols] = u[j1] + (rhs - psum[j1])/lambda[j1]
          }
        }
      }


      # draw dropout time for ongoing and new subjects
      if (!is.null(dropout_fit)) {
        dropoutTime = rep(NA, n1)

        if (tolower(dropout_fit$model) == "exponential") {
          for (j in 1:k) {
            cols = which(w[j,] == 1)
            if (length(cols) > 0) {
              rate = exp(theta3[[j]][i])
              dropoutTime[cols] = rexp(length(cols), rate)
            }
          }
        } else if (tolower(dropout_fit$model) == "weibull") {
          for (j in 1:k) {
            cols = which(w[j,] == 1)
            if (length(cols) > 0) {
              shape = exp(theta3[[j]][i,1])
              scale = exp(theta3[[j]][i,2])
              dropoutTime[cols] = rweibull(length(cols), shape, scale)
            }
          }
        } else if (tolower(dropout_fit$model) == "log-normal") {
          for (j in 1:k) {
            cols = which(w[j,] == 1)
            if (length(cols) > 0) {
              meanlog = theta3[[j]][i,1]
              sdlog = exp(theta3[[j]][i,2])
              dropoutTime[cols] = rlnorm(length(cols), meanlog, sdlog)
            }
          }
        }
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
      index = (i-1)*m + (1:m)
      newEvents[index, "draw"] = i
      newEvents[index, "arrivalTime"] = arrivalTime
      newEvents[index, "treatment"] = treatment
      newEvents[index, "time"] = time
      newEvents[index, "event"] = event
      newEvents[index, "dropout"] = dropout
    }
  }


  # calculate total time since trial start
  newEvents <- newEvents %>%
    mutate(totalTime = .data$arrivalTime + .data$time)

  t1 = min(max(newEvents$totalTime[newEvents$event == 1]), t0 + 365*4)

  # lower and upper percentages
  plower = (1 - pilevel)/2
  pupper = 1 - plower

  # A general quantile method if there are data sets not reaching target_d
  # Find t such that sum(I{D_i(t) < target_d}, {i, 1, nreps}) / nreps = q.
  # This works because {D_i(t) < target_d} = {T_i(target_d) > t},
  # where D_i(t) is the cumulative number of events at time t, and
  # T_i(target_d) is the time to reach target_d for data set i.
  sdf <- function(t, target_d, d0, newEvents) {
    sumdata <- newEvents %>%
      group_by(.data$draw) %>%
      summarize(n = sum(.data$totalTime <= t & .data$event == 1) + d0)
    mean(sumdata$n < target_d)
  }

  # obtain the quantiles
  q = 1 - c(0.5, plower, pupper)
  pred2dy = rep(NA, length(q))
  for (j in 1:length(q)) {
    # check if the quantile can be estimated from observed data
    if (sdf(t1, target_d, d0, newEvents) <= q[j]) {
      pred2dy[j] = uniroot(function(x) sdf(x, target_d, d0, newEvents) - q[j],
                           c(t0, t1))$root
      pred2dy[j] = ceiling(pred2dy[j])
    }
  }


  # set up time points for plotting
  t = c(seq(t0, t1, 30), t1)

  # predicted number of events after data cut
  dfb = data.frame(matrix(nrow=length(t), ncol=4))
  colnames(dfb) = c('time', 'n', 'lower', 'upper')
  for (i in 1:length(t)) {
    # number of events after data cut in each simulated data set
    sumdata <- newEvents %>%
      group_by(.data$draw) %>%
      summarize(n = sum(.data$totalTime <= t[i] & .data$event == 1) + d0)

    # summary across simulated data sets
    dfb[i, 'time'] = t[i]
    dfb[i, c('n', 'lower', 'upper')] =
      quantile(sumdata$n, probs = c(0.5, plower, pupper))
  }


  if (!is.null(df)) {
    pred2dt <- as.Date(pred2dy - 1, origin = trialsdt)

    # observed number of events before data cut
    dfa <- df %>%
      filter(.data$event == 1 | .data$dropout == 1) %>%
      mutate(arrivalTime = as.numeric(.data$randdt - trialsdt + 1)) %>%
      mutate(totalTime = .data$arrivalTime + .data$time) %>%
      arrange(.data$totalTime) %>%
      mutate(time = .data$totalTime,
             n = cumsum(.data$event),
             lower = NA,
             upper = NA) %>%
      select(.data$time, .data$n, .data$lower, .data$upper)


    # add time zero and concatenate events before and after data cut
    dfs <- tibble(time = 0, n = 0, lower = NA, upper = NA) %>%
      bind_rows(dfa) %>%
      bind_rows(dfb) %>%
      mutate(date = as.Date(.data$time - 1, origin = trialsdt)) %>%
      mutate(year = format(.data$date, format = "%Y"))

    if (showplot) {
      # separate data into observed and predicted
      dfa <- dfs %>% filter(is.na(.data$lower))
      dfb <- dfs %>% filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      g2 <- flabel(dfs, trialsdt)

      # plot the enrollment and time to event data with month as x-axis label
      # generate the plot
      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$date, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_step(data=dfa, aes(x=.data$date, y=.data$n), color="black") +
        geom_line(data=dfb, aes(x=.data$date, y=.data$n), color="blue") +
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
    }

    list(predEventDay = pred2dy, predEventDate = pred2dt, pilevel = pilevel,
         newEvents = newEvents, plotdata = dfs)
  } else {
    if (showplot) {
      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$time, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_line(data=dfb, aes(x=.data$time, y=.data$n), color="blue") +
        geom_hline(yintercept = target_d, linetype = 2) +
        scale_x_continuous(name = "Days since randomization",
                           expand = c(0.01, 0.01)) +
        labs(y = "Events", title = "Predicted events") +
        theme_bw()

      print(g1)
    }

    list(predEventDay = pred2dy, pilevel = pilevel, newEvents = newEvents,
         plotdata = dfb)
  }
}
