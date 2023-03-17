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
#' @param newSubjects The enrollment data for new subjects including
#'   \code{draw} and \code{arrivalTime}. By default,
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
#' @param nyears The number of years after the data cut for prediction.
#'   By default, it is set to 4.
#' @param nreps The number of replications for simulation. By default,
#'   it is set to 500.
#' @param showDropout A Boolean variable to control whether or not to
#'   show the number of dropouts. By default, it is set to
#'   \code{FALSE}.
#' @param showOngoing A Boolean variable to control whether or not to
#'   show the number of ongoing subjects. By default, it is set to
#'   \code{FALSE}.
#' @param showplot A Boolean variable to control whether or not to
#'   show the prediction plot. By default, it is set to \code{TRUE}.
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
#' piecewise exponential event model, \code{piecewiseSurvivalTime}
#' should also be included to indicate the location of knots.
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
#' event_pred <- predictEvent(df = interimData2, target_d = 244,
#'                            event_fit = event_fit$event_fit,
#'                            dropout_fit = dropout_fit$dropout_fit,
#'                            pilevel = 0.90, nreps = 200)
#'
#' @export
#'
predictEvent <- function(df = NULL, target_d, newSubjects = NULL,
                         event_fit, dropout_fit = NULL,
                         fixedFollowup = FALSE, followupTime = 365,
                         pilevel = 0.90, nyears = 4, nreps = 500,
                         showDropout = FALSE, showOngoing = FALSE,
                         showplot = TRUE) {

  if (!is.null(df)) erify::check_class(df, "data.frame")
  erify::check_n(target_d)
  if (!is.null(newSubjects)) erify::check_class(newSubjects, "data.frame")
  if (is.null(df) & is.null(newSubjects)) {
    stop("At least one of df and newSubjects must be specified.")
  }

  erify::check_class(event_fit, "list")
  erify::check_content(tolower(event_fit$model),
                       c("exponential", "weibull", "log-normal",
                         "piecewise exponential", "model averaging"))

  if (!is.null(dropout_fit)) {
    erify::check_class(dropout_fit, "list")
    erify::check_content(tolower(dropout_fit$model),
                         c("exponential", "weibull", "log-normal"))
  }

  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_positive(nyears)
  erify::check_n(nreps)
  erify::check_bool(showDropout)
  erify::check_bool(showOngoing)
  erify::check_bool(showplot)

  if (!is.null(df)) {
    df <- dplyr::as_tibble(df)
    names(df) <- tolower(names(df))
    trialsdt = min(df$randdt)
    cutoffdt = df$cutoffdt[1]
    d0 = sum(df$event)
    c0 = sum(df$dropout)
    t0 = as.numeric(cutoffdt - trialsdt + 1)
  } else {
    d0 = 0
    c0 = 0
    t0 = 1
  }
  d1 = target_d - d0  # number of new events

  erify::check_n(d1)

  if (!is.null(df)) {
    n0 = nrow(df)
    ongoingSubjects <- df %>%
      dplyr::filter(.data$event == 0 & .data$dropout == 0) %>%
      dplyr::mutate(arrivalTime = as.numeric(.data$randdt - trialsdt + 1))

    arrivalTimeOngoing <- ongoingSubjects$arrivalTime
    time0 = t0 - arrivalTimeOngoing
    r0 = nrow(ongoingSubjects)
  } else {
    n0 = 0
    time0 = t0
    r0 = 0
  }


  if (!is.null(newSubjects)) {
    n1 = nrow(newSubjects)/nreps

    # lower and upper percentages
    plower = (1 - pilevel)/2
    pupper = 1 - plower

    new1 <- newSubjects %>%
      dplyr::group_by(.data$draw) %>%
      dplyr::filter(dplyr::row_number() == dplyr::n())

    pred1dy <- ceiling(quantile(new1$arrivalTime, c(0.5, plower, pupper)))

    t1 = pred1dy[3] + 30 # extend time to 30 days after

    # future time points at which to predict number of subjects
    t = sort(unique(c(seq(t0, t1, 30), t1, pred1dy)))

    # predicted number of subjects enrolled after data cut
    dfb <- dplyr::tibble(t = t) %>%
      dplyr::cross_join(newSubjects) %>%
      dplyr::group_by(.data$t, .data$draw) %>%
      dplyr::summarise(nenrolled = sum(.data$arrivalTime <= .data$t) + n0,
                       .groups = "drop_last") %>%
      dplyr::summarise(n = quantile(.data$nenrolled, probs = 0.5),
                       lower = quantile(.data$nenrolled, probs = plower),
                       upper = quantile(.data$nenrolled, probs = pupper))

    if (!is.null(df)) {
      # arrival time for subjects already enrolled before data cut
      dfa <- df %>%
        dplyr::arrange(.data$randdt) %>%
        dplyr::mutate(t = as.numeric(.data$randdt - trialsdt + 1),
                      n = dplyr::row_number()) %>%
        dplyr::mutate(lower = NA, upper = NA) %>%
        dplyr::select(.data$t, .data$n, .data$lower, .data$upper)

      # concatenate subjects enrolled before and after data cut
      enroll_pred_df <- dfa %>%
        dplyr::bind_rows(dfb) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(year = format(.data$date, format = "%Y"))
    } else {
      enroll_pred_df <- dfb
    }
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
    newEvents = dplyr::as_tibble(matrix(
      nrow = nreps*m, ncol = 5,
      dimnames = list(NULL, c("draw", "arrivalTime", "time", "event",
                              "dropout"))))

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
        u = event_fit$piecewiseSurvivalTime # left end points of the intervals

        # partial sums of lambda*interval_width
        if (J > 1) {
          psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))
        } else {
          psum = 0
        }

        # find the interval containing time0
        j0 = findInterval(time0, u)
        rhs = psum[j0] + lambda[j0]*(time0 - u[j0]) + rexp(r0)

        # find the interval containing time
        j1 = findInterval(rhs, psum)
        survivalTimeOngoing = u[j1] + (rhs - psum[j1])/lambda[j1]

        if (n1 > 0) {
          rhs = rexp(n1)
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
    newEvents = dplyr::as_tibble(matrix(
      nrow = nreps*m, ncol = 6,
      dimnames = list(NULL, c("draw", "arrivalTime", "treatment", "time",
                              "event", "dropout"))))

    k = event_fit$ngroups # number of treatment groups
    prob = event_fit$prob # randomization probabilities

    for (i in 1:nreps) {
      arrivalTime = newSubjects$arrivalTime[newSubjects$draw == i]

      # draw treatment for each subject
      w = rmultinom(n1, size = 1, prob = prob)  # k x n1 matrix
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
            if (J > 1) {
              psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))
            } else {
              psum = 0
            }
            rhs = rexp(length(cols)) # cumulative hazard
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
    dplyr::mutate(totalTime = .data$arrivalTime + .data$time)

  # upper bound for time axis for plotting
  t1 = min(max(newEvents$totalTime[newEvents$event == 1]), t0 + 365*nyears)

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
      dplyr::group_by(.data$draw) %>%
      dplyr::summarize(n = sum(.data$totalTime <= t & .data$event == 1) + d0)
    mean(sumdata$n < target_d)
  }

  # obtain the quantiles
  q = 1 - c(0.5, plower, pupper)
  pred2dy = rep(NA, length(q))
  for (j in 1:length(q)) {
    # check if the quantile can be estimated from observed data
    if (sdf(t1, target_d, d0, newEvents) <= q[j]) {
      pred2dy[j] = uniroot(function(x) sdf(x, target_d, d0, newEvents) - q[j],
                           c(t0, t1), tol = 1)$root
      pred2dy[j] = ceiling(pred2dy[j])
    }
  }


  # set up time points for plotting
  t = unique(c(seq(t0, t1, 30), t1))

  # predicted number of events after data cut
  dfb = dplyr::tibble(t = t) %>%
    dplyr::cross_join(newEvents) %>%
    dplyr::group_by(.data$t, .data$draw) %>%
    dplyr::summarise(nevents = sum(.data$totalTime <= .data$t &
                                     .data$event == 1) + d0,
                     .groups = "drop_last") %>%
    dplyr::summarise(n = quantile(.data$nevents, probs = 0.5),
                     lower = quantile(.data$nevents, probs = plower),
                     upper = quantile(.data$nevents, probs = pupper))


  if (!is.null(df)) {
    pred2dt <- as.Date(pred2dy - 1, origin = trialsdt)

    # observed number of events before data cut
    dfa <- df %>%
      dplyr::filter(.data$event == 1 | .data$dropout == 1) %>%
      dplyr::mutate(arrivalTime = as.numeric(.data$randdt - trialsdt + 1)) %>%
      dplyr::mutate(totalTime = .data$arrivalTime + .data$time) %>%
      dplyr::arrange(.data$totalTime) %>%
      dplyr::mutate(t = .data$totalTime,
                    n = cumsum(.data$event),
                    lower = NA,
                    upper = NA) %>%
      dplyr::select(.data$t, .data$n, .data$lower, .data$upper)


    # add time zero and concatenate events before and after data cut
    dfs <- dplyr::tibble(t = 0, n = 0, lower = NA, upper = NA) %>%
      dplyr::bind_rows(dfa) %>%
      dplyr::bind_rows(dfb) %>%
      dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
      dplyr::mutate(year = format(.data$date, format = "%Y"))


    if (is.null(newSubjects)) { # plot predicted events only
      # separate data into observed and predicted
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      # plot the enrollment and time to event data with month as x-axis label
      # generate the plot
      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$date,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_step(data=dfa, ggplot2::aes(x=.data$date, y=.data$n),
                           color="black") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$date, y=.data$n),
                           color="blue") +
        ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_date(name = NULL,
                              labels = scales::date_format("%b"),
                              breaks = scales::breaks_width(bw),
                              minor_breaks = NULL,
                              expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Events", title = "Predicted events") +
        ggplot2::theme_bw()

      g2 <- flabel(dfs, trialsdt)

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
    } else { # add predicted subjects plot
      df1 <- enroll_pred_df %>%
        dplyr::mutate(parameter = "subjects")
      df1last <- df1 %>% dplyr::slice(dplyr::n())

      df2 <- dfs %>%
        dplyr::mutate(parameter = "events")
      df2last <- df2 %>% dplyr::slice(dplyr::n())

      # extend enrollment prediction time to event prediction time
      df12 <- df2last %>%
        dplyr::mutate(n = df1last$n,
                      lower = df1last$lower,
                      upper = df1last$upper,
                      parameter = df1last$parameter)

      dfs <- df1 %>%
        dplyr::bind_rows(df12) %>%
        dplyr::bind_rows(df2)

      # separate data into observed and predicted
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfs$date),
                                     max(dfs$date)) %/% months(1)
      bw = fbw(n_months)

      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$date,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper,
                                                    group=.data$parameter),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_step(data=dfa, ggplot2::aes(x=.data$date, y=.data$n,
                                                  group=.data$parameter),
                           color="black") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$date, y=.data$n,
                                                  group=.data$parameter),
                           color="blue") +
        ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_date(name = NULL,
                              labels = scales::date_format("%b"),
                              breaks = scales::breaks_width(bw),
                              minor_breaks = NULL,
                              expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Subjects / Events",
                      title = "Predicted subjects and events") +
        ggplot2::theme_bw()

      g2 <- flabel(dfs, trialsdt)

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
    }

    if (showplot) print(p1)


    if (showDropout) {
      # observed number of dropouts before data cut
      dfc <- df %>%
        dplyr::filter(.data$event == 1 | .data$dropout == 1) %>%
        dplyr::mutate(arrivalTime = as.numeric(.data$randdt-trialsdt+1)) %>%
        dplyr::mutate(totalTime = .data$arrivalTime + .data$time) %>%
        dplyr::arrange(.data$totalTime) %>%
        dplyr::mutate(t = .data$totalTime,
                      n = cumsum(.data$dropout),
                      lower = NA,
                      upper = NA) %>%
        dplyr::select(.data$t, .data$n, .data$lower, .data$upper)


      # predicted number of dropouts after data cut
      dfd = dplyr::tibble(t = t) %>%
        dplyr::cross_join(newEvents) %>%
        dplyr::group_by(.data$t, .data$draw) %>%
        dplyr::summarise(ndropouts = sum(.data$totalTime <= .data$t &
                                           .data$dropout == 1) + c0,
                         .groups = "drop_last") %>%
        dplyr::summarise(n = quantile(.data$ndropouts, probs = 0.5),
                         lower = quantile(.data$ndropouts, probs = plower),
                         upper = quantile(.data$ndropouts, probs = pupper))

      # add time zero and concatenate dropouts before and after data cut
      dft <- dplyr::tibble(t = 0, n = 0, lower = NA, upper = NA) %>%
        dplyr::bind_rows(dfc) %>%
        dplyr::bind_rows(dfd) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(year = format(.data$date, format = "%Y"))


      # separate data into observed and predicted
      dfc <- dft %>% dplyr::filter(is.na(.data$lower))
      dfd <- dft %>% dplyr::filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dft$date),
                                     max(dft$date)) %/% months(1)
      bw = fbw(n_months)

      # generate the plot
      g3 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfd, ggplot2::aes(x=.data$date,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_step(data=dfc, ggplot2::aes(x=.data$date, y=.data$n),
                           color="black") +
        ggplot2::geom_line(data=dfd, ggplot2::aes(x=.data$date, y=.data$n),
                           color="blue") +
        ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
        ggplot2::scale_x_date(name = NULL,
                              labels = scales::date_format("%b"),
                              breaks = scales::breaks_width(bw),
                              minor_breaks = NULL,
                              expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Dropouts", title = "Predicted dropouts") +
        ggplot2::theme_bw()

      g4 <- flabel(dft, trialsdt)

      # stack them together
      p2 <- g3 + g4 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      if (showplot) print(p2)

    }


    # whether to show ongoing subjects
    if (showOngoing) {
      dfx <- df %>%
        dplyr::mutate(arrivalTime = as.numeric(.data$randdt -
                                                 trialsdt + 1)) %>%
        dplyr::mutate(totalTime = .data$arrivalTime + .data$time) %>%
        dplyr::select(.data$arrivalTime, .data$totalTime, .data$event,
                      .data$dropout)

      tx = sort(unique(c(dfx$arrivalTime, dfx$totalTime)))

      # observed data
      dfe <- dplyr::tibble(t = tx) %>%
        dplyr::cross_join(dfx) %>%
        dplyr::group_by(.data$t) %>%
        dplyr::summarise(n = sum(.data$arrivalTime <= .data$t &
                                   .data$totalTime > .data$t),
                         .groups = "drop_last") %>%
        dplyr::mutate(lower = NA, upper = NA) %>%
        dplyr::select(.data$t, .data$n, .data$lower, .data$upper)

      # ongoing and new subjects
      dff = dplyr::tibble(t = t) %>%
        dplyr::cross_join(newEvents) %>%
        dplyr::group_by(.data$t, .data$draw) %>%
        dplyr::summarise(natrisk = sum(.data$arrivalTime <= .data$t &
                                         .data$totalTime > .data$t),
                         .groups = "drop_last") %>%
        dplyr::summarise(n = quantile(.data$natrisk, probs = 0.5),
                         lower = quantile(.data$natrisk, probs = plower),
                         upper = quantile(.data$natrisk, probs = pupper))

      # add time zero and concatenate data before and after data cut
      dfu <- dplyr::tibble(t = 0, n = 0, lower = NA, upper = NA) %>%
        dplyr::bind_rows(dfe) %>%
        dplyr::bind_rows(dff) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(year = format(.data$date, format = "%Y"))


      # separate data into observed and predicted
      dfe <- dfu %>% dplyr::filter(is.na(.data$lower))
      dff <- dfu %>% dplyr::filter(!is.na(.data$lower))

      n_months = lubridate::interval(min(dfu$date),
                                     max(dfu$date)) %/% months(1)
      bw = fbw(n_months)

      # generate the plot
      g5 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dff, ggplot2::aes(x=.data$date,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_step(data=dfe, ggplot2::aes(x=.data$date, y=.data$n),
                           color="black") +
        ggplot2::geom_line(data=dff, ggplot2::aes(x=.data$date, y=.data$n),
                           color="blue") +
        ggplot2::geom_vline(xintercept = cutoffdt, linetype = 2) +
        ggplot2::scale_x_date(name = NULL,
                              labels = scales::date_format("%b"),
                              breaks = scales::breaks_width(bw),
                              minor_breaks = NULL,
                              expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Ongoing subjects",
                      title = "Predicted ongoing subjects") +
        ggplot2::theme_bw()

      g6 <- flabel(dfu, trialsdt)

      # stack them together
      p3 <- g5 + g6 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      if (showplot) print(p3)

    }


    if (showDropout) {
      if (showOngoing) {
        list(event_pred_day = pred2dy, event_pred_date = pred2dt,
             pilevel = pilevel, newEvents = newEvents,
             event_pred_df = dfs, dropout_pred_df = dft,
             ongoing_pred_df = dfu,
             event_pred_plot = p1, dropout_pred_plot = p2,
             ongoing_pred_plot = p3)
      } else {
        list(event_pred_day = pred2dy, event_pred_date = pred2dt,
             pilevel = pilevel, newEvents = newEvents,
             event_pred_df = dfs, dropout_pred_df = dft,
             event_pred_plot = p1, dropout_pred_plot = p2)
      }
    } else {
      if (showOngoing) {
        list(event_pred_day = pred2dy, event_pred_date = pred2dt,
             pilevel = pilevel, newEvents = newEvents,
             event_pred_df = dfs, ongoing_pred_df = dfu,
             event_pred_plot = p1, ongoing_pred_plot = p3)
      } else {
        list(event_pred_day = pred2dy, event_pred_date = pred2dt,
             pilevel = pilevel, newEvents = newEvents,
             event_pred_df = dfs, event_pred_plot = p1)
      }
    }

  } else {
    if (is.null(newSubjects)) {
      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$t,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$t, y=.data$n),
                           color="blue") +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_continuous(name = "Days since randomization",
                                    expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Events", title = "Predicted events") +
        ggplot2::theme_bw()
    } else {
      df1 <- enroll_pred_df %>%
        dplyr::mutate(parameter = "subjects")
      df1last <- df1 %>% dplyr::slice(dplyr::n())

      df2 <- dfb %>%
        dplyr::mutate(parameter = "events")
      df2last <- df2 %>% dplyr::slice(dplyr::n())

      # extend enrollment prediction time to event prediction time
      df12 <- df2last %>%
        dplyr::mutate(n = df1last$n,
                      lower = df1last$lower,
                      upper = df1last$upper,
                      parameter = df1last$parameter)

      dfb <- df1 %>%
        dplyr::bind_rows(df12) %>%
        dplyr::bind_rows(df2)

      g1 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfb, ggplot2::aes(x=.data$t,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper,
                                                    group=.data$parameter),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_line(data=dfb, ggplot2::aes(x=.data$t, y=.data$n,
                                                  group=.data$parameter),
                           color="blue") +
        ggplot2::geom_hline(yintercept = target_d, linetype = 2) +
        ggplot2::scale_x_continuous(name = "Days since randomization",
                                    expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Subjects / Events",
                      title = "Predicted subjects and events") +
        ggplot2::theme_bw()
    }

    if (showplot) print(g1)


    if (showDropout) {
      dfd = dplyr::tibble(t = t) %>%
        dplyr::cross_join(newEvents) %>%
        dplyr::group_by(.data$t, .data$draw) %>%
        dplyr::summarise(ndropouts = sum(.data$totalTime <= .data$t &
                                           .data$dropout == 1) + c0,
                         .groups = "drop_last") %>%
        dplyr::summarise(n = quantile(.data$ndropouts, probs = 0.5),
                         lower = quantile(.data$ndropouts, probs = plower),
                         upper = quantile(.data$ndropouts, probs = pupper))


      g3 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dfd, ggplot2::aes(x=.data$t,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_line(data=dfd, ggplot2::aes(x=.data$t, y=.data$n),
                           color="blue") +
        ggplot2::scale_x_continuous(name = "Days since randomization",
                                    expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Dropouts", title = "Predicted dropouts") +
        ggplot2::theme_bw()

      if (showplot) print(g3)

    }


    if (showOngoing) {
      dff = dplyr::tibble(t = t) %>%
        dplyr::cross_join(newEvents) %>%
        dplyr::group_by(.data$t, .data$draw) %>%
        dplyr::summarise(natrisk = sum(.data$arrivalTime <= .data$t &
                                         .data$totalTime > .data$t),
                         .groups = "drop_last") %>%
        dplyr::summarise(n = quantile(.data$natrisk, probs = 0.5),
                         lower = quantile(.data$natrisk, probs = plower),
                         upper = quantile(.data$natrisk, probs = pupper))


      g5 <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=dff, ggplot2::aes(x=.data$t,
                                                    ymin=.data$lower,
                                                    ymax=.data$upper),
                             alpha=0.5, fill="lightblue") +
        ggplot2::geom_line(data=dff, ggplot2::aes(x=.data$t, y=.data$n),
                           color="blue") +
        ggplot2::scale_x_continuous(name = "Days since randomization",
                                    expand = c(0.01, 0.01)) +
        ggplot2::labs(y = "Ongoing subjects",
                      title = "Predicted ongoing subjects") +
        ggplot2::theme_bw()

      if (showplot) print(g5)

    }


    if (showDropout) {
      if (showOngoing) {
        list(event_pred_day = pred2dy, pilevel = pilevel,
             newEvents = newEvents,
             event_pred_df = dfb, dropout_pred_df = dfd,
             ongoing_pred_df = dff,
             event_pred_plot = g1, dropout_pred_plot = g3,
             ongoing_pred_df = g5)
      } else {
        list(event_pred_day = pred2dy, pilevel = pilevel,
             newEvents = newEvents,
             event_pred_df = dfb, dropout_pred_df = dfd,
             event_pred_plot = g1, dropout_pred_plot = g3)
      }
    } else {
      if (showOngoing) {
        list(event_pred_day = pred2dy, pilevel = pilevel,
             newEvents = newEvents,
             event_pred_df = dfb, ongoing_pred_df = dff,
             event_pred_plot = g1, ongoing_pred_plot = g5)
      } else {
        list(event_pred_day = pred2dy, pilevel = pilevel,
             newEvents = newEvents,
             event_pred_df = dfb,
             event_pred_plot = g1)
      }
    }
  }
}
