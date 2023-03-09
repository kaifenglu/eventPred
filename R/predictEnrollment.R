#' @title Predict enrollment
#' @description Utilizes a pre-fitted enrollment model to generate
#'   enrollment times for new subjects and also provide a prediction
#'   interval for the expected time to reach the enrollment target.
#'
#' @param df The subject-level enrollment data, including
#'   \code{randdt} and \code{cutoffdt}. By default, it is set to
#'   \code{NULL} for enrollment prediction at the design stage.
#' @param target_n The target number of subjects to enroll in the study.
#' @param enroll_fit the pre-fitted enrollment model used to
#'   generate predictions.
#' @param lags The day lags to compute the average enrollment rate to
#'   carry forward for the B-spline enrollment model. By default,
#'   it is set to 30.
#' @param pilevel The prediction interval level. By default,
#'   it is set to 0.90.
#' @param nreps The number of replications for simulation. By default,
#'   it is set to 500.
#' @param showplot A Boolean variable to control whether or not
#'   the prediction plot is displayed. By default, it is set to
#'   \code{TRUE}.
#'
#' @details
#' The \code{enroll_fit} variable can be used for enrollment
#' prediction at the design stage. A piecewise Poisson enrollment
#' model can be parameterized through the time intervals,
#' \code{accrualTime}, and the enrollment rates in the intervals,
#' \code{accrualIntensity}. These are treated as fixed for
#' design-stage enrollment prediction.
#' For the homogeneous Poisson and time-decay models,
#' \code{enroll_fit} is used to specify the prior distribution of
#' model parameters, with a very small variance being used to fix
#' the parameter values. It should be noted that the B-spline model
#' is not appropriate for use during the design stage.
#'
#' @return
#' A list of prediction results, which includes important information
#' such as the median, lower and upper percentiles for the estimated
#' time to reach the target number of subjects, as well as simulated
#' enrollment data for new subjects. Additionally, the data for the
#' prediction plot is also included within the list.
#'
#' @examples
#'
#' # Example 1: Enrollment prediction at analysis stage
#'
#' enroll_fit <- fitEnrollment(
#'   df = observedData, enroll_model = "b-spline", nknots = 1)
#'
#' enroll_pred <- predictEnrollment(
#'   df = observedData, target_n = 480, enroll_fit = enroll_fit,
#'   lags = 30, pilevel = 0.90, nreps = 500)
#'
#' # Example 2: Enrollment prediction at design stage
#'
#' enroll_pred <- predictEnrollment(
#'   target_n = 480,
#'   enroll_fit = list(model = "piecewise poisson",
#'                     accrualTime = seq(0, 8)*30.4375,
#'                     accrualIntensity = 26/9*seq(1, 9)/30.4375),
#'   pilevel = 0.90, nreps = 500)
#'
#'
#' @export
#'
predictEnrollment <- function(df = NULL, target_n, enroll_fit, lags = 30,
                              pilevel = 0.90, nreps = 500,
                              showplot = TRUE) {
  if (!is.null(df)) erify::check_class(df, "data.frame")
  erify::check_n(target_n)
  erify::check_class(enroll_fit, "list")
  erify::check_n(lags, zero = TRUE)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(nreps)
  erify::check_bool(showplot)

  if (!is.null(df)) {
    names(df) <- tolower(names(df))
    trialsdt = min(df$randdt)
    cutoffdt = df$cutoffdt[1]
    n0 = nrow(df)
    t0 = as.numeric(cutoffdt - trialsdt + 1)
    df <- df %>%
      arrange(.data$randdt) %>%
      mutate(time = as.numeric(.data$randdt - trialsdt + 1),
             n = row_number())
  } else {
    n0 = 0
    t0 = 1
  }
  n1 = target_n - n0  # number of new subjects

  erify::check_n(n1)


  if (tolower(enroll_fit$model) == "poisson") {
    # draw parameter from posterior distribution
    theta = rnorm(nreps, mean = enroll_fit$theta,
                  sd = sqrt(enroll_fit$vtheta))

    # draw arrival time for new subjects
    newEnrollment_po <- function(t0, n1, theta, nreps) {
      lambda = exp(theta)
      df = data.frame(matrix(nrow=nreps*n1, ncol=2))
      colnames(df) = c("draw", "arrivalTime")
      for (i in 1:nreps) {
        index = (i-1)*n1 + (1:n1)
        df[index, "draw"] = i
        gapTime = rexp(n1, lambda[i])
        df[index, "arrivalTime"] = cumsum(gapTime) + t0
      }

      df
    }

    newSubjects <- newEnrollment_po(t0, n1, theta, nreps)
  } else if (tolower(enroll_fit$model) == "time-decay") {
    # draw parameter from posterior distribution
    theta = mvtnorm::rmvnorm(nreps, mean = enroll_fit$theta,
                             sigma = enroll_fit$vtheta)

    # mean function of the NHPP
    fmu_td <- function(t, theta) {
      mu = exp(theta[1])
      delta = exp(theta[2])
      mu/delta*(t - 1/delta*(1 - exp(-delta*t)))
    }


    # draw arrival time for new subjects
    newEnrollment_td <- function(t0, n1, theta, nreps) {
      df = data.frame(matrix(nrow=nreps*n1, ncol=2))
      colnames(df) = c("draw", "arrivalTime")
      for (i in 1:nreps) {
        index = (i-1)*n1 + (1:n1)
        df[index, "draw"] = i
        gapmuTime = rexp(n1)
        muTime = cumsum(gapmuTime) + fmu_td(t0, theta[i,])

        # equation to solve for t
        fenroll <- function(t, theta, muTime) {
          fmu_td(t, theta) - muTime
        }

        mu = exp(theta[i,1])
        delta = exp(theta[i,2])
        # find the tangent line with half of maximum slope:
        #   v(t) = mu(ti) + mu/(2*delta)*(t-ti)
        # which lies between mu(t), and then find tmax such that
        #   v(tmax) = muTime, which implies mu(tmax) > muTime
        ti = log(2)/delta
        tmax = (muTime - fmu_td(ti, theta[i,]))*2*delta/mu + ti
        interval = cbind(t0, tmax)

        # draw arrival time
        df[index, "arrivalTime"] = rstpm2::vuniroot(
          fenroll, interval, theta = theta[i,], muTime)$root
      }

      df
    }

    newSubjects <- newEnrollment_td(t0, n1, theta, nreps)
  } else if (tolower(enroll_fit$model) == "b-spline") {
    if (is.null(df)) {
      stop("B-spline enrollment model cannot be used at the design stage.")
    }

    # draw parameter from posterior distribution
    theta = mvtnorm::rmvnorm(nreps, mean = enroll_fit$theta,
                             sigma = enroll_fit$vtheta)

    newEnrollment_bs <- function(t0, n1, theta, x, lags, nreps) {
      lambda = exp(x %*% t(theta))
      # moving average for enrollment rate after t0
      lambdaT = colMeans(lambda[(t0 - lags):t0,])

      df = data.frame(matrix(nrow=nreps*n1, ncol=2))
      colnames(df) = c("draw", "arrivalTime")
      for (i in 1:nreps) {
        index = (i-1)*n1 + (1:n1)
        df[index, "draw"] = i
        gapTime = rexp(n1, lambdaT[i])
        df[index, "arrivalTime"] = cumsum(gapTime) + t0
      }

      df
    }

    newSubjects <- newEnrollment_bs(t0, n1, theta, enroll_fit$x, lags, nreps)
  } else if (tolower(enroll_fit$model) == "piecewise poisson") {
    if (!is.null(df)) {
      stop("Piecewise Poisson model is only used at the design stage.")
    }

    u = enroll_fit$accrualTime
    a = enroll_fit$accrualIntensity

    newEnrollment_pw <- function(t0, n1, u, a, nreps) {
      df = data.frame(matrix(nrow=nreps*n1, ncol=2))
      colnames(df) = c("draw", "arrivalTime")
      J = length(a)
      for (i in 1:nreps) {
        index = (i-1)*n1 + (1:n1)
        df[index, "draw"] = i
        rhs = cumsum(rexp(n1))
        psum = c(0, cumsum(a[1:(J-1)] * diff(u)))
        j1 = findInterval(rhs, psum)
        df[index, "arrivalTime"] = u[j1] + (rhs - psum[j1])/a[j1] + t0
      }

      df
    }

    newSubjects <- newEnrollment_pw(t0, n1, u, a, nreps)
  }


  # lower and upper percentages
  plower = (1 - pilevel)/2
  pupper = 1 - plower

  new1 <- newSubjects %>%
    group_by(.data$draw) %>%
    filter(row_number() == n())

  pred1dy <- ceiling(quantile(new1$arrivalTime, c(0.5, plower, pupper)))

  t1 = pred1dy[3] + 30 # extend time to 30 days after

  # future time points at which to predict number of subjects
  t = c(seq(t0, pred1dy[2], by=30), seq(pred1dy[2], pred1dy[3]), t1)

  # predicted number of subjects enrolled after data cut
  dfb = data.frame(matrix(nrow=length(t), ncol=4))
  colnames(dfb) = c('time', 'n', 'lower', 'upper')
  for (i in 1:length(t)) {
    # number of subjects enrolled after data cut in each simulated data set
    # adding to the number of subjects already enrolled before data cut
    sumdata <- newSubjects %>%
      group_by(.data$draw) %>%
      summarize(n = sum(.data$arrivalTime <= t[i]) + n0)

    # summary across simulated data sets
    dfb[i, 'time'] = t[i]
    dfb[i, c('n', 'lower', 'upper')] =
      quantile(sumdata$n, probs = c(0.5, plower, pupper))
  }


  if (!is.null(df)) {
    pred1dt <- as.Date(pred1dy - 1, origin = trialsdt)

    # arrival time for subjects already enrolled before data cut
    dfa <- df %>%
      mutate(lower = NA, upper = NA) %>%
      dplyr::select(.data$time, .data$n, .data$lower, .data$upper)

    # concatenate subjects enrolled before and after data cut
    dfs <- dfa %>%
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

      # plot the enrollment data with month as x-axis label
      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$date, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_step(data=dfa, aes(x=.data$date, y=.data$n), color="black") +
        geom_line(data=dfb, aes(x=.data$date, y=.data$n), color="blue") +
        geom_vline(xintercept = cutoffdt, linetype = 2) +
        scale_x_date(name = NULL,
                     labels = scales::date_format("%b"),
                     breaks = scales::breaks_width(bw),
                     minor_breaks = NULL,
                     expand = c(0.01, 0.01)) +
        labs(y = "Subjects", title = "Predicted subject enrollment") +
        theme_bw()

      # stack them together
      p1 <- g1 + g2 + patchwork::plot_layout(nrow = 2, heights = c(15, 1))
      print(p1)
    }

    list(predEnrollDay = pred1dy, predEnrollDate = pred1dt, pilevel = pilevel,
         newSubjects = newSubjects, plotdata = dfs)
  } else {
    if (showplot) {
      # plot the enrollment data
      g1 <- ggplot() +
        geom_ribbon(data=dfb,
                    aes(x=.data$time, ymin=.data$lower, ymax=.data$upper),
                    alpha=0.5, fill="lightblue") +
        geom_line(data=dfb, aes(x=.data$time, y=.data$n), color="blue") +
        scale_x_continuous(name = "Days since randomization",
                           expand = c(0.01, 0.01)) +
        labs(y = "Subjects", title = "Predicted subject enrollment") +
        theme_bw()

      print(g1)
    }

    list(predEnrollDay = pred1dy, pilevel = pilevel, newSubjects = newSubjects,
         plotdata = dfb)
  }


}
