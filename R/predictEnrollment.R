#' @title Predict enrollment
#' @description Utilizes a pre-fitted enrollment model to generate
#'   enrollment times for new subjects and provide a prediction
#'   interval for the expected time to reach the enrollment target.
#'
#' @param df The subject-level enrollment data, including \code{trialsdt},
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
#' @param nyears The number of years after the data cut for prediction.
#'   By default, it is set to 4.
#' @param nreps The number of replications for simulation. By default,
#'   it is set to 500.
#' @param showsummary A Boolean variable to control whether or not to
#'   show the prediction summary. By default, it is set to \code{TRUE}.
#' @param showplot A Boolean variable to control whether or not to
#'   show the prediction plot. By default, it is set to \code{TRUE}.
#'
#' @details
#'
#' The \code{enroll_fit} variable can be used for enrollment prediction
#' at the design stage. A piecewise Poisson can be parameterized
#' through the time intervals, \code{accrualTime}, which is
#' treated as fixed, and the enrollment rates in the intervals,
#' \code{accrualIntensity}, the log of which is used as the
#' model parameter. For the homogeneous Poisson, time-decay,
#' and piecewise Poisson models, \code{enroll_fit} is used to
#' specify the prior distribution of model parameters, with
#' a very small variance being used to fix the parameter values.
#' It should be noted that the B-spline model is not appropriate
#' for use during the design stage.
#'
#' @return
#' A list of prediction results, which includes important information
#' such as the median, lower and upper percentiles for the estimated
#' time to reach the target number of subjects, as well as simulated
#' enrollment data for new subjects. The data for the
#' prediction plot is also included within the list.
#'
#' @examples
#'
#' # Enrollment prediction at the design stage
#'
#' enroll_pred <- predictEnrollment(
#'   target_n = 300,
#'   enroll_fit = list(model = "piecewise poisson",
#'                     theta = log(26/9*seq(1, 9)/30.4375),
#'                     vtheta = diag(9)*1e-8,
#'                     accrualTime = seq(0, 8)*30.4375),
#'   pilevel = 0.90, nreps = 100)
#'
#'
#' @export
#'
predictEnrollment <- function(df = NULL, target_n, enroll_fit, lags = 30,
                              pilevel = 0.90, nyears = 4, nreps = 500,
                              showsummary = TRUE, showplot = TRUE) {
  if (!is.null(df)) erify::check_class(df, "data.frame")
  erify::check_n(target_n)

  erify::check_class(enroll_fit, "list")
  erify::check_content(tolower(enroll_fit$model), c(
    "poisson", "time-decay", "b-spline", "piecewise poisson"))

  erify::check_n(lags, zero = TRUE)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(nreps)
  erify::check_bool(showsummary)
  erify::check_bool(showplot)

  if (!is.null(df)) {
    df <- dplyr::as_tibble(df)
    names(df) <- tolower(names(df))
    df$trialsdt <- as.Date(df$trialsdt)
    df$randdt <- as.Date(df$randdt)
    df$cutoffdt <- as.Date(df$cutoffdt)

    trialsdt = df$trialsdt[1]
    cutoffdt = df$cutoffdt[1]
    n0 = nrow(df)
    t0 = as.numeric(cutoffdt - trialsdt + 1)
    df <- df %>%
      dplyr::arrange(.data$randdt) %>%
      dplyr::mutate(t = as.numeric(.data$randdt - trialsdt + 1),
                    n = dplyr::row_number())
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
      df = dplyr::as_tibble(matrix(
        nrow = nreps*n1, ncol = 2,
        dimnames = list(NULL, c("draw", "arrivalTime"))))
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

    # equation to solve for t
    fenroll <- function(t, theta, muTime) {
      fmu_td(t, theta) - muTime
    }


    # draw arrival time for new subjects
    newEnrollment_td <- function(t0, n1, theta, nreps) {
      df = dplyr::as_tibble(matrix(
        nrow = nreps*n1, ncol = 2,
        dimnames = list(NULL, c("draw", "arrivalTime"))))
      for (i in 1:nreps) {
        index = (i-1)*n1 + (1:n1)
        df[index, "draw"] = i
        gapmuTime = rexp(n1)
        muTime = cumsum(gapmuTime) + fmu_td(t0, theta[i,])

        mu = exp(theta[i,1])
        delta = exp(theta[i,2])
        # find the tangent line with half of maximum slope:
        #   v(t) = mu(ti) + mu/(2*delta)*(t-ti)
        # which lies beneath mu(t), and then find tmax such that
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
      t0x = nrow(lambda)  # to account for enrollment pause
      lambdaT = colMeans(lambda[(t0x - lags):t0x,])

      df = dplyr::as_tibble(matrix(
        nrow = nreps*n1, ncol = 2,
        dimnames = list(NULL, c("draw", "arrivalTime"))))
      for (i in 1:nreps) {
        index = (i-1)*n1 + (1:n1)
        df[index, "draw"] = i
        gapTime = rexp(n1, lambdaT[i])
        df[index, "arrivalTime"] = cumsum(gapTime) + t0
      }

      df
    }

    newSubjects <- newEnrollment_bs(t0, n1, theta, x = enroll_fit$x,
                                    lags, nreps)
  } else if (tolower(enroll_fit$model) == "piecewise poisson") {
    # draw parameter from posterior distribution
    if (length(enroll_fit$theta) == 1) {
      theta = matrix(rnorm(nreps, mean = enroll_fit$theta,
                           sd = sqrt(enroll_fit$vtheta)), ncol=1)
    } else {
      theta = mvtnorm::rmvnorm(nreps, mean = enroll_fit$theta,
                               sigma = enroll_fit$vtheta)
    }
    u = enroll_fit$accrualTime

    # mu(t[j]) - mu(t[j-1]) is standard exponential distribution, j=1,...,n1
    newEnrollment_pw <- function(t0, n1, theta, u, nreps) {
      df = dplyr::as_tibble(matrix(
        nrow = nreps*n1, ncol = 2,
        dimnames = list(NULL, c("draw", "arrivalTime"))))
      J = length(u)
      j0 = findInterval(t0, u)
      for (i in 1:nreps) {
        index = (i-1)*n1 + (1:n1)
        df[index, "draw"] = i
        a = exp(theta[i,]) # enrollment rate in each interval
        if (J>1) {
          psum = c(0, cumsum(a[1:(J-1)] * diff(u)))
        } else {
          psum = 0
        }
        rhs =  psum[j0] + a[j0]*(t0 - u[j0]) + cumsum(rexp(n1))
        j1 = findInterval(rhs, psum)
        df[index, "arrivalTime"] = u[j1] + (rhs - psum[j1])/a[j1]
      }

      df
    }

    newSubjects <- newEnrollment_pw(t0, n1, theta, u, nreps)
  }


  # lower and upper percentages
  plower = (1 - pilevel)/2
  pupper = 1 - plower

  # find the arrivalTime of the last subject for each simulated data set
  new1 <- newSubjects %>%
    dplyr::group_by(.data$draw) %>%
    dplyr::slice(dplyr::n())

  pred_day <- ceiling(quantile(new1$arrivalTime, c(0.5, plower, pupper)))

  t1 = t0 + nyears*365 # extend time to nyears after cutoff

  # future time points at which to predict number of subjects
  t = sort(unique(c(seq(t0, t1, 30), t1, pred_day)))
  t = t[t <= t1] # restrict range of x-axis

  # predicted number of subjects enrolled after data cut
  dfb <- dplyr::tibble(t = t) %>%
    dplyr::cross_join(newSubjects) %>%
    dplyr::group_by(.data$t, .data$draw) %>%
    dplyr::summarise(nenrolled = sum(.data$arrivalTime <= .data$t) + n0,
                     .groups = "drop_last") %>%
    dplyr::summarise(n = quantile(.data$nenrolled, probs = 0.5),
                     lower = quantile(.data$nenrolled, probs = plower),
                     upper = quantile(.data$nenrolled, probs = pupper),
                     mean = mean(.data$nenrolled),
                     var = var(.data$nenrolled))


  if (!is.null(df)) {
    pred_date <- as.Date(pred_day - 1, origin = trialsdt)

    str1 <- paste0("Time from cutoff until ", target_n, " subjects: ",
                   pred_date[1] - cutoffdt + 1, " days")
    str2 <- paste0("Median prediction date: ", pred_date[1])
    str3 <- paste0("Prediction interval: ", pred_date[2], ", ", pred_date[3])
    s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")

    # arrival time for subjects already enrolled before data cut
    dfa <- df %>%
      dplyr::mutate(lower = NA, upper = NA, mean = .data$n, var = 0) %>%
      dplyr::select(.data$t, .data$n, .data$lower, .data$upper,
                    .data$mean, .data$var) %>%
      dplyr::group_by(.data$t) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::ungroup()

    # extend observed to cutoff date
    dfa1 <- dfa %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::mutate(t = t0)

    if (max(dfa$t) < t0) {
      dfa <- dfa %>%
        dplyr::bind_rows(dfa1)
    }

    # concatenate subjects enrolled before and after data cut
    dfs <- dfa %>%
      dplyr::bind_rows(dfb) %>%
      dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt))


    # separate data into observed and predicted
    dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
    dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

    # plot the enrollment data with month as x-axis label
    g1 <- plotly::plot_ly() %>%
      plotly::add_ribbons(data = dfb, x = ~date,
                          ymin = ~lower, ymax = ~upper,
                          name = "prediction interval",
                          fill = "tonexty",
                          line = list(width=0)) %>%
      plotly::add_lines(data = dfb, x = ~date, y = ~n,
                        name = "median prediction",
                        line = list(width=1)) %>%
      plotly::add_lines(data = dfa, x = ~date, y = ~n,
                        name = "observed",
                        line = list(shape="hv", width=1)) %>%
      plotly::add_lines(x = rep(cutoffdt, 2), y = range(dfs$n),
                        name = "cutoff",
                        line = list(dash="dash"),
                        showlegend = FALSE) %>%
      plotly::layout(
        annotations = list(x = cutoffdt, y = 0,
                           text = 'cutoff', xanchor = "left",
                           yanchor = "bottom",
                           font = list(size=12),
                           showarrow = FALSE),
        xaxis = list(title = "", zeroline = FALSE),
        yaxis = list(title = "Subjects", zeroline = FALSE),
        legend = list(x = 0, y = 1.05, yanchor = "bottom",
                      orientation = "h"))
  } else {
    str1 <- paste0("Time from trial start until ", target_n, " subjects")
    str2 <- paste0("Median prediction day: ", pred_day[1])
    str3 <- paste0("Prediction interval: ", pred_day[2], ", ", pred_day[3])
    s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")

    g1 <- plotly::plot_ly(dfb, x = ~t) %>%
      plotly::add_ribbons(ymin = ~lower, ymax = ~upper,
                          name = "prediction interval",
                          fill = "tonexty",
                          line = list(width=0)) %>%
      plotly::add_lines(y = ~n, name = "median prediction",
                        line = list(width=1)) %>%
      plotly::layout(xaxis = list(title = "Days since trial start",
                                  zeroline = FALSE),
                     yaxis = list(title = "Subjects", zeroline = FALSE),
                     legend = list(x = 0, y = 1.05, yanchor = "bottom",
                                   orientation = "h"))
  }


  if (showsummary) cat(s1)
  if (showplot) print(g1)

  if (!is.null(df)) {
    list(target_n = target_n, enroll_pred_day = pred_day,
         enroll_pred_date = pred_date,
         pilevel = pilevel, newSubjects = newSubjects,
         enroll_pred_df = dfs,
         enroll_pred_summary = s1, enroll_pred_plot = g1)
  } else {
    list(target_n = target_n, enroll_pred_day = pred_day,
         pilevel = pilevel, newSubjects = newSubjects,
         enroll_pred_df = dfb,
         enroll_pred_summary = s1, enroll_pred_plot = g1)
  }

}
