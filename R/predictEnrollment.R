#' @title Predict enrollment
#' @description Use a fitted enrollment model to generate enrollment time
#'   for new subjects and provide prediction interval for time to reach
#'   the enrollment target.
#'
#' @param target_n Target number of subjects to enroll.
#' @param df Observed subject-level enrollment data.
#' @param fit Enrollment model fit to observed data.
#' @param lags Day lags for averaging enrollment rates from B-spline.
#'   Defaults to 30.
#' @param pilevel Prediction interval level. Defaults to 0.90.
#' @param nreps Number of replications for simulation. Defaults to 500.
#'
#' @return A list of prediction results consisting of the median time to
#'   reach enrollment target, date of median, lower and upper percentiles
#'   to reach enrollment target, and simulated enrollment data for new
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
#'                             fit = fitEnr, lags = 30, pilevel = 0.90,
#'                             nreps = 500)
#'
#' @export
#'
predictEnrollment <- function(target_n, df, fit, lags = 30, pilevel = 0.90,
                              nreps = 500) {
  erify::check_n(target_n)
  erify::check_class(df, "data.frame")
  erify::check_class(fit, "list")
  erify::check_n(lags, zero = TRUE)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(nreps)

  trialsdt = min(df$RANDDT)
  cutoffdt = df$CUTOFFDT[1]
  n0 = nrow(df)
  t0 = as.numeric(cutoffdt - trialsdt + 1)
  n1 = target_n - n0  # number of new subjects

  erify::check_n(n1)


  if (tolower(fit$model) == "poisson") {
    # draw parameter from posterior distribution
    theta = rnorm(nreps, mean = fit$theta, sd = sqrt(fit$vtheta))

    # draw arrival time for new subjects
    newEnrollment_po <- function(t0, n1, theta, nreps) {
      lambda = exp(theta)
      df = data.frame(matrix(nrow=nreps*n1, ncol=2))
      colnames(df) = c("draw", "arrivalTime")
      if (n1 > 0) {
        for (i in 1:nreps) {
          index = (i-1)*n1 + (1:n1)
          df[index, "draw"] = i
          gapTime = rexp(n1, lambda[i])
          df[index, "arrivalTime"] = cumsum(gapTime) + t0
        }
      }

      df
    }

    newSubjects <- newEnrollment_po(t0, n1, theta, nreps)
  } else if (tolower(fit$model) == "time-decay") {
    # draw parameter from posterior distribution
    theta = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)

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
      if (n1 > 0) {
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
      }

      df
    }

    newSubjects <- newEnrollment_td(t0, n1, theta, nreps)
  } else if (tolower(fit$model) == "b-spline") {
    # draw parameter from posterior distribution
    theta = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)

    newEnrollment_bs <- function(t0, n1, theta, x, lags, nreps) {
      lambda = exp(x %*% t(theta))
      # moving average for enrollment rate after t0
      lambdaT = colMeans(lambda[(t0 - lags):t0,])

      df = data.frame(matrix(nrow=nreps*n1, ncol=2))
      colnames(df) = c("draw", "arrivalTime")
      if (n1 > 0) {
        for (i in 1:nreps) {
          index = (i-1)*n1 + (1:n1)
          df[index, "draw"] = i
          gapTime = rexp(n1, lambdaT[i])
          df[index, "arrivalTime"] = cumsum(gapTime) + t0
        }
      }

      df
    }

    newSubjects <- newEnrollment_bs(t0, n1, theta, fit$x, lags, nreps)
  }



  # lower and upper percentages
  plower = (1 - pilevel)/2
  pupper = 1 - plower

  new1 <- newSubjects %>%
    group_by(.data$draw) %>%
    filter(row_number() == n())

  pred1dy <- ceiling(quantile(new1$arrivalTime, c(0.5, plower, pupper)))
  pred1dt <- as.Date(pred1dy - 1, origin = trialsdt)
  r1month <- round((pred1dy[1] - t0 + 1)/30.4375, digits = 1)

  t1 = pred1dy[3] + 30 # extend time to 30 days after

  # future time points at which to predict number of subjects
  t = c(seq(t0, pred1dy[2], by=30), seq(pred1dy[2], pred1dy[3]), t1)

  # arrival time for subjects already enrolled before data cut
  dfa <- df %>%
    mutate(lower = NA, upper = NA) %>%
    dplyr::select(.data$time, .data$n, .data$lower, .data$upper)

  # predicted number of subjects enrolled after data cut
  if (n1 > 0) {
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
  } else {
    dfb <- tibble(time = t1, n = n0, lower = NA, upper = NA)
  }


  # concatenate subjects enrolled before and after data cut
  df1 <- dfa %>%
    bind_rows(dfb) %>%
    mutate(date = as.Date(.data$time - 1, origin = trialsdt)) %>%
    mutate(year = format(.data$date, format = "%Y"))

  # separate data into observed and predicted
  df1a <- df1 %>% filter(is.na(.data$lower))
  df1b <- df1 %>% filter(!is.na(.data$lower))



  n_months = lubridate::interval(min(df1$date), max(df1$date)) %/% months(1)
  bw = fbw(n_months)

  g2 <- flabel(df1, trialsdt)

  # plot the enrollment data with month as x-axis label
  g1 <- ggplot() +
    geom_ribbon(data=df1b,
                aes(x=.data$date, ymin=.data$lower, ymax=.data$upper),
                alpha=0.5, fill="lightblue") +
    geom_step(data=df1a, aes(x=.data$date, y=.data$n), color="black") +
    geom_line(data=df1b, aes(x=.data$date, y=.data$n), color="blue") +
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

  list(remEnrMonth = r1month, predEnrDT = pred1dt, pilevel = pilevel,
       newSubjects = newSubjects)
}
