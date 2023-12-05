#' @title Predict enrollment
#' @description Utilizes a pre-fitted enrollment model to generate
#'   enrollment times for new subjects and provide a prediction
#'   interval for the expected time to reach the enrollment target.
#'
#' @param df The subject-level enrollment data, including \code{trialsdt},
#'   \code{randdt} and \code{cutoffdt}. The data should also include
#'   \code{treatment} coded as 1, 2, and so on, and
#'   \code{treatment_description} for prediction
#'   by treatment group. By default, it is set to \code{NULL}
#'   for enrollment prediction at the design stage.
#' @param target_n The target number of subjects to enroll in the study.
#' @param enroll_fit The pre-fitted enrollment model used to
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
#' @param by_treatment A Boolean variable to control whether or not to
#'   predict enrollment by treatment group. By default,
#'   it is set to \code{FALSE}.
#' @param ngroups The number of treatment groups for enrollment prediction
#'   at the design stage. By default, it is set to 1.
#'   It is replaced with the actual number of
#'   treatment groups in the observed data if \code{df} is not \code{NULL}.
#' @param alloc The treatment allocation in a randomization block.
#'   By default, it is set to \code{NULL}, which yields equal allocation
#'   among the treatment groups.
#' @param treatment_label The treatment labels for treatments in a
#'   randomization block for design stage prediction.
#'   It is replaced with the treatment_description
#'   in the observed data if \code{df} is not \code{NULL}.
#'
#' @details
#' The \code{enroll_fit} variable can be used for enrollment prediction
#' at the design stage. A piecewise Poisson model can be parameterized
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
#' During the enrollment stage, \code{enroll_fit} is the enrollment model
#' fit based on the observed data. The fitted enrollment model is used to
#' generate enrollment times for new subjects.
#'
#' @return
#' A list of prediction results, which includes important information
#' such as the median, lower and upper percentiles for the estimated
#' time to reach the target number of subjects, as well as simulated
#' enrollment data for new subjects. The data for the
#' prediction plot is also included within the list.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Xiaoxi Zhang and Qi Long. Stochastic modeling and prediction for
#' accrual in clinical trials. Stat in Med. 2010; 29:649-658.
#'
#' @examples
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
#' @export
#'
predictEnrollment <- function(df = NULL, target_n, enroll_fit, lags = 30,
                              pilevel = 0.90, nyears = 4, nreps = 500,
                              showsummary = TRUE, showplot = TRUE,
                              by_treatment = FALSE, ngroups = 1,
                              alloc = NULL, treatment_label = NULL) {
  if (!is.null(df)) erify::check_class(df, "data.frame")
  erify::check_n(target_n)

  erify::check_class(enroll_fit, "list")
  erify::check_content(tolower(enroll_fit$model), c(
    "poisson", "time-decay", "b-spline", "piecewise poisson"))

  model = tolower(enroll_fit$model)
  p = length(enroll_fit$theta)
  vtheta = enroll_fit$vtheta

  if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                 ncol(vtheta) != p)) ||
      (p == 1 && length(c(vtheta)) != 1)) {
    stop(paste("Dimensions of vtheta must be compatible with the length",
               "of theta in enroll_fit"))
  }

  if ((model == "poisson" && p != 1) ||
      (model == "time-decay" && p != 2) ||
      (model == "piecewise poisson" &&
       p != length(enroll_fit$accrualTime)) ||
      (model == "b-spline" && p != ncol(enroll_fit$x))) {
    stop(paste("Length of theta must be compatible with model",
               "in enroll_fit"))
  }

  if (model == "piecewise poisson") {
    if (enroll_fit$accrualTime[1] != 0) {
      stop("accrualTime must start with 0 in enroll_fit")
    }
    if (length(enroll_fit$accrualTime) > 1 &&
        any(diff(enroll_fit$accrualTime) <= 0)) {
      stop("accrualTime should be increasing in enroll_fit")
    }
  }



  erify::check_n(lags, zero = TRUE)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_n(nreps)
  erify::check_bool(showsummary)
  erify::check_bool(showplot)
  erify::check_bool(by_treatment)
  erify::check_n(ngroups)

  if (is.null(df)) by_treatment = TRUE

  if (by_treatment) {
    if (!is.null(df)) {
      if (!("treatment_description" %in% names(df))) {
        df <- df %>% dplyr::mutate(
          treatment_description = paste0("Treatment ", .data$treatment))
      }

      treatment_mapping <- df %>%
        dplyr::select(.data$treatment, .data$treatment_description) %>%
        dplyr::arrange(.data$treatment) %>%
        dplyr::group_by(.data$treatment) %>%
        dplyr::slice(dplyr::n())

      ngroups = nrow(treatment_mapping)
      treatment_label = treatment_mapping$treatment_description
    } else if (!is.null(treatment_label)) {
      treatment_mapping <- dplyr::tibble(
        treatment = 1:ngroups, treatment_description = treatment_label)
    } else {
      treatment_mapping <- dplyr::tibble(
        treatment = 1:ngroups,
        treatment_description = paste0("Treatment ", .data$treatment))

      treatment_label = treatment_mapping$treatment_description
    }

    if (is.null(alloc)) {
      alloc = rep(1, ngroups)
    } else {
      if (length(alloc) != ngroups) {
        stop("length of alloc must be equal to the number of treatments")
      }

      if (any(alloc <= 0 | alloc != round(alloc))) {
        stop("elements of alloc must be positive integers")
      }
    }
  } else {
    ngroups = 1
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }

  if (!is.null(treatment_label) && length(treatment_label) != ngroups) {
    stop(paste("length of treatment_label must be equal to",
               "the number of treatments"))
  }


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

    if (any(df$randdt < trialsdt)) {
      stop("randdt must be greater than or equal to trialsdt")
    }

    if (any(df$randdt > cutoffdt)) {
      stop("randdt must be less than or equal to cutoffdt")
    }

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
        # so that tmax > t
        ti = log(2)/delta  # obtained by setting lambda(ti) = mu/(2*delta)
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
      stop("B-spline enrollment model cannot be used at the design stage")
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

    # mu(t[j]) - mu(t[j-1]) is standard exponential distribution, t[0]=t0
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

  # assign usubjid for new subjects
  newSubjects$usubjid <- rep(paste0("Z-", 100000 + (1:n1)), nreps)
  newSubjects$arrivalTime <- pmax(round(newSubjects$arrivalTime), 1)


  if (by_treatment) {
    # add treatment group information
    blocksize = sum(alloc)
    nblocks = ceiling(n1/blocksize)
    m = nblocks*blocksize
    treats = rep(1:ngroups, alloc)
    index = rep(1:n1, nreps) + rep((0:(nreps-1))*m, each=n1)
    newSubjects$treatment = c(replicate(nreps*nblocks, sample(treats)))[index]

    # summary of observed data by treatment
    if (!is.null(df)) {
      newSubjects <- newSubjects %>%
        dplyr::left_join(treatment_mapping, by = "treatment")

      # add overall treatment
      df2 <- df %>%
        dplyr::bind_rows(df %>% dplyr::mutate(
          treatment = 9999, treatment_description = "Overall"))

      sum_by_trt <- df2 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::summarise(n0 = dplyr::n(), .groups = "drop")
    } else if (!is.null(treatment_label)) {
        newSubjects <- newSubjects %>%
          dplyr::left_join(treatment_mapping, by = "treatment")

        sum_by_trt <- dplyr::tibble(
          treatment = c(1:ngroups, 9999),
          treatment_description = c(treatment_label, "Overall"),
          n0 = 0)
      } else {
        newSubjects <- newSubjects %>% dplyr::mutate(
          treatment_description = paste0("Treatment ", .data$treatment))

        sum_by_trt <- dplyr::tibble(treatment = c(1:ngroups, 9999)) %>%
          dplyr::mutate(treatment_description = c(
            paste0("Treatment ", 1:ngroups), "Overall"),
            n0 = 0)
    }
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

  if (!is.null(df)) {
    pred_date <- as.Date(pred_day - 1, origin = trialsdt)

    str1 <- paste0("Time from cutoff until ", target_n, " subjects: ",
                   pred_date[1] - cutoffdt + 1, " days")
    str2 <- paste0("Median prediction date: ", pred_date[1])
    str3 <- paste0("Prediction interval: ", pred_date[2], ", ", pred_date[3])
    s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")
  } else {
    str1 <- paste0("Time from trial start until ", target_n, " subjects")
    str2 <- paste0("Median prediction day: ", pred_day[1])
    str3 <- paste0("Prediction interval: ", pred_day[2], ", ", pred_day[3])
    s1 <- paste0(str1, "\n", str2, "\n", str3, "\n")
  }


  if (!by_treatment) {
    # predicted number of subjects enrolled after data cut
    dfb1 <- dplyr::tibble(t = t) %>%
      dplyr::cross_join(newSubjects) %>%
      dplyr::group_by(.data$t, .data$draw) %>%
      dplyr::summarise(nenrolled = sum(.data$arrivalTime <= .data$t) + n0,
                       .groups = "drop_last") %>%
      dplyr::summarise(n = quantile(.data$nenrolled, probs = 0.5),
                       pilevel = pilevel,
                       lower = quantile(.data$nenrolled, probs = plower),
                       upper = quantile(.data$nenrolled, probs = pupper),
                       mean = mean(.data$nenrolled),
                       var = var(.data$nenrolled)) %>%
      dplyr::ungroup()

    if (!is.null(df)) {
      # day 1
      df0 <- dplyr::tibble(t = 1, n = 0, pilevel = pilevel,
                           lower = NA, upper = NA, mean = 0, var = 0)

      # arrival time for subjects already enrolled before data cut
      dfa1 <- df %>%
        dplyr::mutate(pilevel = pilevel, lower = NA, upper = NA,
                      mean = .data$n, var = 0) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::bind_rows(dplyr::tibble(t = t0, n = n0, pilevel = pilevel,
                                       lower = NA, upper = NA,
                                       mean = n0, var= 0)) %>%
        dplyr::select(.data$t, .data$n, .data$pilevel, .data$lower,
                      .data$upper, .data$mean, .data$var) %>%
        dplyr::group_by(.data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()


      # concatenate subjects enrolled before and after data cut
      dfs <- dfa1 %>%
        dplyr::bind_rows(dfb1) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::arrange(.data$t)

      # separate data into observed and predicted
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      # plot the enrollment data with date as x-axis
      g1 <- plotly::plot_ly() %>%
        plotly::add_ribbons(
          data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
          name = "prediction interval",
          fill = "tonexty", line = list(width=0)) %>%
        plotly::add_lines(
          data = dfb, x = ~date, y = ~n, name = "median prediction",
          line = list(width=2)) %>%
        plotly::add_lines(
          data = dfa, x = ~date, y = ~n, name = "observed",
          line = list(shape="hv", width=2)) %>%
        plotly::add_lines(
          x = rep(cutoffdt, 2), y = c(min(dfa$n), max(dfb$upper)),
          name = "cutoff", line = list(dash="dash"),
          showlegend = FALSE) %>%
        plotly::layout(
          annotations = list(
            x = cutoffdt, y = 0, text = 'cutoff', xanchor = "left",
            yanchor = "bottom", font = list(size=12),
            showarrow = FALSE),
          xaxis = list(title = "", zeroline = FALSE),
          yaxis = list(title = "Subjects", zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = "h"))
    } else {
      # plot the enrollment data with day as x-axis
      g1 <- plotly::plot_ly(dfb1, x = ~t) %>%
        plotly::add_ribbons(
          ymin = ~lower, ymax = ~upper, name = "prediction interval",
          fill = "tonexty", line = list(width=0)) %>%
        plotly::add_lines(
          y = ~n, name = "median prediction", line = list(width=2)) %>%
        plotly::layout(
          xaxis = list(title = "Days since trial start", zeroline = FALSE),
          yaxis = list(title = "Subjects", zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = "h"))
    }

  } else { # by treatment
    # add overall treatment
    newSubjects2 <- newSubjects %>%
      dplyr::bind_rows(newSubjects %>% dplyr::mutate(
        treatment = 9999, treatment_description = "Overall"))

    # predicted number of subjects enrolled by treatment after cutoff
    dfb1 <- dplyr::tibble(t = t) %>%
      dplyr::cross_join(newSubjects2) %>%
      dplyr::group_by(.data$treatment, .data$treatment_description,
                      .data$t, .data$draw) %>%
      dplyr::summarise(nenrolled = sum(.data$arrivalTime <= .data$t),
                       .groups = "drop_last") %>%
      dplyr::left_join(sum_by_trt,
                       by = c("treatment", "treatment_description")) %>%
      dplyr::mutate(nenrolled = .data$nenrolled + .data$n0) %>%
      dplyr::summarise(n = quantile(.data$nenrolled, probs = 0.5),
                       pilevel = pilevel,
                       lower = quantile(.data$nenrolled, probs = plower),
                       upper = quantile(.data$nenrolled, probs = pupper),
                       mean = mean(.data$nenrolled),
                       var = var(.data$nenrolled),
                       .groups = "drop_last") %>%
      dplyr::ungroup()

    if (!is.null(df)) {
      # day 1
      df0 <- sum_by_trt %>%
        dplyr::select(.data$treatment, .data$treatment_description) %>%
        dplyr::mutate(t = 1, n = 0, pilevel = pilevel, lower = NA,
                      upper = NA, mean = 0, var = 0)

      # arrival time for subjects already enrolled before data cut
      dfa1 <- df2 %>%
        dplyr::group_by(.data$treatment, .data$treatment_description) %>%
        dplyr::arrange(.data$randdt) %>%
        dplyr::mutate(t = as.numeric(.data$randdt - trialsdt + 1),
                      n = dplyr::row_number()) %>%
        dplyr::mutate(pilevel = pilevel, lower = NA, upper = NA,
                      mean = .data$n, var = 0) %>%
        dplyr::bind_rows(df0) %>%
        dplyr::bind_rows(sum_by_trt %>%
                           dplyr::mutate(t = t0, n = n0, pilevel = pilevel,
                                         lower = NA, upper = NA,
                                         mean = n0, var = 0)) %>%
        dplyr::select(.data$treatment, .data$treatment_description,
                      .data$t, .data$n, .data$pilevel, .data$lower,
                      .data$upper, .data$mean, .data$var) %>%
        dplyr::group_by(.data$treatment, .data$treatment_description,
                        .data$t) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::ungroup()

      # concatenate subjects enrolled before and after data cut
      dfs <- dfa1 %>%
        dplyr::bind_rows(dfb1) %>%
        dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
        dplyr::mutate(parameter = "Enrollment") %>%
        dplyr::arrange(.data$treatment, .data$t)

      # separate data into observed and predicted
      dfa <- dfs %>% dplyr::filter(is.na(.data$lower))
      dfb <- dfs %>% dplyr::filter(!is.na(.data$lower))

      g <- list()
      for (i in c(9999, 1:ngroups)) {
        dfsi <- dfs %>%
          dplyr::filter(.data$treatment == i)
        dfbi <- dfb %>%
          dplyr::filter(.data$treatment == i)
        dfai <- dfa %>%
          dplyr::filter(.data$treatment == i)

        g[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
          plotly::add_ribbons(
            data = dfbi, x = ~date, ymin = ~lower, ymax = ~upper,
            name = "prediction interval", fill = "tonexty",
            line = list(width=0)) %>%
          plotly::add_lines(
            data = dfbi, x = ~date, y = ~n, name = "median prediction",
            line = list(width=2)) %>%
          plotly::add_lines(
            data = dfai, x = ~date, y = ~n, name = "observed",
            line = list(shape="hv", width=2)) %>%
          plotly::add_lines(
            x = rep(cutoffdt, 2), y = c(min(dfai$n), max(dfbi$upper)),
            name = "cutoff", line = list(dash="dash"),
            showlegend = FALSE) %>%
          plotly::layout(
            xaxis = list(title = "", zeroline = FALSE),
            yaxis = list(title = "Subjects", zeroline = FALSE),
            legend = list(x = 0, y = 1.05, yanchor = "bottom",
                          orientation = "h")) %>%
          plotly::layout(
            annotations = list(
              x = 0.5, y = 1,
              text = paste0("<b>", dfbi$treatment_description[1], "</b>"),
              xanchor = "center", yanchor = "bottom",
              showarrow = FALSE, xref='paper', yref='paper'))

        if (i == 9999) {
          g[[1]] <- g[[1]] %>%
            plotly::layout(
              annotations = list(
                x = cutoffdt, y = 0, text = 'cutoff', xanchor = "left",
                yanchor = "bottom", font = list(size=12),
                showarrow = FALSE))
        }
      }
    } else { # prediction at design stage
      g <- list()
      for (i in c(9999, 1:ngroups)) {
        dfbi <- dfb1 %>%
          dplyr::filter(.data$treatment == i)

        g[[(i+1) %% 9999]] <- dfbi %>%
          plotly::plot_ly(x = ~t) %>%
          plotly::add_ribbons(
            ymin = ~lower, ymax = ~upper, name = "prediction interval",
            fill = "tonexty", line = list(width=0)) %>%
          plotly::add_lines(
            y = ~n, name = "median prediction", line = list(width=2)) %>%
          plotly::layout(
            xaxis = list(title = "Days since trial start", zeroline = FALSE),
            yaxis = list(title = "Subjects", zeroline = FALSE),
            legend = list(x = 0, y = 1.05, yanchor = "bottom",
                          orientation = "h")) %>%
          plotly::layout(
            annotations = list(
              x = 0.5, y = 1,
              text = paste0("<b>", dfbi$treatment_description[1], "</b>"),
              xanchor = "center", yanchor = "bottom",
              showarrow = FALSE, xref='paper', yref='paper'))
      }
    }

    g1 <- plotly::subplot(g, nrows = ngroups + 1, margin = 0.05)
  }


  if (showsummary) cat(s1)
  if (showplot) print(g1)

  if (!is.null(df)) {
    list(target_n = target_n, enroll_pred_day = pred_day,
         enroll_pred_date = pred_date,
         pilevel = pilevel, nyears = nyears, nreps = nreps,
         newSubjects = newSubjects,
         enroll_pred_df = dfs,
         enroll_pred_summary = s1, enroll_pred_plot = g1)
  } else {
    list(target_n = target_n, enroll_pred_day = pred_day,
         pilevel = pilevel, nyears = nyears, nreps = nreps,
         newSubjects = newSubjects,
         enroll_pred_df = dfb1,
         enroll_pred_summary = s1, enroll_pred_plot = g1)
  }

}
