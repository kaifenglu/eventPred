#' @title Fit time-to-event model
#' @description Fits a specified time-to-event model to the event data.
#'
#' @param df The subject-level event data, including \code{time}
#'   and \code{event}. The data should also include \code{treatment}
#'   coded as 1, 2, and so on, and \code{treatment_description}
#'   for fitting the event model by treatment.
#' @param event_model The event model used to analyze the event data
#'   which can be set to one of the following options:
#'   "exponential", "Weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", "model averaging", or "spline".
#'   The model averaging uses the \code{exp(-bic/2)} weighting and
#'   combines Weibull and log-normal models. The spline model of
#'   Royston and Parmar (2002) assumes that a transformation of
#'   the survival function is modeled as a natural cubic spline
#'   function of log time. By default, it is set to "model averaging".
#' @param piecewiseSurvivalTime A vector that specifies the time
#'   intervals for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param k The number of inner knots of the spline. The default
#'   \code{k=0} gives a Weibull, log-logistic or log-normal model,
#'   if \code{scale} is "hazard", "odds", or "normal", respectively.
#'   The knots are chosen as equally-spaced quantiles of the log
#'   uncensored survival times. The boundary knots are chosen as the
#'   minimum and maximum log uncensored survival times.
#' @param scale If "hazard", the log cumulative hazard is modeled
#'   as a spline function. If "odds", the log cumulative odds is
#'   modeled as a spline function. If "normal", -qnorm(S(t)) is
#'   modeled as a spline function.
#' @param showplot A Boolean variable to control whether or not to
#'   show the fitted time-to-event survival curve. By default, it is
#'   set to \code{TRUE}.
#' @param by_treatment A Boolean variable to control whether or not to
#'   fit the time-to-event data by treatment group. By default,
#'   it is set to \code{FALSE}.
#'
#'
#' @return
#' A list of results from the model fit including key information
#' such as the event model, \code{model}, the estimated model parameters,
#' \code{theta}, the covariance matrix, \code{vtheta}, as well as the
#' Akaike Information Criterion, \code{aic}, and
#' Bayesian Information Criterion, \code{bic}.
#'
#' If the piecewise exponential model is used, the location
#' of knots used in the model, \code{piecewiseSurvivalTime}, will
#' be included in the list of results.
#'
#' If the model averaging option is chosen, the weight assigned
#' to the Weibull component is indicated by the \code{w1} variable.
#'
#' If the spline option is chosen, the \code{knots} and \code{scale}
#' will be included in the list of results.
#'
#' When fitting the event model by treatment, the outcome is presented
#' as a list of lists, where each list element corresponds to a
#' specific treatment group.
#'
#' The fitted time-to-event survival curve is also returned.
#'
#' @examples
#'
#' event_fit <- fitEvent(df = interimData2,
#'                       event_model = "piecewise exponential",
#'                       piecewiseSurvivalTime = c(0, 180))
#'
#' @export
#'
fitEvent <- function(df, event_model = "model averaging",
                     piecewiseSurvivalTime = 0,
                     k = 0, scale = "hazard",
                     showplot = TRUE, by_treatment = FALSE) {
  erify::check_class(df, "data.frame")

  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "model averaging", "spline"))

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  if (length(piecewiseSurvivalTime) > 1 &
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  erify::check_n(k, zero = TRUE)
  erify::check_content(tolower(scale), c("hazard", "odds", "normal"))

  erify::check_bool(showplot)
  erify::check_bool(by_treatment)

  df <- dplyr::as_tibble(df)
  names(df) <- tolower(names(df))

  if (by_treatment) {
    ngroups = length(table(df$treatment))

    if (!("treatment_description" %in% names(df))) {
      df <- df %>% dplyr::mutate(
        treatment_description = paste0("Treatment ", .data$treatment))
    }
  } else {
    ngroups = 1
    df <- df %>% dplyr::mutate(treatment = 1)
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }


  # fit by treatment group
  event_fit <- list()
  g1 <- list()

  for (i in 1:ngroups) {
    df1 <- df %>% dplyr::filter(.data$treatment == i)

    n0 = nrow(df1)
    d0 = sum(df1$event)
    ex0 = sum(df1$time)

    erify::check_positive(d0, supplement = paste(
      "The number of events must be positive to fit an event model."))

    kmfit <- survival::survfit(survival::Surv(time, event) ~ 1, data = df1)
    kmdf <- dplyr::tibble(time = kmfit$time, surv = kmfit$surv)
    kmdf <- dplyr::tibble(time = 0, surv = 1) %>%
      dplyr::bind_rows(kmdf)

    if (tolower(event_model) == "exponential") {
      # lambda(t) = lambda
      # S(t) = exp(-lambda*t)

      fit2 <- list(model = 'Exponential',
                   theta = log(d0/ex0),
                   vtheta = 1/d0,
                   aic = -2*(-d0 + d0*log(d0/ex0)) + 2,
                   bic = -2*(-d0 + d0*log(d0/ex0)) + log(n0))

      # fitted survival curve
      dffit2 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = pexp(.data$time, rate = exp(fit2$theta), lower.tail = FALSE))
    } else if (tolower(event_model) == "weibull") {
      # lambda(t) = kappa/lambda*(t/lambda)^(kappa-1)
      # S(t) = exp(-(t/lambda)^kappa)

      reg <- survival::survreg(survival::Surv(time, event) ~ 1,
                               data = df1, dist = "weibull")

      # weibull$shape = 1/reg$scale, weibull$scale = exp(reg$coefficients)
      # we define theta = (log(weibull$scale), -log(weibull$shape))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit2 <- list(model = "Weibull",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 4,
                   bic = -2*reg$loglik[1] + 2*log(n0))

      # fitted survival curve
      dffit2 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = pweibull(.data$time, shape = exp(-fit2$theta[2]),
                        scale = exp(fit2$theta[1]), lower.tail = FALSE))
    } else if (tolower(event_model) == "log-logistic") {
      # S(t) = 1/(1 + (t/lambda)^kappa)
      reg <- survival::survreg(survival::Surv(time, event) ~ 1,
                               data = df1, dist = "loglogistic")

      # llogis$shape = 1/reg$scale, llogis$scale = exp(reg$coefficients)
      # we define theta = (log(llogis$scale), -log(llogis$shape))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit2 <- list(model = "Log-logistic",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 4,
                   bic = -2*reg$loglik[1] + 2*log(n0))

      # fitted survival curve
      dffit2 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = plogis(log(.data$time), location = fit2$theta[1],
                      scale = exp(fit2$theta[2]), lower.tail = FALSE))
    } else if (tolower(event_model) == "log-normal") {
      # S(t) = 1 - Phi((log(t) - meanlog)/sdlog)
      reg <- survival::survreg(survival::Surv(time, event) ~ 1,
                               data = df1, dist = "lognormal")

      # we use parameterization theta = (meanlog, log(sdlog))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit2 <- list(model = "Log-normal",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 4,
                   bic = -2*reg$loglik[1] + 2*log(n0))

      # fitted survival curve
      dffit2 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = plnorm(.data$time, meanlog = fit2$theta[1],
                      sdlog = exp(fit2$theta[2]), lower.tail = FALSE))
    } else if (tolower(event_model) == "piecewise exponential") {
      # lambda(t) = lambda[j] for ucut[j] < t <= ucut[j+1], j = 1,...,J
      # where ucut[1]=0< ucut[2]< ...< ucut[J]< ucut[J+1]=Inf are the knots
      u = piecewiseSurvivalTime[piecewiseSurvivalTime < max(df1$time)]
      ucut = c(u, max(df1$time))
      J = length(u)

      d = rep(NA, J)  # number of events in each interval
      ex = rep(NA, J) # total exposure in each interval
      for (j in 1:J) {
        d[j] = sum((df1$time > ucut[j]) * (df1$time <= ucut[j+1]) *
                     (df1$event == 1))
        ex[j] = sum(pmax(0, pmin(df1$time, ucut[j+1]) - ucut[j]))
      }

      # maximum likelihood estimates and covariance matrix
      if (J > 1) {
        vtheta = diag(1/d)
      } else {
        vtheta = 1/d*diag(1)
      }

      fit2 <- list(model = "Piecewise exponential",
                   theta = log(d/ex),
                   vtheta = vtheta,
                   aic = -2*sum(-d + d*log(d/ex)) + 2*J,
                   bic = -2*sum(-d + d*log(d/ex)) + J*log(n0),
                   piecewiseSurvivalTime = u)

      # fitted survival curve
      time = seq(0, max(df1$time))

      lambda = d/ex
      if (J>1) {
        psum = c(0, cumsum(lambda[1:(J-1)] * diff(u)))
      } else {
        psum = 0
      }
      j = findInterval(time, u)
      m = psum[j] + lambda[j]*(time - u[j])
      surv = exp(-m)

      dffit2 <- dplyr::tibble(time, surv)
    } else if (tolower(event_model) == "model averaging") {
      reg1 <- survival::survreg(survival::Surv(time, event) ~ 1,
                                data = df1, dist = "weibull")
      reg2 <- survival::survreg(survival::Surv(time, event) ~ 1,
                                data = df1, dist = "lognormal")
      aic1 <- -2*reg1$loglik[1] + 4
      aic2 <- -2*reg2$loglik[1] + 4
      bic1 <- -2*reg1$loglik[1] + 2*log(n0)
      bic2 <- -2*reg2$loglik[1] + 2*log(n0)

      w1 = 1/(1 + exp(-0.5*(bic2 - bic1)))


      # model parameters from weibull and log-normal
      theta = c(as.numeric(reg1$coefficients), log(reg1$scale),
                as.numeric(reg2$coefficients), log(reg2$scale))

      # variance-covariance matrix, noting that the covariances
      # between the two sets of parameters are zero as they are estimated
      # from different likelihood functions
      vtheta = as.matrix(Matrix::bdiag(reg1$var, reg2$var))

      # model fit, assuming fixed weight w1
      fit2 <- list(model = "Model averaging",
                   theta = theta,
                   vtheta = vtheta,
                   aic = w1*aic1 + (1-w1)*aic2,
                   bic = w1*bic1 + (1-w1)*bic2,
                   w1 = w1)

      # distribution function for model averaging of Weibull and log-normal
      pmodavg <- function(t, theta, w1, lower.tail = TRUE, log.p = FALSE) {
        shape = exp(-theta[2])
        scale = exp(theta[1])
        meanlog = theta[3]
        sdlog = exp(theta[4])

        p1 = pweibull(pmax(0,t), shape, scale)
        p2 = plnorm(pmax(0,t), meanlog, sdlog)
        p = w1*p1 + (1-w1)*p2

        if (!lower.tail) p = 1 - p
        if (log.p) p = log(p)
        p
      }

      dffit2 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = pmodavg(.data$time, theta, w1, lower.tail = FALSE))
    } else if (tolower(event_model) == "spline") {
      # g(S(t)) = gamma_0 +gamma_1*x +gamma_2*v_1(x) +... +gamma_{m+1}*v_m(x)

      spl <- flexsurv::flexsurvspline(survival::Surv(time, event) ~ 1,
                                      data = df1, k = k, scale = scale,
                                      method = "Nelder-Mead")

      fit2 <- list(model = "Spline",
                   theta = spl$coefficients,
                   vtheta = spl$cov,
                   aic = -2*spl$loglik + 2*(k+2),
                   bic = -2*spl$loglik + (k+2)*log(n0),
                   knots = spl$knots,
                   scale = spl$scale)

      # fitted survival curve
      dffit2 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = flexsurv::psurvspline(.data$time, gamma = spl$coefficients,
                                     knots = spl$knots, scale = spl$scale,
                                     lower.tail = FALSE))
    }


    # plot the survival curve
    if (tolower(fit2$model) == "piecewise exponential") {
      modeltext = paste0(paste0(fit2$model, "("),
                         paste(piecewiseSurvivalTime, collapse = " "), ")")
    } else if (tolower(fit2$model) == "spline") {
      modeltext = paste0(fit2$model, "(k = ", k, ", ", "scale = '",
                         scale, "')")
    } else {
      modeltext = fit2$model
    }

    if (tolower(event_model) == "model averaging") {
      aictext = paste("Weighted AIC:", round(fit2$aic,2))
      bictext = paste("Weighted BIC:", round(fit2$bic,2))
    } else {
      aictext = paste("AIC:", round(fit2$aic,2))
      bictext = paste("BIC:", round(fit2$bic,2))
    }

    fittedEvent <- plotly::plot_ly() %>%
      plotly::add_lines(
        data=kmdf, x=~time, y=~surv, name="Kaplan-Meier",
        line=list(shape="hv")) %>%
      plotly::add_lines(
        data=dffit2, x=~time, y=~surv, name="fitted") %>%
      plotly::layout(
        xaxis = list(title = "Days since randomization", zeroline = FALSE),
        yaxis = list(title = "Survival probability", zeroline = FALSE),
        title = list(text = "Fitted time to event survival curve"),
        annotations = list(
          x = c(0.75, 0.75, 0.75), y = c(0.95, 0.80, 0.65), xref = "paper",
          yref = "paper", text = paste('<i>', c(modeltext, aictext,
                                                bictext), '</i>'),
          xanchor = "left", font = list(size = 14, color = "red"),
          showarrow = FALSE)) %>%
      plotly::hide_legend()

    if (by_treatment && ngroups > 1) {
      fittedEvent <- fittedEvent %>%
        plotly::layout(annotations = list(
          x = 0.5, y = 1,
          text = paste0("<b>", df1$treatment_description[1], "</b>"),
          xanchor = "center", yanchor = "middle", showarrow = FALSE,
          xref='paper', yref='paper'))
    }

    if (by_treatment) {
      fit2$treatment = df1$treatment[1]
      fit2$treatment_description = df1$treatment_description[1]
    }

    event_fit[[i]] = fit2
    g1[[i]] = fittedEvent
  }

  if (!by_treatment) {
    event_fit = fit2
    event_fit_plot = fittedEvent
  } else {
    event_fit_plot <- plotly::subplot(g1, nrows = ngroups, titleX = TRUE,
                                      titleY = TRUE, margin = 0.1)
  }

  if (showplot) print(event_fit_plot)

  list(event_fit = event_fit, event_fit_plot = event_fit_plot)
}
