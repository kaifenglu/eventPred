#' @title Fit time-to-event model
#' @description Fits a specified time-to-event model to the event data.
#'
#' @param df The subject-level event data, including \code{time}
#'   and \code{event}.
#' @param event_model The event model used to analyze the event data
#'   which can be set to one of the
#'   following options: "exponential", "Weibull", "log-normal",
#'   "piecewise exponential", or "model averaging". The model averaging
#'   uses the \code{exp(-bic/2)} weighting and combines Weibull and
#'   log-normal models. By default, it is set to "model
#'   averaging".
#' @param piecewiseSurvivalTime A vector that specifies the time
#'   intervals for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param showplot A Boolean variable to control whether or not to
#'   show the fitted time-to-event survival curve. By default, it is
#'   set to \code{TRUE}.
#'
#'
#' @return
#' A list of results from the model fit including key information
#' such as the event model, \code{model}, the estimated model parameters,
#' \code{theta}, the covariance matrix, \code{vtheta}, as well as the
#' Bayesian Information Criterion, \code{bic}.
#'
#' If the piecewise exponential model is used, the location
#' of knots used in the model, \code{piecewiseSurvivalTime}, will
#' be included in the list of results.
#'
#' If the model averaging option is chosen, the weight assigned
#' to the Weibull component is indicated by the \code{w1} variable.
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
                     piecewiseSurvivalTime = 0, showplot = TRUE) {
  erify::check_class(df, "data.frame")

  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-normal",
                         "piecewise exponential", "model averaging"))

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  if (length(piecewiseSurvivalTime) > 1 &
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  erify::check_bool(showplot)


  df <- dplyr::as_tibble(df)
  names(df) <- tolower(names(df))

  n0 = nrow(df)
  d0 = sum(df$event)
  ex0 = sum(df$time)

  erify::check_positive(d0, supplement = paste(
    "The number of events must be positive to fit an event model."))

  kmfit <- survival::survfit(survival::Surv(time, event) ~ 1, data = df)
  kmdf <- dplyr::tibble(time = kmfit$time, surv = kmfit$surv)
  kmdf <- dplyr::tibble(time = 0, surv = 1) %>%
    dplyr::bind_rows(kmdf)

  if (tolower(event_model) == "exponential") {
    # lambda(t) = lambda
    # S(t) = exp(-lambda*t)

    fit2 <- list(model = 'Exponential',
                 theta = log(d0/ex0),
                 vtheta = 1/d0,
                 bic = -2*(-d0 + d0*log(d0/ex0)) + log(n0))

    # fitted survival curve
    dffit2 <- dplyr::tibble(
      time = seq(0, max(df$time)),
      surv = pexp(.data$time, rate = exp(fit2$theta), lower.tail = FALSE))
  } else if (tolower(event_model) == "weibull") {
    # lambda(t) = kappa/lambda*(t/lambda)^(kappa-1)
    # S(t) = exp(-(t/lambda)^kappa)

    reg <- survival::survreg(survival::Surv(time, event) ~ 1,
                             data = df, dist = "weibull")

    # Note: weibull$shape = 1/reg$scale, weibull$scale = exp(reg$coefficients)
    # we use parameterization theta = (log(weibull$shape), log(weibull$scale))
    # reg$var is for c(reg$coefficients, log(reg$scale))
    lmat <- matrix(c(0, -1, 1, 0), nrow=2, ncol=2, byrow=TRUE)
    fit2 <- list(model = "Weibull",
                 theta = c(log(1/reg$scale), as.numeric(reg$coefficients)),
                 vtheta = lmat %*% reg$var %*% t(lmat),
                 bic = -2*reg$loglik[1] + 2*log(n0))

    # fitted survival curve
    dffit2 <- dplyr::tibble(
      time = seq(0, max(df$time)),
      surv = pweibull(.data$time, shape = exp(fit2$theta[1]),
                      scale = exp(fit2$theta[2]), lower.tail = FALSE))
  } else if (tolower(event_model) == "log-normal") {
    # S(t) = 1 - Phi((log(t) - meanlog)/sdlog)
    reg <- survival::survreg(survival::Surv(time, event) ~ 1,
                             data = df, dist = "lognormal")

    # we use parameterization theta = (meanlog, log(sdlog))
    # reg$var is for c(reg$coefficients, log(reg$scale)) = theta
    fit2 <- list(model = "Log-normal",
                 theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                 vtheta = reg$var,
                 bic = -2*reg$loglik[1] + 2*log(n0))

    # fitted survival curve
    dffit2 <- dplyr::tibble(
      time = seq(0, max(df$time)),
      surv = plnorm(.data$time, meanlog = fit2$theta[1],
                    sdlog = exp(fit2$theta[2]), lower.tail = FALSE))
  } else if (tolower(event_model) == "piecewise exponential") {
    # lambda(t) = lambda[j] for ucut[j] < t <= ucut[j+1], j = 1,...,J
    # where ucut[1]=0 < ucut[2] < ... < ucut[J] < ucut[J+1]=Inf are the knots
    u = piecewiseSurvivalTime[piecewiseSurvivalTime < max(df$time)]
    ucut = c(u, max(df$time))
    J = length(u)

    d = rep(NA, J)  # number of events in each interval
    ex = rep(NA, J) # total exposure in each interval
    for (j in 1:J) {
      d[j] = sum((df$time > ucut[j]) * (df$time <= ucut[j+1]) *
                   (df$event == 1))
      ex[j] = sum(pmax(0, pmin(df$time, ucut[j+1]) - ucut[j]))
    }

    # maximum likelihood estimates and covariance matrix
    if (J > 1) {
      fit2 <- list(model = "Piecewise exponential",
                   theta = log(d/ex),
                   vtheta = diag(1/d),
                   bic = -2*sum(-d + d*log(d/ex)) + J*log(n0),
                   piecewiseSurvivalTime = u)
    } else {
      fit2 <- list(model = "Piecewise exponential",
                   theta = log(d/ex),
                   vtheta = 1/d*diag(1),
                   bic = -2*sum(-d + d*log(d/ex)) + J*log(n0),
                   piecewiseSurvivalTime = u)
    }

    # fitted survival curve
    time = seq(0, max(df$time))

    surv = 0
    for (j in 1:J) {
      exj = pmax(0, pmin(time, ucut[j+1]) - ucut[j])
      surv = surv + exp(fit2$theta[j]) * exj
    }
    surv = exp(-surv)

    dffit2 <- dplyr::tibble(time, surv)
  } else if (tolower(event_model) == "model averaging") {
    reg1 <- survival::survreg(survival::Surv(time, event) ~ 1,
                              data = df, dist = "weibull")
    reg2 <- survival::survreg(survival::Surv(time, event) ~ 1,
                              data = df, dist = "lognormal")
    bic1 <- -2*reg1$loglik[1] + 2*log(n0)
    bic2 <- -2*reg2$loglik[1] + 2*log(n0)

    w1 = 1/(1 + exp(-0.5*(bic2 - bic1)))


    # model parameters from weibull and log-normal
    theta = c(log(1/reg1$scale), as.numeric(reg1$coefficients),
              as.numeric(reg2$coefficients), log(reg2$scale))

    # variance-covariance matrix, noting that the covariances
    # between the two sets of parameters are zero as they are estimated
    # from different likelihood functions
    lmat <- matrix(c(0, -1, 1, 0), nrow=2, ncol=2, byrow=TRUE)
    vtheta = as.matrix(Matrix::bdiag(lmat %*% reg1$var %*% t(lmat), reg2$var))

    # model fit, assuming fixed weight w1
    fit2 <- list(model = "Model averaging",
                 theta = theta,
                 vtheta = vtheta,
                 bic = w1*bic1 + (1-w1)*bic2,
                 w1 = w1)

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

    dffit2 <- dplyr::tibble(
      time = seq(0, max(df$time)),
      surv = pmodavg(.data$time, theta, w1, lower.tail = FALSE))
  }


  # plot the survival curve
  if (tolower(event_model) == "model averaging") {
    bictext = paste("Weighted BIC:", round(fit2$bic,2))
  } else {
    bictext = paste("BIC:", round(fit2$bic,2))
  }

  fittedEvent <- plotly::plot_ly() %>%
    plotly::add_lines(data=kmdf, x=~time, y=~surv, name="Kaplan-Meier",
                      line=list(shape="hv")) %>%
    plotly::add_lines(data=dffit2, x=~time, y=~surv, name="fitted") %>%
    plotly::layout(
      xaxis = list(title = "Days since randomization", zeroline = FALSE),
      yaxis = list(title = "Survival probability", zeroline = FALSE),
      title = list(text = "Fitted time to event survival curve"),
      annotations = list(
        x = c(0.75, 0.75), y = c(0.95, 0.85),
        xref = "paper", yref = "paper",
        text = paste('<i>', c(fit2$model, bictext), '</i>'),
        xanchor = "left",
        font = list(size = 14, color = "red"),
        showarrow = FALSE
      )) %>%
    plotly::hide_legend()

  if (showplot) print(fittedEvent)

  list(event_fit = fit2, event_fit_plot = fittedEvent)
}
