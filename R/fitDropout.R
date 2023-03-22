#' @title Fit time-to-dropout model
#' @description Fits a specified time-to-dropout model to the dropout data.
#'
#' @param df The subject-level dropout data, including \code{time} and
#'   \code{dropout}.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options: "exponential",
#'   "Weibull", or "log-normal". By default, it is set to "Weibull".
#' @param showplot A Boolean variable to control whether or not to
#'   show the fitted time-to-dropout survival curve. By default, it is
#'   set to \code{TRUE}.
#'
#' @return A list of results from the model fit including key information
#'   such as the dropout model, \code{model}, the estimated model parameters,
#'   \code{theta}, the covariance matrix, \code{vtheta}, as well as the
#'   Bayesian Information Criterion, \code{bic}.
#'   The fitted time-to-dropout curve is also returned.
#'
#' @examples
#'
#' dropout_fit <- fitDropout(df = interimData2, dropout_model = "exponential")
#'
#' @export
#'
fitDropout <- function(df, dropout_model = "weibull", showplot = TRUE) {
  erify::check_class(df, "data.frame")

  erify::check_content(tolower(dropout_model),
                       c("exponential", "weibull", "log-normal"))

  df <- dplyr::as_tibble(df)
  names(df) <- tolower(names(df))
  df$randdt <- as.Date(df$randdt)
  df$cutoffdt <- as.Date(df$cutoffdt)

  n0 = nrow(df)
  d0 = sum(df$dropout)
  ex0 = sum(df$time)

  erify::check_positive(d0, supplement = paste(
    "The number of dropouts must be positive to fit a dropout model."))

  kmfit <- survival::survfit(survival::Surv(time, dropout) ~ 1, data = df)
  kmdf <- dplyr::tibble(time = kmfit$time, surv = kmfit$surv)
  kmdf <- dplyr::tibble(time = 0, surv = 1) %>%
    dplyr::bind_rows(kmdf)

  if (tolower(dropout_model) == "exponential") {
    # lambda(t) = lambda
    # S(t) = exp(-lambda*t)

    fit3 <- list(model = 'Exponential',
                 theta = log(d0/ex0),
                 vtheta = 1/d0,
                 bic = -2*(-d0 + d0*log(d0/ex0)) + log(n0))

    # fitted survival curve
    dffit3 <- dplyr::tibble(
      time = seq(0, max(df$time)),
      surv = pexp(.data$time, rate = exp(fit3$theta), lower.tail = FALSE))

  } else if (tolower(dropout_model) == "weibull") {
    # lambda(t) = kappa/lambda*(t/lambda)^(kappa-1)
    # S(t) = exp(-(t/lambda)^kappa)

    reg <- survival::survreg(survival::Surv(time, dropout) ~ 1,
                             data = df, dist = "weibull")

    # Note: weibull$shape = 1/reg$scale, weibull$scale = exp(reg$coefficients)
    # we use parameterization theta = (log(weibull$shape), log(weibull$scale))
    # reg$var is for c(reg$coefficients, log(reg$scale)) = lmat %*% theta
    lmat <- matrix(c(0, -1, 1, 0), nrow=2, ncol=2, byrow=TRUE)
    fit3 <- list(model = "Weibull",
                 theta = c(log(1/reg$scale), as.numeric(reg$coefficients)),
                 vtheta = lmat %*% reg$var %*% t(lmat),
                 bic = -2*reg$loglik[1] + 2*log(n0))

    # fitted survival curve
    dffit3 <- dplyr::tibble(
      time = seq(0, max(df$time)),
      surv = pweibull(.data$time, shape = exp(fit3$theta[1]),
                      scale = exp(fit3$theta[2]), lower.tail = FALSE))
  } else if (tolower(dropout_model) == "log-normal") {
    # S(t) = 1 - Phi((log(t) - meanlog)/sdlog)
    reg <- survival::survreg(survival::Surv(time, dropout) ~ 1,
                             data = df, dist = "lognormal")

    # we use parameterization theta = (meanlog, log(sdlog))
    # reg$var is for c(reg$coefficients, log(reg$scale)) = theta
    fit3 <- list(model = "Log-normal",
                 theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                 vtheta = reg$var,
                 bic = -2*reg$loglik[1] + 2*log(n0))

    # fitted survival curve
    dffit3 <- dplyr::tibble(
      time = seq(0, max(df$time)),
      surv = plnorm(.data$time, meanlog = fit3$theta[1],
                    sdlog = exp(fit3$theta[2]), lower.tail = FALSE))
  }


  bictext = paste("BIC:", round(fit3$bic,2))

  # plot the survival curve
  fittedDropout <- plotly::plot_ly() %>%
    plotly::add_lines(data=kmdf, x=~time, y=~surv, name="Kaplan-Meier",
                      line=list(shape="hv")) %>%
    plotly::add_lines(data=dffit3, x=~time, y=~surv, name="fitted") %>%
    plotly::layout(
      xaxis = list(title = "Days since randomization", zeroline = FALSE),
      yaxis = list(title = "Survival probability"),
      title = list(text = "Fitted time to dropout survival curve"),
      annotations = list(
        x = c(0.75, 0.75), y = c(0.95, 0.88),
        xref = "paper", yref = "paper",
        text = paste('<i>', c(fit3$model, bictext), '</i>'),
        xanchor = "left",
        font = list(size = 14, color = "red"),
        showarrow = FALSE
      )) %>%
    plotly::hide_legend()

  if (showplot) print(fittedDropout)

  list(dropout_fit = fit3, dropout_fit_plot = fittedDropout)

}
