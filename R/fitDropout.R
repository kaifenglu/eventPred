#' @title Fit time-to-dropout model
#' @description Fits a specified time-to-dropout model to the dropout data.
#'
#' @param df The subject-level dropout data, including \code{time} and
#'   \code{dropout}. The data should also include \code{treatment}
#'   coded as 1, 2, and so on, and \code{treatment_description}
#'   for fitting the dropout model by treatment.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options: "exponential",
#'   "Weibull", "log-logistic", "log-normal", "piecewise exponential",
#'   or "spline". The spline model of Royston and Parmer (2022) assumes
#'   that a transformation of the survival function is modeled as a
#'   natural cubic spline function of log time.
#'   By default, it is set to "exponential".
#' @param piecewiseDropoutTime A vector that specifies the time
#'   intervals for the piecewise exponential dropout distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param k_dropout The number of inner knots of the spline. The default
#'   \code{k_dropout=0} gives a Weibull, log-logistic or log-normal model,
#'   if \code{scale_dropout} is "hazard", "odds", or "normal", respectively.
#'   The knots are chosen as equally-spaced quantiles of the log
#'   uncensored survival times. The boundary knots are chosen as the
#'   minimum and maximum log uncensored survival times.
#' @param scale_dropout If "hazard", the log cumulative hazard is modeled
#'   as a spline function. If "odds", the log cumulative odds is
#'   modeled as a spline function. If "normal", -qnorm(S(t)) is
#'   modeled as a spline function.
#' @param showplot A Boolean variable to control whether or not to
#'   show the fitted time-to-dropout survival curve. By default, it is
#'   set to \code{TRUE}.
#' @param by_treatment A Boolean variable to control whether or not to
#'   fit the time-to-dropout data by treatment group. By default,
#'   it is set to \code{FALSE}.
#'
#' @return A list of results from the model fit including key information
#' such as the dropout model, \code{model}, the estimated model parameters,
#' \code{theta}, the covariance matrix, \code{vtheta}, as well as the
#' Akaike Information Criterion, \code{aic},
#' and Bayesian Information Criterion, \code{bic}.
#'
#' If the piecewise exponential model is used, the location
#' of knots used in the model, \code{piecewiseDropoutTime}, will
#' be included in the list of results.
#'
#' If the spline option is chosen, the \code{knots} and
#' \code{scale} will be included in the list of results.
#'
#' When fitting the dropout model by treatment, the outcome is presented
#' as a list of lists, where each list element corresponds to a
#' specific treatment group.
#'
#' The fitted time-to-dropout survival curve is also returned.
#'
#' @examples
#'
#' dropout_fit <- fitDropout(df = interimData2, dropout_model = "exponential")
#'
#' @export
#'
fitDropout <- function(df, dropout_model = "exponential",
                       piecewiseDropoutTime = 0,
                       k_dropout = 0, scale_dropout = "hazard",
                       showplot = TRUE, by_treatment = FALSE) {
  erify::check_class(df, "data.frame")

  erify::check_content(tolower(dropout_model),
                       c("exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "spline"))

  if (piecewiseDropoutTime[1] != 0) {
    stop("piecewiseDropoutTime must start with 0");
  }
  if (length(piecewiseDropoutTime) > 1 &
      any(diff(piecewiseDropoutTime) <= 0)) {
    stop("piecewiseDropoutTime should be increasing")
  }

  erify::check_n(k_dropout, zero = TRUE)
  erify::check_content(tolower(scale_dropout), c("hazard", "odds", "normal"))

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
  dropout_fit <- list()
  g1 <- list()

  for (i in 1:ngroups) {
    df1 <- df %>% dplyr::filter(.data$treatment == i)

    n0 = nrow(df1)
    c0 = sum(df1$dropout)
    ex0 = sum(df1$time)

    kmfit <- survival::survfit(survival::Surv(time, dropout) ~ 1, data = df1)
    kmdf <- dplyr::tibble(time = kmfit$time, surv = kmfit$surv)
    kmdf <- dplyr::tibble(time = 0, surv = 1) %>%
      dplyr::bind_rows(kmdf)

    if (tolower(dropout_model) == "exponential") {
      erify::check_positive(c0, supplement = paste(
        "The number of dropouts must be >= 1 to fit an exponential model."))

      # lambda(t) = lambda
      # S(t) = exp(-lambda*t)

      fit3 <- list(model = 'Exponential',
                   theta = log(c0/ex0),
                   vtheta = 1/c0,
                   aic = -2*(-c0 + c0*log(c0/ex0)) + 2,
                   bic = -2*(-c0 + c0*log(c0/ex0)) + log(n0))

      # fitted survival curve
      dffit3 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = pexp(.data$time, rate = exp(fit3$theta), lower.tail = FALSE))
    } else if (tolower(dropout_model) == "weibull") {
      erify::check_positive(c0 - 1, supplement = paste(
        "The number of dropouts must be >= 2 to fit a Weibull model."))

      # lambda(t) = kappa/lambda*(t/lambda)^(kappa-1)
      # S(t) = exp(-(t/lambda)^kappa)

      reg <- survival::survreg(survival::Surv(time, dropout) ~ 1,
                               data = df1, dist = "weibull")

      # weibull$shape = 1/reg$scale, weibull$scale = exp(reg$coefficients)
      # we define theta = c(log(weibull$scale), -log(weibull$shape))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit3 <- list(model = "Weibull",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 4,
                   bic = -2*reg$loglik[1] + 2*log(n0))

      # fitted survival curve
      dffit3 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = pweibull(.data$time, shape = exp(-fit3$theta[2]),
                        scale = exp(fit3$theta[1]), lower.tail = FALSE))
    } else if (tolower(dropout_model) == "log-logistic") {
      erify::check_positive(c0 - 1, supplement = paste(
        "The number of dropouts must be >= 2 to fit a log-logistic model."))

      # S(t) = 1/(1 + (t/lambda)^kappa)

      reg <- survival::survreg(survival::Surv(time, dropout) ~ 1,
                               data = df1, dist = "loglogistic")

      # llogis$shape = 1/reg$scale, llogis$scale = exp(reg$coefficients)
      # we define theta = (log(llogis$scale), -log(llogis$shape))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit3 <- list(model = "Log-logistic",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 4,
                   bic = -2*reg$loglik[1] + 2*log(n0))

      # fitted survival curve
      dffit3 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = plogis(log(.data$time), location = fit3$theta[1],
                      scale = exp(fit3$theta[2]), lower.tail = FALSE))
    } else if (tolower(dropout_model) == "log-normal") {
      erify::check_positive(c0 - 1, supplement = paste(
        "The number of dropouts must be >= 2 to fit a log-normal model."))

      # S(t) = 1 - Phi((log(t) - meanlog)/sdlog)

      reg <- survival::survreg(survival::Surv(time, dropout) ~ 1,
                               data = df1, dist = "lognormal")

      # we use parameterization theta = (meanlog, log(sdlog))
      # reg$var is for c(reg$coefficients, log(reg$scale)) = theta
      fit3 <- list(model = "Log-normal",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 4,
                   bic = -2*reg$loglik[1] + 2*log(n0))

      # fitted survival curve
      dffit3 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = plnorm(.data$time, meanlog = fit3$theta[1],
                      sdlog = exp(fit3$theta[2]), lower.tail = FALSE))
    } else if (tolower(dropout_model) == "piecewise exponential") {
      # lambda(t) = lambda[j] for ucut[j] < t <= ucut[j+1], j = 1,...,J
      # where ucut[1]=0< ucut[2]< ...< ucut[J]< ucut[J+1]=Inf are the knots

      u = piecewiseDropoutTime[piecewiseDropoutTime < max(df1$time)]
      ucut = c(u, max(df1$time))
      J = length(u)

      d = rep(NA, J)  # number of events in each interval
      ex = rep(NA, J) # total exposure in each interval
      for (j in 1:J) {
        d[j] = sum((df1$time > ucut[j]) * (df1$time <= ucut[j+1]) *
                     (df1$dropout == 1))
        ex[j] = sum(pmax(0, pmin(df1$time, ucut[j+1]) - ucut[j]))
      }

      if (any(d == 0)) {
        stop(paste("The number of dropouts must be >= 1 in each interval",
                   "to fit a piecewise exponential model."))
      }

      # maximum likelihood estimates and covariance matrix
      if (J > 1) {
        vtheta = diag(1/d)
      } else {
        vtheta = 1/d*diag(1)
      }

      fit3 <- list(model = "Piecewise exponential",
                   theta = log(d/ex),
                   vtheta = vtheta,
                   aic = -2*sum(-d + d*log(d/ex)) + 2*J,
                   bic = -2*sum(-d + d*log(d/ex)) + J*log(n0),
                   piecewiseDropoutTime = u)

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

      dffit3 <- dplyr::tibble(time, surv)
    } else if (tolower(dropout_model) == "spline") {
      erify::check_positive(c0 - k_dropout - 1, supplement = paste(
        "The number of dropouts must be >=", k_dropout + 2,
        "to fit a spline model."))

      # g(S(t)) = gamma_0 +gamma_1*x +gamma_2*v_1(x) +... +gamma_{m+1}*v_m(x)

      spl <- flexsurv::flexsurvspline(survival::Surv(time, dropout) ~ 1,
                                      data = df1, k = k_dropout,
                                      scale = scale_dropout,
                                      method = "Nelder-Mead")

      fit3 <- list(model = "Spline",
                   theta = spl$coefficients,
                   vtheta = spl$cov,
                   aic = -2*spl$loglik + 2*(k_dropout+2),
                   bic = -2*spl$loglik + (k_dropout+2)*log(n0),
                   knots = spl$knots,
                   scale = spl$scale)

      # fitted survival curve
      dffit3 <- dplyr::tibble(
        time = seq(0, max(df1$time)),
        surv = flexsurv::psurvspline(.data$time, gamma = spl$coefficients,
                                     knots = spl$knots, scale = spl$scale,
                                     lower.tail = FALSE))
    }



    # plot the survival curve
    if (tolower(fit3$model) == "piecewise exponential") {
      modeltext = paste0(paste0(fit3$model, "("),
                         paste(piecewiseDropoutTime, collapse = " "), ")")
    } else if (tolower(fit3$model) == "spline") {
      modeltext = paste0(fit3$model, "(k = ", k_dropout, ", ", "scale = '",
                         scale_dropout, "')")
    } else {
      modeltext = fit3$model
    }

    aictext = paste("AIC:", round(fit3$aic,2))
    bictext = paste("BIC:", round(fit3$bic,2))

    # plot the survival curve
    fittedDropout <- plotly::plot_ly() %>%
      plotly::add_lines(
        data=kmdf, x=~time, y=~surv, name="Kaplan-Meier",
        line=list(shape="hv")) %>%
      plotly::add_lines(
        data=dffit3, x=~time, y=~surv, name="fitted") %>%
      plotly::layout(
        xaxis = list(title = "Days since randomization", zeroline = FALSE),
        yaxis = list(title = "Survival probability", zeroline = FALSE),
        title = list(text = "Fitted time to dropout survival curve"),
        annotations = list(
          x = c(0.75, 0.75, 0.75), y = c(0.95, 0.80, 0.65), xref = "paper",
          yref = "paper", text = paste('<i>', c(modeltext, aictext,
                                                bictext), '</i>'),
          xanchor = "left", font = list(size = 14, color = "red"),
          showarrow = FALSE)) %>%
      plotly::hide_legend()

    if (by_treatment && ngroups > 1) {
      fittedDropout <- fittedDropout %>%
        plotly::layout(annotations = list(
          x = 0.5, y = 1,
          text = paste0("<b>", df1$treatment_description[1], "</b>"),
          xanchor = "center", yanchor = "middle", showarrow = FALSE,
          xref='paper', yref='paper'))
    }

    if (by_treatment) {
      fit3$treatment = df1$treatment[1]
      fit3$treatment_description = df1$treatment_description[1]
    }

    dropout_fit[[i]] = fit3
    g1[[i]] = fittedDropout
  }

  if (!by_treatment) {
    dropout_fit = fit3
    dropout_fit_plot = fittedDropout
  } else {
    dropout_fit_plot <- plotly::subplot(g1, nrows = ngroups, titleX = TRUE,
                                        titleY = TRUE, margin = 0.1)
  }

  if (showplot) print(dropout_fit_plot)

  list(dropout_fit = dropout_fit, dropout_fit_plot = dropout_fit_plot)
}
