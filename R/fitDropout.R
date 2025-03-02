#' @title Fit time-to-dropout model
#' @description Fits a specified time-to-dropout model to the dropout data.
#'
#' @param df The subject-level dropout data, including \code{time} and
#'   \code{dropout}. The data should also include \code{treatment}
#'   coded as 1, 2, and so on, and \code{treatment_description}
#'   for fitting the dropout model by treatment.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options:
#'   "exponential", "Weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", "model averaging", or "spline".
#'   The model averaging uses the \code{exp(-bic/2)} weighting and
#'   combines Weibull and log-normal models. The spline model of
#'   Royston and Parmar (2002) assumes that a transformation of
#'   the survival function is modeled as a natural cubic spline
#'   function of log time. By default, it is set to "exponential".
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
#' @param covariates The names of baseline covariates from the input
#'   data frame to include in the dropout model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#'
#' @return A list of results from the model fit including key information
#' such as the dropout model, \code{model}, the estimated model parameters,
#' \code{theta}, the covariance matrix, \code{vtheta}, as well as the
#' Akaike Information Criterion, \code{aic}, and
#' Bayesian Information Criterion, \code{bic}.
#'
#' If the piecewise exponential model is used, the location
#' of knots used in the model, \code{piecewiseDropoutTime}, will
#' be included in the list of results.
#'
#' If the model averaging option is chosen, the weight assigned
#' to the Weibull component is indicated by the \code{w1} variable.
#'
#' If the spline option is chosen, the \code{knots} and \code{scale}
#' will be included in the list of results.
#'
#' When fitting the dropout model by treatment, the outcome is presented
#' as a list of lists, where each list element corresponds to a
#' specific treatment group.
#'
#' The fitted time-to-dropout survival curve is also returned.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Patrick Royston and Mahesh K. B. Parmar. Flexible parametric
#' proportional-hazards and proportional-odds models for censored
#' survival data, with application to prognostic modelling and
#' estimation of treatment effects. Stat in Med. 2002; 21:2175-2197.
#'
#' @examples
#'
#' dropout_fit <- fitDropout(df = interimData2,
#'                           dropout_model = "exponential")
#'
#' @export
#'
fitDropout <- function(df, dropout_model = "exponential",
                       piecewiseDropoutTime = 0,
                       k_dropout = 0, scale_dropout = "hazard",
                       showplot = TRUE, by_treatment = FALSE,
                       covariates = NULL) {

  erify::check_class(df, "data.frame")

  erify::check_content(tolower(dropout_model),
                       c("exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "model averaging", "spline"))

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

  # construct the formula for survival analysis
  if (!is.null(covariates)) {
    if (!all(covariates %in% colnames(df))) {
      stop("All covariates must exist in df")
    }

    covariates <- tolower(covariates)
    xnames = paste(covariates, collapse = "+")
    formula = as.formula(paste("survival::Surv(time, dropout) ~", xnames))
  } else {
    formula = survival::Surv(time, dropout) ~ 1
  }


  setDT(df)
  setnames(df, tolower(names(df)))

  if (by_treatment) {
    ngroups = df[, uniqueN(get("treatment"))]

    if (!("treatment_description" %in% names(df))) {
      df[, `:=`(treatment_description = paste("Treatment", get("treatment")))]
    }
  } else {
    ngroups = 1
    df[, `:=`(treatment = 1)]
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }


  # fit by treatment group
  dropout_fit <- list()

  for (i in 1:ngroups) {
    df1 <- df[get("treatment") == i]

    n0 = nrow(df1)
    c0 = df1[, sum(get("dropout"))]
    ex0 = df1[, sum(get("time"))]

    x = model.matrix(formula, df1)
    q = ncol(x) - 1

    kmfit <- survival::survfit(survival::Surv(time, dropout) ~ 1, data = df1)
    kmdf <- data.table(time = kmfit$time, surv = kmfit$surv)
    df0 <- data.table(time = 0, surv = 1)
    kmdf <- rbindlist(list(df0, kmdf), use.names = TRUE)

    if (tolower(dropout_model) == "exponential") {
      erify::check_positive(c0 - q, supplement = paste(
        "The number of dropouts must be >=", q + 1,
        "to fit an exponential model."))

      # lambda(t) = lambda
      # S(t) = exp(-lambda*t)

      reg <- survival::survreg(formula, data = df1, dist = "exponential")

      # use the more common parameterization for the exponential distribution
      fit3 <- list(model = 'Exponential',
                   theta = -as.numeric(reg$coefficients),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 2*(q+1),
                   bic = -2*reg$loglik[1] + (q+1)*log(n0))

      # fitted survival curve
      rate = exp(as.numeric(x %*% fit3$theta))

      dffit3 <- data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(pexp(t, rate, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "weibull") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a Weibull model."))

      # lambda(t) = kappa/lambda*(t/lambda)^(kappa-1)
      # S(t) = exp(-(t/lambda)^kappa)

      reg <- survival::survreg(formula, data = df1, dist = "weibull")

      # weibull$shape = 1/reg$scale, weibull$scale = exp(reg$coefficients)
      # we define theta = c(log(weibull$scale), -log(weibull$shape))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit3 <- list(model = "Weibull",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 2*(q+2),
                   bic = -2*reg$loglik[1] + (q+2)*log(n0))

      # fitted survival curve
      shape = exp(-fit3$theta[q+2])
      scale = exp(as.numeric(x %*% fit3$theta[1:(q+1)]))

      dffit3 <- data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(pweibull(t, shape, scale, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "log-logistic") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a log-logistic model."))

      # S(t) = 1/(1 + (t/lambda)^kappa)

      reg <- survival::survreg(formula, data = df1, dist = "loglogistic")

      # llogis$shape = 1/reg$scale, llogis$scale = exp(reg$coefficients)
      # we define theta = (log(llogis$scale), -log(llogis$shape))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit3 <- list(model = "Log-logistic",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 2*(q+2),
                   bic = -2*reg$loglik[1] + (q+2)*log(n0))

      # fitted survival curve
      location = as.numeric(x %*% fit3$theta[1:(q+1)])
      scale = exp(fit3$theta[q+2])

      dffit3 <- data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(plogis(log(t), location, scale, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "log-normal") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a log-normal model."))

      # S(t) = 1 - Phi((log(t) - meanlog)/sdlog)

      reg <- survival::survreg(formula, data = df1, dist = "lognormal")

      # we use parameterization theta = (meanlog, log(sdlog))
      # reg$var is for theta = c(reg$coefficients, log(reg$scale))
      fit3 <- list(model = "Log-normal",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2*reg$loglik[1] + 2*(q+2),
                   bic = -2*reg$loglik[1] + (q+2)*log(n0))

      # fitted survival curve
      meanlog = as.numeric(x %*% fit3$theta[1:(q+1)])
      sdlog = exp(fit3$theta[q+2])

      dffit3 <- data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(plnorm(t, meanlog, sdlog, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "piecewise exponential") {
      # lambda(t) = lambda[j] for ucut[j] <= t < ucut[j+1], j = 1,...,J
      # where ucut[1]=0< ucut[2] < ... < ucut[J] < ucut[J+1]=Inf are the knots
      J = length(piecewiseDropoutTime)

      erify::check_positive(c0 - J - q + 1, supplement = paste(
        "The number of dropouts must be >=", J + q,
        "to fit a piecewise exponential model."))

      # maximum likelihood estimates and covariance matrix
      fit3 <- pwexpreg(df1$time, df1$dropout, J, piecewiseDropoutTime, q, x)

      # fitted survival curve
      time = seq(0, max(df1$time))

      surv = purrr::map(1:n0, function(l)
        ppwexp(time, fit3$theta, J, fit3$piecewiseDropoutTime,
               q, x[l,], lower.tail = FALSE))
      surv = apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)

      dffit3 <- data.table(time, surv)
    } else if (tolower(dropout_model) == "model averaging") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a model averaging model."))

      reg1 <- survival::survreg(formula, data = df1, dist = "weibull")
      reg2 <- survival::survreg(formula, data = df1, dist = "lognormal")
      aic1 <- -2*reg1$loglik[1] + 2*(q+2)
      aic2 <- -2*reg2$loglik[1] + 2*(q+2)
      bic1 <- -2*reg1$loglik[1] + (q+2)*log(n0)
      bic2 <- -2*reg2$loglik[1] + (q+2)*log(n0)

      w1 = 1/(1 + exp(-0.5*(bic2 - bic1)))

      # model parameters from weibull and log-normal
      theta = c(as.numeric(reg1$coefficients), log(reg1$scale),
                as.numeric(reg2$coefficients), log(reg2$scale))

      # variance-covariance matrix, noting that the covariances
      # between the two sets of parameters are zero as they are estimated
      # from different likelihood functions
      vtheta = as.matrix(Matrix::bdiag(reg1$var, reg2$var))

      # model fit, assuming fixed weight w1
      fit3 <- list(model = "Model averaging",
                   theta = theta,
                   vtheta = vtheta,
                   aic = w1*aic1 + (1-w1)*aic2,
                   bic = w1*bic1 + (1-w1)*bic2,
                   w1 = w1)

      # fitted survival curve
      time = seq(0, max(df1$time))

      surv = purrr::map(1:n0, function(l)
        pmodavg(time, fit3$theta, w1, q, x[l,], lower.tail = FALSE))
      surv = apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)

      dffit3 <- data.table(time, surv)
    } else if (tolower(dropout_model) == "spline") {
      erify::check_positive(c0 - k_dropout - q - 1, supplement = paste(
        "The number of dropouts must be >=", k_dropout + q + 2,
        "to fit a spline model."))

      # g(S(t)) = gamma_0 +gamma_1*x +gamma_2*v_1(x) +... +gamma_{k+1}*v_k(x)

      spl <- flexsurv::flexsurvspline(
        formula, data = df1, k = k_dropout, scale = scale_dropout,
        method = "Nelder-Mead")

      fit3 <- list(model = "Spline",
                   theta = as.numeric(spl$coefficients),
                   vtheta = spl$cov,
                   aic = -2*spl$loglik + 2*(k_dropout+q+2),
                   bic = -2*spl$loglik + (k_dropout+q+2)*log(n0),
                   knots = spl$knots,
                   scale = spl$scale)

      # fitted survival curve
      time = seq(0, max(df1$time))

      if (q > 0) {
        xbeta = as.numeric(as.matrix(x[,-1]) %*%
                             fit3$theta[(k_dropout+3):(k_dropout+q+2)])

        surv = purrr::map(1:n0, function(l)
          flexsurv::psurvspline(
            time, gamma = fit3$theta[1:(k_dropout+2)], knots = fit3$knots,
            scale = fit3$scale, offset = xbeta[l], lower.tail = FALSE))
        surv = apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
      } else {
        surv = flexsurv::psurvspline(
          time, gamma = fit3$theta, knots = fit3$knots, scale = fit3$scale,
          lower.tail = FALSE)
      }

      dffit3 <- data.table(time, surv)
    }


    # plot the survival curve
    if (tolower(fit3$model) == "piecewise exponential") {
      modeltext = paste0(paste0(fit3$model, "("),
                         paste(piecewiseDropoutTime, collapse = ","), ")")
    } else if (tolower(fit3$model) == "spline") {
      modeltext = paste0(fit3$model, "(k = ", k_dropout, ", ", "scale = '",
                         scale_dropout, "')")
    } else {
      modeltext = fit3$model
    }

    if (tolower(dropout_model) == "model averaging") {
      aictext = paste("Weighted AIC:",
                      formatC(fit3$aic, format = "f", digits = 2))
      bictext = paste("Weighted BIC:",
                      formatC(fit3$bic, format = "f", digits = 2))
    } else {
      aictext = paste("AIC:", formatC(fit3$aic, format = "f", digits = 2))
      bictext = paste("BIC:", formatC(fit3$bic, format = "f", digits = 2))
    }

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
          x = c(0.7, 0.7, 0.7), y = c(0.95, 0.80, 0.65), xref = "paper",
          yref = "paper", text = paste('<i>', c(modeltext, aictext,
                                                bictext), '</i>'),
          xanchor = "left", font = list(size = 14, color = "red"),
          showarrow = FALSE)) %>%
      plotly::hide_legend()

    if (by_treatment) {
      fittedDropout <- fittedDropout %>%
        plotly::layout(annotations = list(
          x = 0.5, y = 1,
          text = paste0("<b>", df1$treatment_description[1], "</b>"),
          xanchor = "center", yanchor = "middle", showarrow = FALSE,
          xref = 'paper', yref = 'paper'))

      fit3$treatment = df1$treatment[1]
      fit3$treatment_description = df1$treatment_description[1]
    }

    dropout_fit[[i]] = list(fit = fit3, fit_plot = fittedDropout,
                            kmdf = kmdf, dffit = dffit3,
                            text = c(modeltext, aictext, bictext))
  }

  # ensure that the sub plots share the same x axis range
  if (by_treatment) {
    x_range = range(df$time)
    for (i in 1:ngroups) {
      dropout_fit[[i]]$fit_plot <- dropout_fit[[i]]$fit_plot %>%
        plotly::layout(xaxis = list(range = x_range))
    }
  } else {
    dropout_fit = list(fit = fit3, fit_plot = fittedDropout,
                       kmdf = kmdf, dffit = dffit3,
                       text = c(modeltext, aictext, bictext))
  }

  if (showplot) {
    if (by_treatment) {
      for (i in 1:ngroups) {
        print(dropout_fit[[i]]$fit_plot)
      }
    } else {
      print(dropout_fit$fit_plot)
    }
  }

  dropout_fit
}
