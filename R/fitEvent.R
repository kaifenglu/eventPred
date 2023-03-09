#' @title Fit time-to-event model
#' @description Fits a specified time-to-event model to the event data.
#'
#' @param df The input data set which includes information on
#'   \code{time} and \code{event}.
#' @param event_model The event model which specifies the type of event
#'   model to be used in the analysis and can be set to one of the
#'   following options: exponential, Weibull, log-normal,
#'   piecewise exponential, or model averaging, which uses the
#'   \code{exp(-bic)} weighting and combines Weibull and
#'   log-normal models. If not provided, defaults to model
#'   averaging.
#' @param npieces The number of pieces for the piecewise exponential
#'   event model. If not provided, defaults to 3.
#'
#' @return
#' A list of results from the model fit including key information
#' such as the event model, \code{model}, the estimated model parameters,
#' \code{theta}, the covariance matrix, \code{vtheta}, as well as the
#' Bayesian Information Criterion, \code{bic}.
#'
#' If the piecewise exponential model was used, additional
#' variables will be included in the list of results, such as the
#' number of pieces specified, \code{npieces}, and the location
#' of knots used in the model, \code{knots}.
#'
#' If the model averaging option is chosen, the weight assigned
#' to the Weibull component is indicated by the \code{w1} variable.
#'
#' @examples
#'
#' event_fit <- fitEvent(df = observedData,
#'                       event_model = "piecewise exponential", npieces = 3)
#'
#' @export
#'
fitEvent <- function(df, event_model = "model averaging", npieces = 3) {
  erify::check_class(df, "data.frame")
  erify::check_content(tolower(event_model),
                       c("exponential", "weibull", "log-normal",
                         "piecewise exponential", "model averaging"))
  erify::check_n(npieces)

  names(df) <- tolower(names(df))
  n0 = nrow(df)
  d0 = sum(df$event)
  ex0 = sum(df$time)

  kmfit <- survival::survfit(survival::Surv(time, event) ~ 1, data = df)
  kmdf <- tibble(time = kmfit$time, surv = kmfit$surv)
  kmdf <- tibble(time = 0, surv = 1) %>%
    bind_rows(kmdf)

  if (tolower(event_model) == "exponential") {
    # lambda(t) = lambda
    # S(t) = exp(-lambda*t)

    fit2 <- list(model = 'Exponential',
                 theta = log(d0/ex0),
                 vtheta = 1/d0,
                 bic = -2*(-d0 + d0*log(d0/ex0)) + log(n0))

    # fitted survival curve
    dffit2 <- tibble(
      time = seq(0, max(df$time)),
      surv = pexp(.data$time, rate = exp(fit2$theta), lower.tail = FALSE))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit2, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to event curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit2$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit2$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEvent <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEvent)
  } else if (tolower(event_model) == "weibull") {
    # lambda(t) = kappa/lambda*(t/lambda)^(kappa-1)
    # S(t) = exp(-(t/lambda)^kappa)

    reg <- survival::survreg(survival::Surv(time, event) ~ 1,
                             data = df, dist = "weibull")

    # Note: weibull$shape = 1/reg$scale, weibull$scale = exp(reg$coefficients)
    # we use parameterization theta = (log(weibull$shape), log(weibull$scale))
    # reg$var is for c(reg$coefficients, log(reg$scale)) = lmat %*% theta
    lmat <- matrix(c(0, -1, 1, 0), nrow=2, ncol=2, byrow=TRUE)
    fit2 <- list(model = "Weibull",
                 theta = c(log(1/reg$scale), as.numeric(reg$coefficients)),
                 vtheta = lmat %*% reg$var %*% t(lmat),
                 bic = -2*reg$loglik[1] + 2*log(n0))

    # fitted survival curve
    dffit2 <- tibble(
      time = seq(0, max(df$time)),
      surv = pweibull(.data$time, shape = exp(fit2$theta[1]),
                      scale = exp(fit2$theta[2]), lower.tail = FALSE))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit2, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to event curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit2$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit2$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEvent <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEvent)
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
    dffit2 <- tibble(
      time = seq(0, max(df$time)),
      surv = plnorm(.data$time, meanlog = fit2$theta[1],
                    sdlog = exp(fit2$theta[2]), lower.tail = FALSE))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit2, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to event curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit2$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit2$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEvent <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEvent)
  } else if (tolower(event_model) == "piecewise exponential") {
    # lambda(t) = lambda[j] for tau[j-1] < t <= tau[j], j = 1,...,J
    # where tau[0] = 0 < tau[1] < ... < tau[J-1] < tau[J] = Inf are the knots

    # number of inner knots
    J = npieces

    # number of events in each interval
    d = diff(c(0, round(d0/J*(1:J))))

    # cumulative number of events at the end of the first J-1 intervals
    dcut = cumsum(d[1:(J-1)])

    # ordered event times
    u = sort(df$time[df$event==1])

    # end points of the time intervals
    ucut = c(0, (u[dcut] + u[dcut+1])/2, max(df$time))

    # total exposure within each interval
    fex <- function(t, J, ucut) {
      ex = rep(NA, J)
      for (j in 1:J) {
        ex[j] = pmax(0, pmin(t, ucut[j+1]) - ucut[j])
      }
      ex
    }

    ex = rowSums(sapply(df$time, fex, J, ucut))

    # maximum likelihood estimates and covariance matrix
    fit2 <- list(model = "Piecewise exponential",
                 theta = log(d/ex),
                 vtheta = diag(1/d),
                 bic = -2*sum(-d + d*log(d/ex)) + J*log(n0),
                 npieces = J,
                 knots = ucut[2:J])

    # fitted survival curve
    dffit2 <- tibble(
      time = seq(0, max(df$time)),
      surv = exp(-colSums(exp(fit2$theta) * sapply(.data$time, fex, J, ucut))))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit2, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to event curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit2$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit2$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEvent <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEvent)
  } else if (tolower(event_model) == "model averaging") {
    reg1 <- survival::survreg(survival::Surv(time, event) ~ 1,
                              data = df, dist = "weibull")
    reg2 <- survival::survreg(survival::Surv(time, event) ~ 1,
                              data = df, dist = "lognormal")
    bic1 <- -2*reg1$loglik[1] + 2*log(n0)
    bic2 <- -2*reg2$loglik[1] + 2*log(n0)

    w1 = exp(-0.5*bic1)/(exp(-0.5*bic1) + exp(-0.5*bic2))

    # log-likelihood for model averaging of Weibull and log-normal
    llmodavg <- function(theta, w1, df) {
      shape = exp(theta[1])
      scale = exp(theta[2])
      meanlog = theta[3]
      sdlog = exp(theta[4])

      f1 = dweibull(df$time, shape, scale)
      f2 = dlnorm(df$time, meanlog, sdlog)
      f = w1*f1 + (1-w1)*f2

      s1 = pweibull(df$time, shape, scale, lower.tail = FALSE)
      s2 = plnorm(df$time, meanlog, sdlog, lower.tail = FALSE)
      s = w1*s1 + (1-w1)*s2

      l = df$event*log(f) + (1-df$event)*log(s)
      sum(l)
    }

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
                 bic = -2*llmodavg(theta, w1, df) + 4*log(n0),
                 w1 = w1)

    # distribution function for model averaging of Weibull and log-normal
    pmodavg <- function(t, theta, w1, lower.tail = TRUE) {
      shape = exp(theta[1])
      scale = exp(theta[2])
      meanlog = theta[3]
      sdlog = exp(theta[4])

      p1 = pweibull(pmax(0,t), shape, scale)
      p2 = plnorm(pmax(0,t), meanlog, sdlog)
      p = w1*p1 + (1-w1)*p2

      if (!lower.tail) p = 1 - p
      p
    }


    dffit2 <- tibble(
      time = seq(0, max(df$time)),
      surv = pmodavg(.data$time, theta, w1, lower.tail = FALSE))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit2, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to event curve") +
      theme_bw()


    grob1 <- grid::grobTree(grid::textGrob(
      fit2$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit2$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEvent <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEvent)
  }

  fit2

}
