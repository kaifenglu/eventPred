#' @title Fit time-to-dropout model
#' @description Fits a specified time-to-dropout model to the dropout data.
#'
#' @param df The subject-level dropout data, including \code{time} and
#'   \code{dropout}.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of three options: "exponential", "Weibull", or
#'   "log-normal". By default, it is set to "Weibull".
#'
#' @return A list of results from the model fit including key information
#'   such as the dropout model, \code{model}, the estimated model parameters,
#'   \code{theta}, the covariance matrix, \code{vtheta}, as well as the
#'   Bayesian Information Criterion, \code{bic}.
#'
#' @examples
#'
#' dropout_fit <- fitDropout(df = observedData, dropout_model = "exponential")
#'
#' @export
#'
fitDropout <- function(df, dropout_model = "weibull") {
  erify::check_class(df, "data.frame")
  erify::check_content(tolower(dropout_model),
                       c("exponential", "weibull", "log-normal"))

  names(df) <- tolower(names(df))
  n0 = nrow(df)
  d0 = sum(df$dropout)
  ex0 = sum(df$time)

  kmfit <- survival::survfit(survival::Surv(time, dropout) ~ 1, data = df)
  kmdf <- tibble(time = kmfit$time, surv = kmfit$surv)
  kmdf <- tibble(time = 0, surv = 1) %>%
    bind_rows(kmdf)

  if (tolower(dropout_model) == "exponential") {
    # lambda(t) = lambda
    # S(t) = exp(-lambda*t)

    fit3 <- list(model = 'Exponential',
                 theta = log(d0/ex0),
                 vtheta = 1/d0,
                 bic = -2*(-d0 + d0*log(d0/ex0)) + log(n0))

    # fitted survival curve
    dffit3 <- tibble(
      time = seq(0, max(df$time)),
      surv = pexp(.data$time, rate = exp(fit3$theta), lower.tail = FALSE))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit3, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to dropout curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit3$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit3$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedDropout <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedDropout)
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
    dffit3 <- tibble(
      time = seq(0, max(df$time)),
      surv = pweibull(.data$time, shape = exp(fit3$theta[1]),
                      scale = exp(fit3$theta[2]), lower.tail = FALSE))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit3, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to dropout curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit3$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit3$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedDropout <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedDropout)
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
    dffit3 <- tibble(
      time = seq(0, max(df$time)),
      surv = plnorm(.data$time, meanlog = fit3$theta[1],
                    sdlog = exp(fit3$theta[2]), lower.tail = FALSE))

    p1 <- ggplot() +
      geom_step(data = kmdf, aes(x = .data$time, y = .data$surv)) +
      geom_line(data = dffit3, aes(x = .data$time, y = .data$surv),
                color="blue") +
      labs(x = "Days since randomization",
           y = "Survival probability",
           title = "Fitted time to dropout curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit3$model, x=0.75, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit3$bic, 2)), x=0.75, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedDropout <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedDropout)
  }

  fit3

}
