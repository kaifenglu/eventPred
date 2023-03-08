#' @title Fit enrollment model
#' @description Fit a specified enrollment model to the enrollment data.
#'
#' @param df Input data set containing \code{RANDDT}.
#' @param enroll_model Enrollment model: Poisson (homogeneous), time-decay,
#'   or B-spline.
#' @param nknots Number of inner knots for B-spline. Defaults to 1.
#'
#' @details For time-decay model, the mean function is
#' \code{mu(t) = (mu/delta) (t - (1/delta)(1 - exp(-delta*t)))}
#' and the rate function is
#' \code{lambda(t) = (mu/delta) (1 - exp(-delta*t))}
#' For B-spline model, the daily enrollment rate is approximated using
#' the exponential of a B-spline function:
#' \code{lambda(t) = exp(B(t)*theta)}
#' where B(t) are the B-spline basis functions.
#'
#' @return A list of results from model fit including: \code{model},
#'   \code{theta}, \code{vtheta}, \code{bic}, and the design marix
#'   \code{x} for B-spline enrollment model.
#'
#' @examples
#'
#' observed <- summarizeObserved(df = observedData,
#'                               to_predict = "enrollment and event")
#'
#' fitEnr <- fitEnrollment(df = observed$adsl, enroll_model = "b-spline",
#'                         nknots = 1)
#'
#' @export
#'
fitEnrollment <- function(df, enroll_model, nknots = 1) {
  erify::check_class(df, "data.frame")

  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay", "b-spline"))

  erify::check_n(nknots)

  trialsdt = min(df$RANDDT)
  cutoffdt = df$CUTOFFDT[1]
  n0 = nrow(df)
  t0 = as.numeric(cutoffdt - trialsdt + 1)

  df1 <- df %>%
    mutate(time = as.numeric(.data$RANDDT - trialsdt + 1))


  if (tolower(enroll_model) == "poisson") {
    # a(t) = lambda
    # mu(t) = lambda*t
    fit1 <- list(model = 'Poisson',
                 theta = log(n0/t0),
                 vtheta = 1/n0,
                 bic = -2*(-n0 + n0*log(n0/t0)) + log(n0))

    dffit1 <- tibble(
      time = seq(1, t0),
      n = exp(fit1$theta)*.data$time)

    p1 <- ggplot() +
      geom_step(data = df1, aes(x = .data$time, y = .data$n)) +
      geom_line(data = dffit1, aes(x = .data$time, y = .data$n),
                color="blue") +
      labs(x = "Days since trial start",
           y = "Subjects",
           title = "Fitted enrollment curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit1$model, x=0.05, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit1$bic, 2)), x=0.05, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEnroll <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEnroll)
  } else if (tolower(enroll_model) == "time-decay") {
    # a(t) = mu/delta*(1 - exp(-delta*t))
    # mu(t) = mu/delta*(t - 1/delta*(1 - exp(-delta*t)))
    llik_td <- function(theta, t, df) {
      mu = exp(theta[1])
      delta = exp(theta[2])
      a1 = -mu/delta*(t - 1/delta*(1 - exp(-delta*t)))
      a2 = sum(log(mu/delta) + log(1 - exp(-delta*df$time)))
      a1 + a2
    }

    theta <- c(log(n0/t0), 0)
    opt1 <- optim(theta, llik_td, gr = NULL, t = t0, df = df1,
                  control = c(fnscale = -1))  # maximization
    fit1 <- list(model = "Time-decay",
                 theta = opt1$par,
                 vtheta = solve(-optimHess(opt1$par, llik_td, gr = NULL,
                                           t = t0, df = df1)),
                 bic = -2*llik_td(opt1$par, t = t0, df = df1) + 2*log(n0))

    # mean function of the NHPP
    fmu_td <- function(t, theta) {
      mu = exp(theta[1])
      delta = exp(theta[2])
      mu/delta*(t - 1/delta*(1 - exp(-delta*t)))
    }

    dffit1 <- tibble(
      time = seq(1, t0),
      n = fmu_td(.data$time, fit1$theta))

    p1 <- ggplot() +
      geom_step(data = df1, aes(x = .data$time, y = .data$n)) +
      geom_line(data = dffit1, aes(x = .data$time, y = .data$n),
                color="blue") +
      labs(x = "Days since trial start",
           y = "Subjects",
           title = "Fitted enrollment curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit1$model, x=0.05, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit1$bic, 2)), x=0.05, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEnroll <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEnroll)
  } else if (tolower(enroll_model) == "b-spline") {
    # a(t) = exp(theta' bs(t))
    # mu(t) = sum(a(u), {u,1,t})

    # number of inner knots
    K = nknots

    days = seq(1, t0)
    n = sapply(days, function(i) sum(df1$time == i))

    # design matrix for cubic B-spline
    x = splines::bs(days, df=K+4, intercept=1)

    # log-likelihood with B-spline fit for log(lambda(t))
    llik_bs <- function(theta, n, x) {
      lambda = exp(as.vector(x %*% theta)) # daily enrollment rate
      -sum(lambda) + sum(n*log(lambda))
    }

    # maximum likelihood estimation with initial value from OLS
    theta = as.vector(solve(t(x) %*% x, t(x) %*% log(pmax(n, 0.1))))
    opt1 <- optim(theta, llik_bs, gr = NULL, n = n, x = x,
                  control = c(fnscale = -1))
    fit1 <- list(model = "B-spline",
                 theta = opt1$par,
                 vtheta = solve(-optimHess(opt1$par, llik_bs, gr = NULL,
                                           n = n, x = x)),
                 bic = -2*opt1$value + (K+4)*log(n0),
                 x = x)


    # mean function of the NHPP, assuming t <= t0
    fmu_bs <- function(t, theta, x) {
      lambda = exp(as.vector(x %*% theta))
      lambdasum = cumsum(lambda)
      lambdasum[t]
    }

    dffit1 <- tibble(
      time = seq(1, t0),
      n = fmu_bs(.data$time, fit1$theta, x))

    p1 <- ggplot() +
      geom_step(data = df1, aes(x = .data$time, y = .data$n)) +
      geom_line(data = dffit1, aes(x = .data$time, y = .data$n),
                color="blue") +
      labs(x = "Days since trial start",
           y = "Subjects",
           title = "Fitted enrollment curve") +
      theme_bw()

    grob1 <- grid::grobTree(grid::textGrob(
      fit1$model, x=0.05, y=0.95, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    grob2 <- grid::grobTree(grid::textGrob(
      paste("BIC:", round(fit1$bic, 2)), x=0.05, y=0.88, hjust=0,
      gp=grid::gpar(col="red", fontsize=11, fontface="italic")))

    fittedEnroll <- p1 + annotation_custom(grob1) + annotation_custom(grob2)
    print(fittedEnroll)
  }

  fit1
}
