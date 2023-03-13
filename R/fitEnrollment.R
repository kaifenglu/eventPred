#' @title Fit enrollment model
#' @description Fits a specified enrollment model to the enrollment data.
#'
#' @param df The subject-level enrollment data, including
#'   \code{randdt} and \code{cutoffdt}.
#' @param enroll_model The enrollment model which can be specified as
#'   "Poisson", "Time-decay", "B-spline", or
#'   "Piecewise Poisson". By default, it is set to "B-spline".
#' @param nknots The number of inner knots for the B-spline enrollment
#'   model. By default, it is set to 1.
#' @param accrualTime The accrual time intervals for the piecewise Poisson
#'   model. Must start with 0, e.g., c(0, 3) breaks the time axis into
#'   2 accrual intervals: [0, 3) and [3, Inf). By default, it is set to 0.
#'
#' @details
#' For the time-decay model, the mean function is
#' \code{mu(t) = (mu/delta) (t - (1/delta)(1 - exp(-delta*t)))}
#' and the rate function is
#' \code{lambda(t) = (mu/delta) (1 - exp(-delta*t))}.
#' For the B-spline model, the daily enrollment rate is approximated as
#' \code{lambda(t) = exp(B(t)*theta)},
#' where \code{B(t)} represents the B-spline basis functions.
#'
#' @return
#' A list of results from the model fit including key information
#' such as the enrollment model, \code{model}, the estimated model
#' parameters, \code{theta}, the covariance matrix, \code{vtheta}, and
#' the Bayesian information criterion, \code{bic}, as well as
#' the design matrix \code{x} for the B-spline enrollment model, and
#' \code{accrualTime} for the piecewise Poisson enrollment model.
#'
#' @examples
#'
#' enroll_fit <- fitEnrollment(df = observedData, enroll_model = "b-spline",
#'                             nknots = 1)
#'
#' @export
#'
fitEnrollment <- function(df, enroll_model = "b-spline", nknots = 1,
                          accrualTime = 0) {
  erify::check_class(df, "data.frame")

  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay",
                         "b-spline", "piecewise poisson"))

  erify::check_n(nknots)

  df <- dplyr::as_tibble(df)
  names(df) <- tolower(names(df))
  trialsdt = min(df$randdt)
  cutoffdt = df$cutoffdt[1]
  n0 = nrow(df)
  t0 = as.numeric(cutoffdt - trialsdt + 1)

  df1 <- df %>%
    dplyr::arrange(.data$randdt) %>%
    dplyr::mutate(time = as.numeric(.data$randdt - trialsdt + 1),
                  n = dplyr::row_number())

  # fit enrollment model
  if (tolower(enroll_model) == "poisson") {
    # lambda(t) = lambda
    # mu(t) = lambda*t
    fit1 <- list(model = 'Poisson',
                 theta = log(n0/t0),
                 vtheta = 1/n0,
                 bic = -2*(-n0 + n0*log(n0/t0)) + log(n0))

    dffit1 <- dplyr::tibble(
      time = seq(1, t0),
      n = exp(fit1$theta)*.data$time)
  } else if (tolower(enroll_model) == "time-decay") {
    # lambda(t) = mu/delta*(1 - exp(-delta*t))
    # mu(t) = mu/delta*(t - 1/delta*(1 - exp(-delta*t)))
    # mean function of the NHPP
    fmu_td <- function(t, theta) {
      mu = exp(theta[1])
      delta = exp(theta[2])
      mu/delta*(t - 1/delta*(1 - exp(-delta*t)))
    }

    llik_td <- function(theta, t, df) {
      mu = exp(theta[1])
      delta = exp(theta[2])
      a1 = -mu/delta*(t - 1/delta*(1 - exp(-delta*t)))
      a2 = sum(log(mu/delta) + log(1 - exp(-delta*df$time)))
      a1 + a2
    }

    # slope in the last 1/4 enrollment time interval
    beta = (n0 - df1$n[df1$time >= 3/4*t0][1])/(1/4*t0)
    mu0 = 2*n0/t0^2
    delta0 = mu0/beta
    theta <- c(log(mu0), log(delta0))
    opt1 <- optim(theta, llik_td, gr = NULL, t = t0, df = df1,
                  control = c(fnscale = -1))  # maximization
    fit1 <- list(model = "Time-decay",
                 theta = opt1$par,
                 vtheta = solve(-optimHess(opt1$par, llik_td, gr = NULL,
                                           t = t0, df = df1)),
                 bic = -2*llik_td(opt1$par, t = t0, df = df1) + 2*log(n0))

    dffit1 <- dplyr::tibble(
      time = seq(1, t0),
      n = fmu_td(.data$time, fit1$theta))
  } else if (tolower(enroll_model) == "b-spline") {
    # lambda(t) = exp(theta' bs(t))
    # mu(t) = sum(lambda(u), {u,1,t})

    # number of inner knots
    K = nknots

    days = seq(1, t0)
    n = as.numeric(table(factor(df1$time, levels = days)))

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

    dffit1 <- dplyr::tibble(
      time = seq(1, t0),
      n = fmu_bs(.data$time, fit1$theta, x))
  } else if (tolower(enroll_model) == "piecewise poisson") {
    # truncate the time intervals by data cut
    u = accrualTime[accrualTime < t0]
    u2 = c(u, t0)

    # number of enrolled subjects in each interval
    factors <- cut(df1$time, breaks = u2)
    n = table(factors)

    # length of each interval
    t = diff(u2)

    # constant enrollment rate in each interval
    fit1 <- list(model = 'Piecewise Poisson',
                 theta = log(n/t),
                 vtheta = diag(1/n),
                 accrualTime = u,
                 bic = -2*sum(-n + n*log(n/t)) + length(u)*log(n0))

    lambda = n/t
    psum = c(0, cumsum(n))  # cumulative enrollment by end of interval
    time = seq(1, t0) # find the time interval for each day
    j = findInterval(time, u, left.open = TRUE)
    m = psum[j] + lambda[j]*(time - u[j]) # cumulative enrollment by day

    dffit1 <- dplyr::tibble(
      time = time,
      n = m)
  }

  # plot the enrollment curve
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_step(data = df1, ggplot2::aes(x = .data$time,
                                                y = .data$n)) +
    ggplot2::geom_line(data = dffit1, ggplot2::aes(x = .data$time,
                                                   y = .data$n),
                       color = "blue") +
    ggplot2::labs(x = "Days since trial start",
                  y = "Subjects",
                  title = "Fitted enrollment curve") +
    ggplot2::theme_bw()

  grob1 <- grid::grobTree(grid::textGrob(
    label = fit1$model, x=0.05, y=0.95, hjust=0,
    gp = grid::gpar(col="red", fontsize=11, fontface="italic")))

  grob2 <- grid::grobTree(grid::textGrob(
    label = paste("BIC:", round(fit1$bic,2)), x=0.05, y=0.88, hjust=0,
    gp = grid::gpar(col="red", fontsize=11, fontface="italic")))

  fittedEnroll <- p1 + ggplot2::annotation_custom(grob1) +
    ggplot2::annotation_custom(grob2)
  print(fittedEnroll)


  fit1
}
