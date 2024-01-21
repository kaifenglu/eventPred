#' @title Fit enrollment model
#' @description Fits a specified enrollment model to the enrollment data.
#'
#' @param df The subject-level enrollment data, including \code{trialsdt},
#'   \code{randdt} and \code{cutoffdt}.
#' @param enroll_model The enrollment model which can be specified as
#'   "Poisson", "Time-decay", "B-spline", or
#'   "Piecewise Poisson". By default, it is set to "B-spline".
#' @param nknots The number of inner knots for the B-spline enrollment
#'   model. By default, it is set to 0.
#' @param accrualTime The accrual time intervals for the piecewise Poisson
#'   model. Must start with 0, e.g., c(0, 30) breaks the time axis into
#'   2 accrual intervals: [0, 30) and [30, Inf). By default, it is set to 0.
#' @param showplot A Boolean variable to control whether or not to
#'   show the fitted enrollment curve. By default, it is set to \code{TRUE}.
#'
#' @details
#' For the time-decay model, the mean function is
#' \code{mu(t) = mu/delta*(t - 1/delta*(1 - exp(-delta*t)))}
#' and the rate function is
#' \code{lambda(t) = mu/delta*(1 - exp(-delta*t))}.
#' For the B-spline model, the daily enrollment rate is approximated as
#' \code{lambda(t) = exp(B(t)*theta)},
#' where \code{B(t)} represents the B-spline basis functions.
#'
#' @return
#' A list of results from the model fit including key information
#' such as the enrollment model, \code{model}, the estimated model
#' parameters, \code{theta}, the covariance matrix, \code{vtheta},
#' the Akaike Information Criterion, \code{aic}, and
#' the Bayesian Information Criterion, \code{bic}, as well as
#' the design matrix \code{x} for the B-spline enrollment model, and
#' \code{accrualTime} for the piecewise Poisson enrollment model.
#'
#' The fitted enrollment curve is also returned.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Xiaoxi Zhang and Qi Long. Stochastic modeling and prediction for
#' accrual in clinical trials. Stat in Med. 2010; 29:649-658.
#'
#' @examples
#'
#' enroll_fit <- fitEnrollment(df = interimData1, enroll_model = "b-spline",
#'                             nknots = 1)
#'
#' @export
#'
fitEnrollment <- function(df, enroll_model = "b-spline", nknots = 0,
                          accrualTime = 0, showplot = TRUE) {
  erify::check_class(df, "data.frame")

  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay",
                         "b-spline", "piecewise poisson"))

  erify::check_n(nknots, zero = TRUE)

  if (accrualTime[1] != 0) {
    stop("accrualTime must start with 0");
  }
  if (length(accrualTime) > 1 && any(diff(accrualTime) <= 0)) {
    stop("accrualTime should be increasing")
  }

  erify::check_bool(showplot)


  df <- dplyr::as_tibble(df)
  names(df) <- tolower(names(df))
  df$trialsdt <- as.Date(df$trialsdt)
  df$randdt <- as.Date(df$randdt)
  df$cutoffdt <- as.Date(df$cutoffdt)

  trialsdt = df$trialsdt[1]
  cutoffdt = df$cutoffdt[1]
  n0 = nrow(df)
  # up to the last randomization date to account for enrollment completion
  t0 = as.numeric(max(df$randdt) - trialsdt + 1)

  erify::check_positive(n0, supplement = paste(
    "The number of subjects must be positive to fit an enrollment model."))

  if (any(df$randdt < trialsdt)) {
    stop("randdt must be greater than or equal to trialsdt.")
  }

  if (any(df$randdt > cutoffdt)) {
    stop("randdt must be less than or equal to cutoffdt.")
  }

  df1 <- df %>%
    dplyr::arrange(.data$randdt) %>%
    dplyr::mutate(t = as.numeric(.data$randdt - trialsdt + 1),
                  n = dplyr::row_number())

  df1u <- df1 %>%
    dplyr::group_by(.data$randdt) %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$t, .data$n)

  # add day 1
  df0 <- dplyr::tibble(t = 1, n = 0)
  df1u <- df0 %>%
    dplyr::bind_rows(df1u)

  # fit enrollment model
  if (tolower(enroll_model) == "poisson") {
    # lambda(t) = lambda
    # mu(t) = lambda*t
    fit1 <- list(model = 'Poisson',
                 theta = log(n0/t0),
                 vtheta = 1/n0,
                 aic = -2*(-n0 + n0*log(n0/t0)) + 2,
                 bic = -2*(-n0 + n0*log(n0/t0)) + log(n0))

    dffit1 <- dplyr::tibble(
      t = seq(1, t0),
      n = exp(fit1$theta)*.data$t)
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
      a2 = sum(log(mu/delta) + log(1 - exp(-delta*df$t)))
      a1 + a2
    }

    # slope in the last 1/4 "active" enrollment time interval
    beta = (n0 - df1u$n[df1u$t >= 3/4*t0][1])/(1/4*t0)
    mu0 = 2*n0/t0^2  # Taylor expansion of mu(t) to t^2
    delta0 = mu0/beta  # beta is the asymptotic slope
    theta <- c(log(mu0), log(delta0))
    opt1 <- optim(theta, llik_td, gr = NULL, t = t0, df = df1,
                  control = c(fnscale = -1), hessian = TRUE)  # maximization
    fit1 <- list(model = "Time-decay",
                 theta = opt1$par,
                 vtheta = solve(-opt1$hessian),
                 aic = -2*opt1$value + 4,
                 bic = -2*opt1$value + 2*log(n0))

    dffit1 <- dplyr::tibble(
      t = seq(1, t0),
      n = fmu_td(.data$t, fit1$theta))
  } else if (tolower(enroll_model) == "b-spline") {
    # lambda(t) = exp(theta' bs(t))
    # mu(t) = sum(lambda(u), {u,1,t})

    # number of inner knots
    K = nknots

    days = seq(1, t0)
    n = as.numeric(table(factor(df1$t, levels = days)))

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
                  control = c(fnscale = -1), hessian = TRUE)
    fit1 <- list(model = "B-spline",
                 theta = opt1$par,
                 vtheta = solve(-opt1$hessian),
                 aic = -2*opt1$value + 2*(K+4),
                 bic = -2*opt1$value + (K+4)*log(n0),
                 x = x)

    # mean function of the NHPP, assuming t <= t0
    fmu_bs <- function(t, theta, x) {
      lambda = exp(as.vector(x %*% theta))
      lambdasum = cumsum(lambda)
      lambdasum[t]
    }

    dffit1 <- dplyr::tibble(
      t = seq(1, t0),
      n = fmu_bs(.data$t, fit1$theta, x))
  } else if (tolower(enroll_model) == "piecewise poisson") {
    # truncate the time intervals by data cut
    u = accrualTime[accrualTime < t0]
    u2 = c(u, t0)

    # number of enrolled subjects in each interval
    factors <- cut(df1$t, breaks = u2)
    n = as.numeric(table(factors))

    # length of each interval
    t = diff(u2)

    # covariance matrix
    if (length(u) > 1) {
      vtheta = diag(1/n)
    } else {
      vtheta = 1/n*diag(1)
    }

    # constant enrollment rate in each interval
    fit1 <- list(model = 'Piecewise Poisson',
                 theta = log(n/t),
                 vtheta = vtheta,
                 aic = -2*sum(-n + n*log(n/t)) + 2*length(u),
                 bic = -2*sum(-n + n*log(n/t)) + length(u)*log(n0),
                 accrualTime = u)

    lambda = n/t
    psum = c(0, cumsum(n))  # cumulative enrollment by end of interval
    time = seq(1, t0) # find the time interval for each day
    j = findInterval(time, u)
    m = psum[j] + lambda[j]*(time - u[j]) # cumulative enrollment by day

    dffit1 <- dplyr::tibble(
      t = time,
      n = m)
  }


  # plot the survival curve
  if (tolower(fit1$model) == "piecewise poisson") {
    modeltext = paste0(paste0(fit1$model, "("),
                       paste(accrualTime, collapse = " "), ")")
  } else if (tolower(fit1$model) == "b-spline") {
    modeltext = paste0(fit1$model, "(nknots = ", nknots, ")")
  } else {
    modeltext = fit1$model
  }

  aictext = paste("AIC:", formatC(fit1$aic, format = "f", digits = 2))
  bictext = paste("BIC:", formatC(fit1$bic, format = "f", digits = 2))

  # plot the enrollment curve
  fittedEnroll <- plotly::plot_ly() %>%
    plotly::add_lines(
      data=df1u, x=~t, y=~n, name="observed", line=list(shape="hv")) %>%
    plotly::add_lines(data=dffit1, x=~t, y=~n, name="fitted") %>%
    plotly::layout(
      xaxis = list(title = "Days since trial start", zeroline = FALSE),
      yaxis = list(title = "Subjects", zeroline = FALSE),
      title = list(text = "Fitted enrollment curve"),
      annotations = list(
        x = c(0.05, 0.05, 0.05), y = c(0.95, 0.80, 0.65),
        xref = "paper", yref = "paper",
        text = paste('<i>', c(modeltext, aictext, bictext), '</i>'),
        xanchor = "left",
        font = list(size = 14, color = "red"), showarrow = FALSE)) %>%
    plotly::hide_legend()

  if (showplot) print(fittedEnroll)

  list(enroll_fit = fit1, enroll_fit_plot = fittedEnroll)
}
