#' @title Distribution function for piecewise exponential regression
#' @description Obtains the distribution function value for piecewise
#' exponential regression.
#'
#' @param t The vector of time points.
#' @param theta The parameter vector consisting of gamma for log
#'   piecewise hazards and beta for regression coefficients.
#' @param J The number of time intervals.
#' @param tcut A vector that specifies the time intervals
#'   for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2 event
#'   intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param q The number of elements in the vector of covariates.
#' @param x The vector of covariates (excluding the intercept).
#' @param lower.tail logical; if TRUE (default), probabilities are
#'   the distribution function, otherwise, the survival function.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return The probabilities p = P(T <= t | X = x).
#'
#' @keywords internal
#'
#' @export
#'
ppwexp <- function(t, theta, J, tcut, q = 0, x = 0, lower.tail = TRUE,
                   log.p = FALSE) {

  lambda = exp(theta[1:J]) # hazard rates in the intervals

  # partial sums of lambda*interval_width
  if (J > 1) {
    psum = c(0, cumsum(lambda[1:(J-1)] * diff(tcut[1:J])))
  } else {
    psum = 0
  }

  # find the interval containing time
  j1 = findInterval(t, tcut)
  cumhaz = psum[j1] + lambda[j1]*(t - tcut[j1])

  if (q > 0) {
    xbeta = sum(x * theta[(J+1):(J+q)])
    cumhaz = cumhaz*exp(xbeta)
  }

  p = 1 - exp(-cumhaz)

  if (!lower.tail) p = 1 - p
  if (log.p) p = log(p)
  p
}


#' @title Quantile function for piecewise exponential regression
#' @description Obtains the quantile function value for piecewise
#' exponential regression.
#'
#' @param p The vector of probabilities.
#' @param theta The parameter vector consisting of gamma for log
#'   piecewise hazards and beta for regression coefficients.
#' @param J The number of time intervals.
#' @param tcut A vector that specifies the endpoints of time intervals
#'   for the baseline piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2 event
#'   intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param q The number of elements in the vector of covariates.
#' @param x The vector of covariates (excluding the intercept).
#' @param lower.tail logical; if TRUE (default), probabilities are
#'   the distribution function, otherwise, the survival function.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return The quantiles t such that P(T <= t | X = x) = p.
#'
#' @keywords internal
#'
#' @export
#'
qpwexp <- function(p, theta, J, tcut, q = 0, x = 0, lower.tail = TRUE,
                   log.p = FALSE) {

  lambda = exp(theta[1:J]) # hazard rates in the intervals

  # partial sums of lambda*interval_width
  if (J > 1) {
    psum = c(0, cumsum(lambda[1:(J-1)] * diff(tcut[1:J])))
  } else {
    psum = 0
  }

  if (log.p) p = exp(p)
  if (!lower.tail) p = 1 - p

  cumhaz0 = -log(1-p)

  if (q > 0) {
    xbeta = sum(x * theta[(J+1):(J+q)])
    cumhaz0 = cumhaz0*exp(-xbeta)
  }

  # find the interval containing time
  j1 = findInterval(cumhaz0, psum)
  t = tcut[j1] + (cumhaz0 - psum[j1])/lambda[j1]
  t
}


#' @title Log likelihood for piecewise exponential regression
#' @description Obtains the log likelihood value for piecewise
#' exponential regression.
#'
#' @param theta The parameter vector consisting of gamma for log
#'   piecewise hazards and beta for regression coefficients.
#' @param time The survival time.
#' @param event The event indicator.
#' @param J The number of time intervals.
#' @param tcut A vector that specifies the endpoints of time intervals
#'   for the baseline piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2 event
#'   intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param q The number of columns of the covariates matrix.
#' @param x The covariates matrix (excluding the intercept).
#'
#' @return The log likelihood value for piecewise exponential regression.
#'
#' @keywords internal
#'
#' @export
#'
llik_pwexp <- function(theta, time, event, J, tcut, q, x) {
  if (length(tcut) == J) tcut = c(tcut, Inf)

  d = 0
  ex = 0
  for (j in 1:J) {
    d = d + theta[j] * (time > tcut[j]) * (time <= tcut[j+1])
    ex = ex + exp(theta[j]) * pmax(0, pmin(time, tcut[j+1]) - tcut[j])
  }

  if (q == 0) {
    sum(event*d - ex)
  } else {
    xbeta = as.numeric(x %*% theta[(J+1):(J+q)])
    sum(event * (d + xbeta) - ex * exp(xbeta))
  }
}


#' @title Piecewise exponential regression
#' @description Obtains the maximum likelihood estimates for piecewise
#' exponential regression.
#'
#' @param data The input data frame. It should contain \code{time},
#'   \code{event}, and variables specified in \code{covariates}.
#' @param tcut A vector that specifies the time intervals
#'   for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2 event
#'   intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param event The variable name indicating event, e.g. "event" or
#'   "dropout".
#' @param covariates The names of baseline covariates in
#'   the input data frame, e.g., c("age", "sex").
#'   Factor variables need to declared in the input data frame.
#'
#' @return The maximum likelihood estimates and the associated
#' covariance matrix, AIC and BIC.
#'
#' @export
#'
pwexpreg <- function(data, tcut = 0, event = "event", covariates = NULL) {
  if (!is.null(covariates)) {
    xnames = paste(covariates, collapse = "+")
    formula = as.formula(paste("~ -1 +", xnames))
    x = model.matrix(formula, data)
    q = ncol(x)
  } else {
    x = 0
    q = 0
  }

  u = tcut[tcut < max(data$time)]
  J = length(u)
  ucut = c(u, max(data$time))

  d = rep(0, J)  # number of events in each interval
  for (j in 1:J) {
    d[j] = sum((data$time > ucut[j]) * (data$time <= ucut[j+1]) * data$event)
  }

  if (any(d == 0)) {
    stop(paste("The number of", paste0(event, "s"),
               "must be >= 1 in each interval",
               "to fit a piecewise exponential model."))
  }


  theta0 = rep(0, J+q)

  opt1 <- optim(theta0, llik_pwexp, gr = NULL,
                data$time, data[, event], J, ucut, q, x,
                control = c(fnscale = -1), hessian = TRUE)

  n = nrow(data)

  fit <- list(model = "Piecewise exponential",
              theta = opt1$par,
              vtheta = solve(-opt1$hessian),
              aic = -2*opt1$value + 2*(J+q),
              bic = -2*opt1$value + (J+q)*log(n))

  if (tolower(event) == "event") {
    fit$piecewiseSurvivalTime = u
  } else {
    fit$piecewiseDropoutTime = u
  }

  fit
}


#' @title Distribution function for model averaging of Weibull and log-normal
#' @description Obtains the distribution function value for model-averaging
#' of Weibull and log-normal regression.
#'
#' @param t The vector of time points.
#' @param theta The parameter vector consisting of the accelerate failure
#'   time (AFT) regression coefficients and the logrithm of the AFT
#'   regression scale parameter for the Weibull and log-normal distributions.
#' @param w1 The weight for the Weibull component distribution.
#' @param q The number of elements in the vector of covariates.
#' @param x The vector of covariates (excluding the intercept).
#' @param lower.tail logical; if TRUE (default), probabilities are
#'   the distribution function, otherwise, the survival function.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return The probabilities p = P(T <= t | X = x).
#'
#' @keywords internal
#'
#' @export
#'
pmodavg <- function(t, theta, w1, q, x, lower.tail = TRUE, log.p = FALSE) {
  shape = exp(-theta[q+2])
  scale = exp(sum(c(1,x) * theta[1:(q+1)]))
  meanlog = sum(c(1,x) * theta[(q+3):(2*q+3)])
  sdlog = exp(theta[2*q+4])

  p1 = pweibull(t, shape, scale)
  p2 = plnorm(t, meanlog, sdlog)
  p = w1*p1 + (1-w1)*p2

  if (!lower.tail) p = 1 - p
  if (log.p) p = log(p)
  p
}

