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
#' @param q The number of elements in the vector of covariates
#'   (excluding the intercept).
#' @param x The vector of covariates (including the intercept).
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
    xbeta = sum(x[-1] * theta[(J+1):(J+q)])
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
#' @param q The number of elements in the vector of covariates
#'   (excluding the intercept).
#' @param x The vector of covariates (including the intercept).
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
    xbeta = sum(x[-1] * theta[(J+1):(J+q)])
    cumhaz0 = cumhaz0*exp(-xbeta)
  }

  # find the interval containing time
  j1 = findInterval(cumhaz0, psum)
  t = tcut[j1] + (cumhaz0 - psum[j1])/lambda[j1]
  t
}


#' @title Profile log likelihood for piecewise exponential regression
#' @description Obtains the profile log likelihood value for piecewise
#' exponential regression.
#'
#' @param beta The regression coefficients with respect the covariates.
#' @param time The survival time.
#' @param event The event indicator.
#' @param J The number of time intervals.
#' @param tcut A vector that specifies the endpoints of time intervals
#'   for the baseline piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2 event
#'   intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param q The number of columns of the covariates matrix
#'   (excluding the intercept).
#' @param x The covariates matrix (including the intercept).
#'
#' @return The profile log likelihood value for piecewise exponential
#' regression.
#'
#' @keywords internal
#'
#' @export
#'
pllik_pwexp <- function(beta, time, event, J, tcut, q, x) {
  if (length(tcut) == J) tcut = c(tcut, Inf)

  xbeta = as.numeric(x[,-1] %*% beta)

  d = rep(0,J)
  ex0 = rep(0,J)
  for (j in 1:J) {
    d[j] = sum(event * (time >= tcut[j]) * (time < tcut[j+1]))
    ex0[j] = sum(pmax(0, pmin(time, tcut[j+1]) - tcut[j])*exp(xbeta))
  }

  sum(event*xbeta) - sum(d*log(ex0))
}


#' @title Piecewise exponential regression
#' @description Obtains the maximum likelihood estimates for piecewise
#' exponential regression.
#'
#' @param time The survival time.
#' @param event The event indicator.
#' @param J The number of time intervals.
#' @param tcut A vector that specifies the endpoints of time intervals
#'   for the baseline piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2 event
#'   intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param q The number of columns of the covariates matrix
#'   (exluding the intercept).
#' @param x The covariates matrix (including the intercept).
#'
#' @return The maximum likelihood estimates and the associated
#' covariance matrix, AIC and BIC.
#'
#' @export
#'
pwexpreg <- function(time, event, J, tcut, q, x) {
  # get the variable name as a character string
  if (grepl("dropout", deparse(substitute(event)), ignore.case = TRUE)) {
    variable_name = "dropout"
  } else {
    variable_name = "event"
  }

  if (length(tcut) == J) tcut = c(tcut, Inf)

  n = length(time)

  d = rep(0,J)  # number of events in each interval
  exmat = matrix(0,n,J)
  for (j in 1:J) {
    d[j] = sum(event * (time >= tcut[j]) * (time < tcut[j+1]))
    exmat[,j] = pmax(0, pmin(time, tcut[j+1]) - tcut[j])
  }

  if (any(d == 0)) {
    stop(paste("The number of", paste0(variable_name, "s"),
               "must be >= 1 in each interval",
               "to fit a piecewise exponential model."))
  }

  if (q == 0) { # piecewise exponential with no covariates
    ex = colSums(exmat)  # exposure in each interval

    theta = log(d/ex)
    if (J > 1) {
      vtheta = diag(1/d)
    } else {
      vtheta = 1/d*diag(1)
    }

    llik = sum(theta*d - exp(theta)*ex)

    fit <- list(model = "Piecewise exponential",
                theta = theta,
                vtheta = vtheta,
                aic = -2*llik + 2*J,
                bic = -2*llik + J*log(n))
  } else { # piecewise exponential with covariates
    beta0 = rep(0,q)
    if (q > 1) {
      opt1 <- optim(beta0, pllik_pwexp, gr = NULL,
                    time, event, J, tcut, q, x,
                    control = c(fnscale = -1), hessian = TRUE)
    } else {
      sdx = sd(x[,-1])
      opt1 <- optim(beta0, pllik_pwexp, gr = NULL,
                    time, event, J, tcut, q, x,
                    method = "Brent",
                    lower = beta0 - 6*sdx, upper = beta0 + 6*sdx,
                    control = c(fnscale = -1), hessian = TRUE)
    }

    beta = opt1$par
    vbeta = solve(-opt1$hessian)

    xbeta = as.numeric(x[,-1] %*% beta)
    exbeta = exp(xbeta)

    ex0 = rep(0,J)
    for (j in 1:J) {
      ex0[j] = sum(exmat[,j]*exbeta)
    }

    ex1 = matrix(0,J,q)
    for (j in 1:J) {
      for (k in 1:q) {
        ex1[j,k] = sum(exmat[,j]*x[,k+1]*exbeta)
      }
    }

    rx1 = matrix(0,J,q)
    for (j in 1:J) {
      rx1[j,] = ex1[j,]/ex0[j]
    }

    vtheta = matrix(0, J+q, J+q)
    for (j in 1:J) {
      for (k in 1:J) {
        vtheta[j,k] = 1/d[j]*(j==k) + rx1[j,] %*% vbeta %*% rx1[k,]
      }
    }

    for (j in 1:J) {
      for (k in 1:q) {
        vtheta[j,J+k] = -rx1[j,] %*% vbeta[,k]
        vtheta[J+k,j] = vtheta[j,J+k]
      }
    }

    for (k1 in 1:q) {
      for (k2 in 1:q) {
        vtheta[J+k1,J+k2] = vbeta[k1,k2]
      }
    }

    alpha = log(d/ex0)
    theta = c(alpha, beta)
    llik = sum(alpha*d) + sum(event*xbeta) - sum(exp(alpha)*ex0)

    fit <- list(model = "Piecewise exponential",
                theta = theta,
                vtheta = vtheta,
                aic = -2*llik + 2*(J+q),
                bic = -2*llik + (J+q)*log(n))
  }

  if (variable_name == "event") {
    fit$piecewiseSurvivalTime = tcut[1:J]
  } else {
    fit$piecewiseDropoutTime = tcut[1:J]
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
#' @param q The number of elements in the vector of covariates
#'   (excluding the intercept).
#' @param x The vector of covariates (including the intercept).
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
  scale = exp(sum(x * theta[1:(q+1)]))
  meanlog = sum(x * theta[(q+3):(2*q+3)])
  sdlog = exp(theta[2*q+4])

  p1 = pweibull(t, shape, scale)
  p2 = plnorm(t, meanlog, sdlog)
  p = w1*p1 + (1-w1)*p2

  if (!lower.tail) p = 1 - p
  if (log.p) p = log(p)
  p
}

