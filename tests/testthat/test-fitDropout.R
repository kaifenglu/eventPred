testthat::test_that("exponential dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "exponential",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(as.numeric(fit1$fit$vtheta), 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))

  b <- list(theta = -9.248655,
            vtheta = 0.0833333,
            aic = 247.9677,
            bic = 251.6715)

  testthat::expect_equal(a, b)
})


testthat::test_that("weibull dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "weibull",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(8.737596, -0.189425),
            vtheta = matrix(c(0.4338363, 0.1552454,
                              0.1552454, 0.0639657),
                            2, 2, byrow = TRUE),
            aic = 249.4456,
            bic = 256.8532)

  testthat::expect_equal(a, b)
})


testthat::test_that("log-logistic dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "log-logistic",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(8.682621, -0.203227),
            vtheta = matrix(c(0.4236892, 0.1525186,
                              0.1525186, 0.0637534),
                            2, 2, byrow = TRUE),
            aic = 249.5160,
            bic = 256.9236)

  testthat::expect_equal(a, b)
})


testthat::test_that("log-normal dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "log-normal",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(9.806748, 0.764730),
            vtheta = matrix(c(0.8172429, 0.1943497,
                              0.1943497, 0.0515406),
                            2, 2, byrow = TRUE),
            aic = 250.8187,
            bic = 258.2262)

  testthat::expect_equal(a, b)
})


testthat::test_that("piecewise exponential dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "piecewise exponential",
    piecewiseDropoutTime = c(0, 180),
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(-9.350537, -9.193537),
            vtheta = matrix(c(0.250, 0.000,
                              0.000, 0.125),
                            2, 2, byrow = TRUE),
            aic = 249.9009,
            bic = 257.3085)

  testthat::expect_equal(a, b)
})


testthat::test_that("model averaging dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "model averaging",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4),
            w1 = round(fit1$fit$w1, 6))

  b <- list(theta = c(8.737596, -0.189425, 9.806748, 0.764730),
            vtheta = matrix(c(0.4338363, 0.1552454, 0.0000000, 0.0000000,
                              0.1552454, 0.0639657, 0.0000000, 0.0000000,
                              0.0000000, 0.0000000, 0.8172429, 0.1943497,
                              0.0000000, 0.0000000, 0.1943497, 0.0515406),
                            4, 4, byrow = TRUE),
            aic = 249.9053,
            bic = 257.3129,
            w1 = 0.665195)

  testthat::expect_equal(a, b)
})


testthat::test_that("spline dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "spline",
    k = 1, scale = "hazard",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4),
            knots = round(as.numeric(fit1$fit$knots), 6))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(-8.324153, 0.669893, -0.158284),
            vtheta = matrix(c( 6.6251945, -1.3644615, -0.2816630,
                              -1.3644615,  0.3023492,  0.0738798,
                              -0.2816630,  0.0738798,  0.0250222),
                            3, 3, byrow = TRUE),
            aic = 250.5939,
            bic = 261.7053,
            knots = c(2.302585, 6.032853, 6.501290))

  testthat::expect_equal(a, b)
})


testthat::test_that("cox dropout model", {
  fit1 <- fitDropout(
    df = interimData2, dropout_model = "cox",
    m = 5, showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 7),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4),
            piecewiseDropoutTime = fit1$fit$piecewiseDropoutTime)

  b <- list(theta = c(
    -7.986165,-9.952230,-9.205730,-8.337109,-9.787066,-10.217495,-7.965546,
    -7.135687,-9.019785,-8.931155,-8.359837,-8.217439),
    vtheta = diag(12),
    aic = 258.2305,
    bic = 302.6759,
    piecewiseDropoutTime = c(
      0, 10, 93, 135, 153, 239, 408, 426, 434, 491, 553, 601, 666))

  testthat::expect_equal(a, b)
})


testthat::test_that("cox dropout model validates m_dropout using distinct dropout times", {
  df <- interimData2
  df$time[df$dropout == 1] <- 20

  testthat::expect_error(
    fitDropout(
      df = df,
      dropout_model = "cox",
      m_dropout = 2,
      showplot = FALSE,
      by_treatment = FALSE),
    "distinct observed dropout times"
  )
})
