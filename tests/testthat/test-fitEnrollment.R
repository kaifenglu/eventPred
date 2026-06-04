testthat::test_that("poisson enrollment model", {
  fit1 <- fitEnrollment(
    df = interimData1, enroll_model = "poisson",
    showplot = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))

  b <- list(theta = -0.482757,
            vtheta = 0.004464,
            aic = 666.2750,
            bic = 669.6867)

  testthat::expect_equal(a, b)
})


testthat::test_that("time-decay enrollment model", {
  fit1 <- fitEnrollment(
    df = interimData1, enroll_model = "time-decay",
    showplot = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))

  b <- list(theta = c(-5.394735, -5.958500),
            vtheta = matrix(c(0.043052, 0.144519,
                              0.144519, 0.541253),
                            2, 2, byrow = TRUE),
            aic = 609.0973,
            bic = 615.9206)

  testthat::expect_equal(a, b)
})


testthat::test_that("spline enrollment model", {
  fit1 <- fitEnrollment(
    df = interimData1, enroll_model = "b-spline",
    nknots = 1, showplot = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))

  b <- list(theta = c(-2.410888, -0.782993, -0.854621, 0.260845, -0.000514),
            vtheta = matrix(c(
              0.406455, -0.214008,  0.166808, -0.066814,  0.022651,
              -0.214008,  0.293015, -0.264660,  0.113540, -0.039881,
              0.166808, -0.264660,  0.354398, -0.177616,  0.067117,
              -0.066814,  0.113540, -0.177616,  0.131576, -0.058527,
              0.022651, -0.039881,  0.067117, -0.058527,  0.063791),
              5, 5, byrow = TRUE),
            aic = 614.6706,
            bic = 631.7288)

  testthat::expect_equal(a, b)
})


testthat::test_that("piecewise poisson enrollment model", {
  fit1 <- fitEnrollment(
    df = interimData1, enroll_model = "piecewise poisson",
    accrualTime = c(0, 140), showplot = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))

  b <- list(theta = c(-1.180442, -0.208675),
            vtheta = matrix(c(0.023256, 0.000000,
                              0.000000, 0.005525),
                            2, 2, byrow = TRUE),
            aic = 629.0583,
            bic = 635.8816)

  testthat::expect_equal(a, b)
})



