testthat::test_that("event prediction at design stage", {
  set.seed(1000)

  enroll_prior <- list(
    model = "piecewise poisson",
    accrualTime = seq(0, 8)*30.4375,
    theta = log(26/9*seq(1, 9)/30.4375),
    vtheta = diag(9)*1e-8)

  pred1 <- predictEnrollment(
    target_n = 300,
    enroll_fit = enroll_prior,
    pilevel = 0.90, nreps = 200,
    showsummary = FALSE, showplot = FALSE)

  event_prior <- list(
    model = "piecewise exponential",
    theta = log(c(0.0533, 0.0309)/30.4375),
    vtheta = diag(2)*1e-8,
    piecewiseSurvivalTime = c(0,6)*30.4375)

  dropout_prior <- list(
    model = "exponential",
    theta = log((-log(1-0.05)/12)/30.4375),
    vtheta = 1e-8)

  pred2 <- predictEvent(
    target_d = 200,
    newSubjects = pred1$newSubjects,
    event_fit = event_prior,
    dropout_fit = dropout_prior,
    pilevel = 0.90, nreps = 200,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred2$event_pred_day)
  b <- c(1349, 1208, 1511)

  testthat::expect_equal(a, b)
})


testthat::test_that("event prediction during enrollment stage", {
  set.seed(2000)

  enroll_fits <- fitEnrollment(
    df = interimData1, enroll_model = "b-spline",
    nknots = 1, showplot = FALSE)

  pred1 <- predictEnrollment(
    df = interimData1,
    target_n = 300,
    enroll_fit = enroll_fits$fit,
    lags = 30,
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)


  event_fits <- fitEvent(
    df = interimData1,
    event_model = "weibull",
    showplot = FALSE)

  dropout_fits <- fitDropout(
    df = interimData1,
    dropout_model = "exponential",
    showplot = FALSE)

  pred2 <- predictEvent(
    df = interimData1, target_d = 200,
    newSubjects = pred1$newSubjects,
    event_fit = event_fits$fit,
    dropout_fit = dropout_fits$fit,
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred2$event_pred_day)
  b <- c(1157, 835, 1930)

  testthat::expect_equal(a, b)
})


testthat::test_that("event prediction after enrollment completion", {
  set.seed(3000)

  event_fits <- fitEvent(
    df = interimData2,
    event_model = "piecewise exponential",
    piecewiseSurvivalTime = c(0, 140, 352),
    showplot = FALSE)

  dropout_fits <- fitDropout(
    df = interimData2,
    dropout_model = "exponential",
    showplot = FALSE)

  pred1 <- predictEvent(
    df = interimData2, target_d = 200,
    event_fit = event_fits$fit,
    dropout_fit = dropout_fits$fit,
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$event_pred_day)
  b <- c(1093, 1047, 1168)

  testthat::expect_equal(a, b)
})


testthat::test_that("cox event/dropout prediction handles oversized tail windows", {
  set.seed(7100)

  event_fits <- fitEvent(
    df = interimData2,
    event_model = "cox",
    m = 3,
    showplot = FALSE)

  dropout_fits <- fitDropout(
    df = interimData2,
    dropout_model = "cox",
    m_dropout = 3,
    showplot = FALSE)

  pred <- predictEvent(
    df = interimData2,
    target_d = 200,
    event_fit = event_fits$fit,
    m = 999,
    dropout_fit = dropout_fits$fit,
    m_dropout = 999,
    pilevel = 0.90,
    nreps = 5,
    showsummary = FALSE,
    showplot = FALSE)

  testthat::expect_true(is.list(pred))
  testthat::expect_true(all(c("event_pred_day", "newEvents") %in% names(pred)))
})
