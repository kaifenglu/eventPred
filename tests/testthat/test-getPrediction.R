testthat::test_that("event prediction at design stage", {
  set.seed(1000)

  enroll_prior <- list(
    model = "piecewise poisson",
    accrualTime = seq(0, 8)*30.4375,
    theta = log(26/9*seq(1, 9)/30.4375),
    vtheta = diag(9)*1e-8)

  event_prior <- list(
    model = "piecewise exponential",
    theta = log(c(0.0533, 0.0309)/30.4375),
    vtheta = diag(2)*1e-8,
    piecewiseSurvivalTime = c(0,6)*30.4375)

  dropout_prior <- list(
    model = "exponential",
    theta = log((-log(1-0.05)/12)/30.4375),
    vtheta = 1e-8)

  pred1 <- getPrediction(
    df = NULL,
    to_predict = "enrollment and event",
    target_n = 300,
    target_d = 200,
    enroll_prior = enroll_prior,
    event_prior = event_prior,
    dropout_prior = dropout_prior,
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$event_pred$event_pred_day)
  b <- c(1352, 1180, 1527)

  testthat::expect_equal(a, b)
})


testthat::test_that("design-stage enrollment respects multi-arm allocation", {
  set.seed(4000)

  pred1 <- getPrediction(
    df = NULL,
    to_predict = "enrollment only",
    target_n = 12,
    enroll_prior = list(
      model = "poisson",
      theta = log(0.2),
      vtheta = 1e-8),
    ngroups = 2,
    alloc = c(2, 1),
    treatment_label = c("Control", "Experimental"),
    pilevel = 0.90, nreps = 3,
    showsummary = FALSE, showplot = FALSE)

  counts <- pred1$subject_data[
    , as.list(table(factor(treatment_description,
                           levels = c("Control", "Experimental")))),
    by = "draw"]

  testthat::expect_equal(counts$Control, rep(8L, 3))
  testthat::expect_equal(counts$Experimental, rep(4L, 3))
})


testthat::test_that("event prediction during enrollment phase", {
  set.seed(2000)

  pred1 <- getPrediction(
    df = interimData1,
    to_predict = "enrollment and event",
    target_n = 300,
    target_d = 200,
    enroll_model = "time-decay",
    event_model = "weibull",
    dropout_model = "exponential",
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$event_pred$event_pred_day)
  b <- c(1231, 835, 2018)

  testthat::expect_equal(a, b)
})


testthat::test_that("event prediction after enrollment completion", {
  set.seed(3000)

  pred1 <- getPrediction(
    df = interimData2, to_predict = "event only",
    target_d = 200,
    event_model = "cox", m = 5,
    dropout_model = "exponential",
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$event_pred$event_pred_day)
  b <- c(1027, 1003, 1099)

  testthat::expect_equal(a, b)
})


make_covariate_data <- function(df) {
  df$age <- 50 + seq_len(nrow(df)) %% 12
  df$sex <- factor(rep(c("F", "M"), length.out = nrow(df)))
  df
}


testthat::test_that("enrollment and event prediction assembles prior outputs with covariates", {
  set.seed(4000)

  df <- make_covariate_data(interimData1)

  enroll_prior <- list(
    model = "time-decay",
    theta = c(log(0.8), log(0.02)),
    vtheta = diag(c(0.5, 0.5)))

  event_prior <- list(
    model = "exponential",
    theta = log(1/400),
    vtheta = 0.25)

  dropout_prior <- list(
    model = "exponential",
    theta = log(1/800),
    vtheta = 0.25)

  event_prior_with_covariates <- list(
    model = "exponential",
    theta = c(log(1/400), 0.02, -0.15),
    vtheta = diag(c(0.25, 0.01, 0.01)))

  dropout_prior_with_covariates <- list(
    model = "exponential",
    theta = c(log(1/800), 0.01, -0.10),
    vtheta = diag(c(0.25, 0.01, 0.01)))

  pred <- getPrediction(
    df = df,
    to_predict = "enrollment and event",
    target_n = 230,
    target_d = 200,
    enroll_model = "time-decay",
    enroll_prior = enroll_prior,
    event_model = "exponential",
    event_prior = event_prior,
    covariates_event = c("age", "sex"),
    event_prior_with_covariates = event_prior_with_covariates,
    dropout_model = "exponential",
    dropout_prior = dropout_prior,
    nreps = 5,
    showsummary = FALSE,
    showplot = FALSE,
    generate_plot = FALSE)

  testthat::expect_equal(pred$stage, "Real-time before enrollment completion")
  testthat::expect_true(all(c(
    "enroll_prior", "enroll_post",
    "event_fit", "event_fit_with_covariates",
    "event_prior", "event_post",
    "event_prior_with_covariates", "event_post_with_covariates",
    "dropout_fit", "dropout_prior", "dropout_post"
  ) %in% names(pred)))
  testthat::expect_false(any(c(
    "dropout_fit_with_covariates",
    "dropout_prior_with_covariates",
    "dropout_post_with_covariates"
  ) %in% names(pred)))
  testthat::expect_true(all(c("age", "sex") %in% names(pred$subject_data)))
})


testthat::test_that("event only prediction uses covariate-specific outputs", {
  set.seed(5000)

  df <- make_covariate_data(interimData2)

  event_prior_with_covariates <- list(
    model = "exponential",
    theta = c(log(1/400), 0.02, -0.15),
    vtheta = diag(c(0.25, 0.01, 0.01)))

  dropout_prior_with_covariates <- list(
    model = "exponential",
    theta = c(log(1/800), 0.01, -0.10),
    vtheta = diag(c(0.25, 0.01, 0.01)))

  pred <- getPrediction(
    df = df,
    to_predict = "event only",
    target_d = 200,
    event_model = "exponential",
    covariates_event = c("age", "sex"),
    event_prior_with_covariates = event_prior_with_covariates,
    dropout_model = "exponential",
    covariates_dropout = c("age", "sex"),
    dropout_prior_with_covariates = dropout_prior_with_covariates,
    nreps = 5,
    showsummary = FALSE,
    showplot = FALSE,
    generate_plot = FALSE)

  testthat::expect_equal(pred$stage, "Real-time after enrollment completion")
  testthat::expect_true(all(c(
    "event_fit_with_covariates", "event_prior_with_covariates",
    "event_post_with_covariates", "dropout_fit_with_covariates",
    "dropout_prior_with_covariates", "dropout_post_with_covariates"
  ) %in% names(pred)))
  testthat::expect_false(any(c(
    "event_fit", "event_prior", "event_post",
    "dropout_fit", "dropout_prior", "dropout_post"
  ) %in% names(pred)))
  testthat::expect_true(all(c("age", "sex") %in% names(pred$subject_data)))
})


testthat::test_that("design-stage multi-arm allocation is reflected in simulated subjects", {
  set.seed(6000)

  enroll_prior <- list(
    model = "poisson",
    theta = log(0.2),
    vtheta = 1e-8)

  event_prior <- list(
    list(model = "exponential", theta = log(1/200), vtheta = 1e-8),
    list(model = "exponential", theta = log(1/300), vtheta = 1e-8))

  dropout_prior <- list(
    list(model = "exponential", theta = log(1/1000), vtheta = 1e-8),
    list(model = "exponential", theta = log(1/1000), vtheta = 1e-8))

  pred <- getPrediction(
    df = NULL,
    to_predict = "enrollment and event",
    target_n = 8,
    target_d = 3,
    enroll_prior = enroll_prior,
    event_prior = event_prior,
    dropout_prior = dropout_prior,
    ngroups = 2,
    alloc = c(1, 3),
    treatment_label = c("Control", "Experimental"),
    nreps = 4,
    showsummary = FALSE,
    showplot = FALSE,
    generate_plot = FALSE,
    fix_parameter = TRUE)

  count_by_draw <- with(
    pred$enroll_pred$newSubjects,
    table(draw, treatment_description))

  testthat::expect_equal(pred$stage, "Design stage")
  testthat::expect_equal(unname(count_by_draw[, "Control"]), rep(2, 4))
  testthat::expect_equal(unname(count_by_draw[, "Experimental"]), rep(6, 4))
  testthat::expect_true(all(c("event_prior", "dropout_prior") %in% names(pred)))
})


testthat::test_that("design-stage event-only prediction is rejected", {
  event_prior <- list(
    model = "exponential",
    theta = log(1/200),
    vtheta = 1e-8)

  testthat::expect_error(
    getPrediction(
      df = NULL,
      to_predict = "event only",
      target_d = 3,
      event_prior = event_prior,
      nreps = 2,
      showsummary = FALSE,
      showplot = FALSE,
      generate_plot = FALSE),
    "Design-stage event prediction requires enrollment prediction"
  )
})
