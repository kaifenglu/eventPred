testthat::test_that("exponential event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "exponential",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(as.numeric(fit1$fit$vtheta), 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))

  b <- list(theta = -6.524076,
            vtheta = 0.005464,
            aic = 2755.8117,
            bic = 2759.5155)

  testthat::expect_equal(a, b)
})


testthat::test_that("weibull event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "weibull",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(6.550712, 0.092555),
            vtheta = matrix(c(0.007007, 0.001375,
                              0.001375, 0.004380),
                            2, 2, byrow = TRUE),
            aic = 2755.7842,
            bic = 2763.1918)

  testthat::expect_equal(a, b)
})


testthat::test_that("log-logistic event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "log-logistic",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(6.128323, -0.093219),
            vtheta = matrix(c(0.009133, 0.001056,
                              0.001056, 0.004095),
                            2, 2, byrow = TRUE),
            aic = 2760.6224,
            bic = 2768.0300)

  testthat::expect_equal(a, b)
})


testthat::test_that("log-normal event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "log-normal",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(6.116114, 0.510694),
            vtheta = matrix(c(0.011718, 0.001850,
                              0.001850, 0.003156),
                            2, 2, byrow = TRUE),
            aic = 2769.5140,
            bic = 2776.9216)

  testthat::expect_equal(a, b)
})


testthat::test_that("piecewise exponential event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "piecewise exponential",
    piecewiseSurvivalTime = c(0, 180),
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(-6.446372, -6.572499),
            vtheta = matrix(c(0.013699, 0.000000,
                              0.000000, 0.009091),
                            2, 2, byrow = TRUE),
            aic = 2757.1200,
            bic = 2764.5276)

  testthat::expect_equal(a, b)
})


testthat::test_that("model averaging event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "model averaging",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4),
            w1 = round(fit1$fit$w1, 6))

  b <- list(theta = c(6.550712, 0.092555, 6.116114, 0.510694),
            vtheta = matrix(c(0.007007, 0.001375, 0.000000, 0.000000,
                              0.001375, 0.004380, 0.000000, 0.000000,
                              0.000000, 0.000000, 0.011718, 0.001850,
                              0.000000, 0.000000, 0.001850, 0.003156),
                            4, 4, byrow = TRUE),
            aic = 2755.7985,
            bic = 2763.2061,
            w1 = 0.998957)

  testthat::expect_equal(a, b)
})


testthat::test_that("spline event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "spline",
    k = 1, scale = "hazard",
    showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4),
            knots = round(as.numeric(fit1$fit$knots), 6))
  dimnames(a$vtheta) <- NULL

  b <- list(theta = c(-5.997161, 0.919534, 0.000487),
            vtheta = matrix(c( 0.522662, -0.136872, -0.006925,
                              -0.136872,  0.038342,  0.002105,
                              -0.006925,  0.002105,  0.000128),
                            3, 3, byrow = TRUE),
            aic = 2757.7835,
            bic = 2768.8948,
            knots = c(0.000000, 5.480639, 6.845880))

  testthat::expect_equal(a, b)
})


testthat::test_that("cox event model", {
  fit1 <- fitEvent(
    df = interimData2, event_model = "cox",
    m = 5, showplot = FALSE, by_treatment = FALSE)

  a <- list(theta = round(fit1$fit$theta, 6),
            vtheta = round(fit1$fit$vtheta, 6),
            aic = round(fit1$fit$aic, 4),
            bic = round(fit1$fit$bic, 4),
            piecewiseSurvivalTime = fit1$fit$piecewiseSurvivalTime)

  b <- list(theta = c(
    -5.703782,-6.393591,-5.002264,-6.383507,-6.380123,-6.376727,-6.775366,
    -6.771936,-5.669881,-5.666427,-7.272398,-5.657735,-5.652489,-6.747587,
    -5.645447,-5.641907,-4.943423,-6.324359,-6.320768,-6.874971,-6.715383,
    -6.997593,-4.910812,-5.598422,-5.998313,-5.871169,-6.674561,-5.568335,
    -6.253829,-6.943122,-6.467307,-5.545177,-5.539297,-7.927324,-5.525453,
    -5.521461,-6.210600,-5.513429,-5.509388,-6.603944,-7.110696,-6.595781,
    -5.493061,-5.488938,-7.094235,-6.579251,-7.779049,-6.165418,-7.543273,
    -6.558198,-6.148468,-6.144186,-7.234177,-5.438079,-7.379632,-5.429346,
    -5.424950,-7.029973,-6.802395,-6.510258,-6.793466,-6.095825,-5.398163,
    -7.339538,-5.792246,-4.684438,-6.756932,-7.930925,-6.459904,-5.356586,
    -6.450470,-6.445720,-6.035481,-7.129298,-7.635304,-4.627415,-6.003881,
    -7.249215,-4.602661,-5.981414,-7.229114,-7.475339,-5.273000,-6.874709,
    -5.950643,-6.350886,-5.247024,-6.628041,-6.335054,-5.924256,-6.612041,
    -5.220356,-7.412160,-5.612208,-5.198497,-6.984716,-6.285998,-6.280396,
    -6.562444,-6.556778,-5.164786,-6.950815,-6.251904,-7.226936,-5.141664,
    -5.538309,-6.730412,-7.191429,-6.204558,-5.793014,-8.312626,-7.720462,
    -7.371489,-6.161207,-7.352441,-6.429719,-5.036953,-6.416732,-6.273338,
    -7.408531,-6.102559,-4.300670,-6.369901,-7.742402,-6.349139,-4.955827,
    -6.047372,-6.040255,-4.934474,-6.025866,-7.477604,-4.897840,-6.276643,
    -7.919356,-5.505332,-5.493045,-4.779123,-4.770685,-6.938021,-5.365976,
    -5.337538,-6.654153,-4.564348,-5.929589,-7.162397,-6.995766,-5.023881,
    -7.637234,-7.022868,-5.010635,-4.969813,-6.182085,-5.634790,-7.049255,
    -5.768321,-6.664409,-5.023881),
    vtheta = diag(c(
      1.000000,1.000000,0.500001,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,0.500002,1.000000,1.000000,
      1.000000,1.000000,0.500002,1.000000,1.000000,0.500002,1.000000,
      0.500002,0.500002,1.000000,0.500002,0.333336,1.000000,0.333337,
      1.000000,1.000000,0.500002,1.000000,0.500002,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,0.500003,0.500003,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,0.500003,0.500003,
      1.000000,0.500003,1.000000,1.000000,1.000000,1.000000,0.500003,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,0.500004,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      0.500004,0.500004,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,0.500005,
      1.000000,1.000000,0.500006,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,0.500008,1.000000,1.000000,0.500011,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
      1.000000,1.000000,1.000000)),
    aic = 2930.3483,
    bic = 3511.8421,
    piecewiseSurvivalTime = c(
      0,   1,   3,   4,   6,   8,  10,  13,  16,  17,
      18,  23,  25,  26,  29,  30,  31,  32,  34,
      36,  43,  46,  54,  55,  56,  59,  63,  66,  69,
      71,  75,  80,  81,  83,  94,  95,  96,  98,
      99, 100, 103, 108, 111, 112, 113, 118, 121, 131,
      133, 141, 144, 146, 148, 154, 155, 162, 163,
      164, 169, 173, 176, 180, 182, 183, 190, 193, 194,
      198, 211, 214, 215, 218, 221, 223, 229, 239,
      240, 244, 251, 252, 254, 261, 270, 271, 281, 283,
      286, 287, 291, 294, 296, 300, 301, 310, 313,
      314, 320, 323, 326, 330, 334, 335, 341, 344, 352,
      353, 356, 366, 374, 377, 379, 404, 418, 428,
      431, 441, 445, 446, 450, 457, 468, 471, 472, 476,
      492, 496, 497, 500, 503, 504, 507, 520, 521,
      525, 547, 549, 553, 554, 555, 574, 576, 578, 586,
      587, 591, 606, 620, 622, 656, 678, 681, 684,
      695, 702, 750, 766, 864, 940))

  testthat::expect_equal(a, b)
})


testthat::test_that("cox event model validates m using distinct event times", {
  df <- interimData2
  df$time[df$event == 1] <- 10

  testthat::expect_error(
    fitEvent(
      df = df,
      event_model = "cox",
      m = 2,
      showplot = FALSE,
      by_treatment = FALSE),
    "distinct observed event times"
  )
})
