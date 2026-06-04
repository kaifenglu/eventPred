testthat::test_that("summarize observed data during enrollment", {
  observed1 <- summarizeObserved(
    df = interimData1,
    to_predict = "enrollment and event",
    showplot = FALSE,
    by_treatment = FALSE)

  a <- data.frame(trialsdt = observed1$trialsdt,
                  cutoffdt = observed1$cutoffdt,
                  t0 = observed1$t0,
                  n0 = observed1$n0,
                  d0 = observed1$d0,
                  c0 = observed1$c0,
                  r0 = observed1$r0)

  b <- data.frame(trialsdt = as.Date("2018-03-01"),
                  cutoffdt = as.Date("2019-03-01"),
                  t0 = 366, n0 = 224,
                  d0 = 42, c0 = 1, r0 = 181)

  testthat::expect_equal(as.numeric(a$cutoffdt - a$trialsdt + 1), a$t0)
  testthat::expect_equal(a$n0, a$d0 + a$c0 + a$r0)
  testthat::expect_equal(a, b)
})


testthat::test_that("summarize observed data after enrollment", {
  observed1 <- summarizeObserved(
    df = interimData2,
    to_predict = "event only",
    showplot = FALSE,
    by_treatment = FALSE)

  a <- data.frame(trialsdt = observed1$trialsdt,
                  cutoffdt = observed1$cutoffdt,
                  t0 = observed1$t0,
                  n0 = observed1$n0,
                  d0 = observed1$d0,
                  c0 = observed1$c0,
                  r0 = observed1$r0)

  b <- data.frame(trialsdt = as.Date("2018-03-01"),
                  cutoffdt = as.Date("2020-10-21"),
                  t0 = 966, n0 = 300,
                  d0 = 183, c0 = 12, r0 = 105)

  testthat::expect_equal(as.numeric(a$cutoffdt - a$trialsdt + 1), a$t0)
  testthat::expect_equal(a$n0, a$d0 + a$c0 + a$r0)
  testthat::expect_equal(a, b)
})

