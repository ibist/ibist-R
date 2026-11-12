test_that("Katz risk-ratio intervals reproduce published SAS results", {
  ## SAS/STAT 14.2 FREQ, Output 41.5.3, p. 2861.
  tab <- matrix(c(11, 4, 2, 6), nrow = 2, byrow = TRUE)
  result <- rr.2x2.ci(tab, method = "katz")
  expect_equal(round(unname(result$estimate), 4), 2.9333)
  expect_equal(round(as.numeric(result$conf.int), 4), c(0.8502, 10.1204))
  other <- rr.2x2.ci(tab[, 2:1], method = "katz")
  expect_equal(round(as.numeric(other$conf.int), 4), c(0.1403, 0.9009))
})

test_that("Koopman limits satisfy an independently optimized Pearson score", {
  x <- c(38, 267)
  n <- c(100, 984)
  tab <- cbind(x, n - x)
  result <- rr.2x2.ci(tab)
  for (ratio in as.numeric(result$conf.int)) {
    loglik <- function(q) {
      sum(dbinom(x, n, c(ratio * q, q), log = TRUE))
    }
    fit <- optimize(loglik, c(0, min(1, 1 / ratio)), maximum = TRUE,
                    tol = 1e-12)
    p <- c(ratio * fit$maximum, fit$maximum)
    statistic <- sum((x - n * p)^2 / (n * p * (1 - p)))
    expect_equal(statistic, qchisq(0.95, 1), tolerance = 1e-6)
  }
  expect_equal(unname(result$estimate), 0.38 / (267 / 984))
})

test_that("risk ratios preserve row-reversal and one-sided relationships", {
  for (x in list(c(0, 0), c(0, 5), c(3, 0), c(3, 5), c(2, 4))) {
    tab <- cbind(x, c(3, 5) - x)
    result <- rr.2x2.ci(tab, conf.level = 0.9)
    opposite <- rr.2x2.ci(tab[2:1, ], conf.level = 0.9)
    expect_equal(as.numeric(result$conf.int),
                 1 / rev(as.numeric(opposite$conf.int)), tolerance = 1e-7)
    lower <- rr.2x2.ci(tab, alternative = "greater")
    upper <- rr.2x2.ci(tab, alternative = "less")
    expect_equal(as.numeric(lower$conf.int), c(result$conf.int[1], Inf))
    expect_equal(as.numeric(upper$conf.int), c(0, result$conf.int[2]))
  }
  empty <- rr.2x2.ci(matrix(c(0, 3, 0, 5), 2, byrow = TRUE))
  expect_true(is.na(empty$estimate))
  expect_equal(as.numeric(empty$conf.int), c(0, Inf))
})

test_that("risk-ratio methods validate table orientation and zero counts", {
  expect_error(rr.2x2.ci(c(1, 2, 3, 4)), "2 x 2")
  expect_error(rr.2x2.ci(matrix(1, 2, 3)), "2 x 2")
  expect_error(rr.2x2.ci(matrix(c(0, 0, 2, 3), 2, byrow = TRUE)),
               "integer counts")
  expect_error(rr.2x2.ci(matrix(c(0.5, 2, 3, 4), 2)), "integer counts")
  tab <- matrix(c(0, 3, 2, 5), 2, byrow = TRUE)
  expect_error(rr.2x2.ci(tab, method = "katz"), "positive event counts")
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  result <- rr.2x2.ci(tab, method = c("katz", "koopman"))
  expect_equal(rownames(result$conf.int), c("katz", "koopman"))
  expect_s3_class(result, "ibist_ci")
})
