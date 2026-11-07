test_that("Pearson wraps chisq.test with explicit correction defaults", {
  tab <- matrix(c(12, 8, 5, 15), nrow = 2, byrow = TRUE,
                dimnames = list(group = c("A", "B"), outcome = c("yes", "no")))
  for (correct in c(FALSE, TRUE)) {
    actual <- indep.test(tab, correct = correct)
    expected <- stats::chisq.test(tab, correct = correct)
    expect_identical(actual, expected)
  }
  expect_identical(indep.test(tab), indep.test(tab, method = "pearson"))

  balanced <- matrix(10, nrow = 2, ncol = 2)
  balanced[1, 1] <- 11
  expect_equal(unname(indep.test(balanced, correct = TRUE)$statistic), 0)
})

test_that("likelihood ratio matches a hand calculation with a zero cell", {
  tab <- matrix(c(0, 20, 10, 10), nrow = 2, byrow = TRUE)
  result <- indep.test(tab, method = "lr")
  # Expected counts are (5, 15) in each row; 0 log(0/E) contributes zero.
  expected.statistic <- 2 * (20 * log(4 / 3) + 10 * log(2) +
                              10 * log(2 / 3))
  expect_s3_class(result, "htest")
  expect_identical(names(result$statistic), "G-squared")
  expect_equal(unname(result$statistic), expected.statistic,
               tolerance = 1e-12)
  expect_equal(unname(result$parameter), 1)
  # A chi-square with 1 df is the square of a standard normal variable.
  expect_equal(unname(result$p.value),
               2 * stats::pnorm(-sqrt(expected.statistic)), tolerance = 1e-12)
  expect_equal(result$expected, matrix(c(5, 15, 5, 15), 2, byrow = TRUE))
  expect_identical(result$data.name, "tab")
  expect_match(paste(capture.output(print(result)), collapse = " "),
               "Likelihood-ratio G-squared")
})

test_that("RxC likelihood ratio agrees with an independent Poisson model", {
  tab <- matrix(c(109, 68, 14, 72, 79, 37, 9, 42, 139), nrow = 3,
                byrow = TRUE,
                dimnames = list(size = c("S", "M", "L"),
                                risk = c("low", "medium", "high")))
  cells <- as.data.frame(as.table(tab))
  fit <- stats::glm(Freq ~ size + risk, family = stats::poisson(), data = cells)
  result <- indep.test(tab, method = "lr")
  expect_equal(unname(result$statistic), stats::deviance(fit),
               tolerance = 1e-10)
  expect_equal(unname(result$parameter), stats::df.residual(fit))
  # For 4 df the chi-square survival function is exp(-x/2) * (1 + x/2).
  # Retain tiny upper-tail probabilities rather than rounding them to zero.
  g2 <- unname(result$statistic)
  expect_gt(result$p.value, 0)
  expect_equal(unname(result$p.value), exp(-g2 / 2) * (1 + g2 / 2),
               tolerance = 1e-12)
  expect_equal(result$expected, matrix(stats::fitted(fit), nrow = 3,
                                      dimnames = dimnames(tab)),
               tolerance = 1e-10)
  expect_identical(result$observed, tab)
  for (method in c("pearson", "lr")) {
    original <- indep.test(tab, method = method)
    transposed <- indep.test(t(tab), method = method)
    permuted <- indep.test(tab[c(3, 1, 2), c(2, 3, 1)], method = method)
    expect_equal(transposed$statistic, original$statistic)
    expect_equal(permuted$p.value, original$p.value)
  }
})

test_that("breast screening matches the validated book comparison", {
  tab <- matrix(c(1366, 383, 14205, 4983), nrow = 2)
  result <- indep.test(tab, method = "lr")
  # Also checked against DescTools 0.99.60 GTest(correct = "none").
  expect_equal(unname(result$statistic), 14.3788943329454, tolerance = 1e-10)
  expect_equal(unname(result$p.value), 0.00014946825536589, tolerance = 1e-10)
})

test_that("independence gives zero and sparse tables warn for either method", {
  independent <- matrix(c(10, 20, 30, 20, 40, 60), 2, byrow = TRUE)
  sparse <- matrix(c(0, 2, 3, 0), 2)
  for (method in c("pearson", "lr")) {
    result <- indep.test(independent, method = method)
    expect_equal(unname(result$statistic), 0)
    expect_equal(unname(result$p.value), 1)
    expect_warning(result <- indep.test(sparse, method = method),
                   "Chi-squared approximation may be incorrect")
    expect_true(is.finite(result$statistic))
    expect_true(result$p.value >= 0 && result$p.value <= 1)
  }
})

test_that("independence tests reject invalid counts, margins, and options", {
  tab <- matrix(c(12, 8, 5, 15), 2)
  invalid <- list(1:4, matrix(1:4, nrow = 1), array(1, c(2, 2, 2)),
                  matrix("a", 2, 2), matrix(TRUE, 2, 2))
  for (bad in invalid) {
    expect_error(indep.test(bad), "two-dimensional table")
  }
  for (value in c(-1, 0.5, NA, NaN, Inf)) {
    bad <- tab
    bad[1, 1] <- value
    expect_error(indep.test(bad), "non-negative integer counts")
  }
  for (bad in list(matrix(0, 2, 2), rbind(c(0, 0), c(1, 2)),
                   cbind(c(0, 0), c(1, 2)), matrix(1e308, 2, 2))) {
    expect_error(indep.test(bad), "finite total and positive row and column")
  }
  for (correct in list(NA, 1, "yes", logical(), c(TRUE, FALSE))) {
    expect_error(indep.test(tab, correct = correct), "must be TRUE or FALSE")
  }
  expect_error(indep.test(tab, method = "fisher"), "arg")
  expect_error(indep.test(tab, method = "lr", correct = TRUE),
               "only available for Pearson")
  expect_error(indep.test(matrix(10, 2, 3), correct = TRUE),
               "only available for Pearson")
})
