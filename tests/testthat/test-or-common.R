test_that("common odds ratio test matches Fisher one-sided tests", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)

  expect_equal(
    or.common.test(list(tab), alternative = "greater")$p.value,
    fisher.test(tab, alternative = "greater")$p.value,
    tolerance = 1e-12
  )
  expect_equal(
    or.common.test(list(tab), alternative = "less")$p.value,
    fisher.test(tab, alternative = "less")$p.value,
    tolerance = 1e-12
  )
})

test_that("common odds ratio test combines strata by direct convolution", {
  tabs <- list(
    matrix(c(2, 3, 1, 4), nrow = 2, byrow = TRUE),
    matrix(c(1, 2, 3, 2), nrow = 2, byrow = TRUE)
  )
  result <- or.common.test(tabs, alternative = "greater")

  probabilities <- numeric()
  totals <- integer()
  for (x1 in 0:3) {
    for (x2 in 0:3) {
      probabilities <- c(probabilities, NA_real_)
      totals <- c(totals, x1 + x2)
      probabilities[length(probabilities)] <-
        dhyper(x1, 5, 5, 3) * dhyper(x2, 3, 5, 4)
    }
  }
  expect_equal(result$p.value, sum(probabilities[totals >= 3]),
               tolerance = 1e-12)
})

test_that("common odds ratio test reproduces Jung's published example", {
  # Li et al. (1979), as reproduced in Jung (2014), Table 2.
  tabs <- list(
    matrix(c(10, 1, 12, 1), nrow = 2, byrow = TRUE),
    matrix(c(9, 0, 11, 1), nrow = 2, byrow = TRUE),
    matrix(c(8, 0, 7, 3), nrow = 2, byrow = TRUE)
  )

  result <- or.common.test(tabs, alternative = "greater")

  expect_equal(result$p.value, 0.1563, tolerance = 5e-5)
})

test_that("common odds ratio confidence interval inverts exact tails", {
  tabs <- list(
    matrix(c(7, 3, 2, 8), nrow = 2, byrow = TRUE),
    matrix(c(4, 6, 1, 9), nrow = 2, byrow = TRUE)
  )
  ci <- or.common.ci(tabs)
  convolution_ci <- or.common.ci(tabs, algorithm = "convolution")

  expect_s3_class(ci, "ibist_ci")
  expect_equal(ci$conf.int, convolution_ci$conf.int, tolerance = 1e-10)
  expect_true(ci$conf.int[1, "lower"] < ci$estimate)
  expect_true(ci$estimate < ci$conf.int[1, "upper"])
  expect_equal(
    or.common.test(tabs, theta = ci$conf.int[1, "lower"],
                   alternative = "greater")$p.value,
    0.025,
    tolerance = 2e-6
  )
  expect_equal(
    or.common.test(tabs, theta = ci$conf.int[1, "upper"],
                   alternative = "less")$p.value,
    0.025,
    tolerance = 2e-6
  )
})

test_that("common odds ratio functions validate input", {
  expect_error(or.common.test(matrix(1:4, nrow = 2)), "list")
  expect_error(or.common.test(list(matrix(1:6, nrow = 2))), "2 x 2")
  expect_error(or.common.test(list(matrix(1:4, nrow = 2)), theta = 0),
               "positive")
  expect_error(or.common.ci(list(matrix(1:4, nrow = 2)), conf.level = 1),
               "conf.level")
})
