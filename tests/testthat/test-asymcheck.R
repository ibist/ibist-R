test_that("independence pre, attached, and post checks agree", {
  tab <- matrix(c(12, 8, 5, 15), 2, byrow = TRUE)
  pre <- asymcheck(indep.test, x = tab)
  fit <- indep.test(tab, asymcheck = TRUE)
  expect_equal(pre, fit$diagnostics)
  expect_equal(pre, asymcheck(fit))
  expect_equal(pre, asymcheck(tab, method = "pearson"))
  expect_equal(pre$details$expected, outer(c(20, 20), c(17, 23)) / 40)
  expect_true(pre$condition.met)
  expect_identical(indep.test(tab), stats::chisq.test(tab, correct = FALSE))
  expect_null(indep.test(tab)$diagnostics)
  expect_false(asymcheck(indep.test, x = tab, method = "fisher")$applicable)
})

test_that("screening criteria include their exact boundaries", {
  tab <- rbind(c(1, 5, 5, 5, 5), c(1, 5, 5, 5, 5))
  check <- asymcheck(tab, method = "pearson")
  expect_equal(check$checks$value, c(1, 0.8))
  expect_true(check$condition.met)
  expect_false(asymcheck(matrix(c(0, 1, 1, 20), 2),
                         method = "pearson")$condition.met)
})

test_that("MH matches the legacy example and direct hypergeometric variance", {
  tab <- array(c(88, 469, 48, 653, 27, 480, 34, 1324), c(2, 2, 2))
  pre <- asymcheck(tab, method = "mantel-haenszel")
  legacy <- mantelhaen.asymcheck(tab)
  expect_equal(pre$details$variance, legacy$variance)
  expect_equal(pre$details$variance.sum, 41.6369241, tolerance = 1e-8)
  simple <- matrix(1, 2, 2)
  expect_equal(asymcheck(list(simple), method = "mantel-haenszel")$
                 details$variance.sum, 1 / 3)
  expect_equal(pre, asymcheck(stats::mantelhaen.test, x = tab))
  expect_false(asymcheck(stats::mantelhaen.test, x = tab,
                         exact = TRUE)$applicable)
  expect_false(asymcheck(list(matrix(0, 2, 2)),
                         method = "mantel-haenszel")$condition.met)
})

test_that("rate checks preserve settings without arbitrary cutoffs", {
  methods <- c("wald", "exact", "wald", "log")
  pre <- asymcheck(rate.1s.ci, x = 12, exposure = 100, method = methods,
                   side = "lower", correct = FALSE)
  fit <- rate.1s.ci(12, 100, method = methods, side = "lower",
                     correct = FALSE, asymcheck = TRUE)
  expect_equal(pre, fit$diagnostics)
  expect_equal(pre, asymcheck(fit))
  expect_identical(pre$method, methods)
  expect_equal(pre$checks$value[c(2, 6, 8)], rep(1 / sqrt(12), 3))
  expect_true(is.na(pre$condition.met))
  expect_false(asymcheck(rate.1s.ci, x = 0, method = "log")$applicable)
  expect_false(asymcheck(rate.1s.ci, x = 0, method = "wald",
                         correct = FALSE)$applicable)
  expect_true(asymcheck(rate.1s.ci, x = 0, method = "wald")$applicable)
  local_mocked_bindings(ci_wald = function(...) stop("inference ran"))
  expect_s3_class(asymcheck(rate.1s.ci, x = 12, method = "wald"),
                   "ibist_asymcheck")
})

test_that("unsupported inputs and overrides fail clearly", {
  expect_error(asymcheck(mean, x = 1), "unsupported")
  expect_error(asymcheck(matrix(1, 2, 2)), "explicit method")
  expect_error(asymcheck(indep.test, x = matrix(0, 2, 2)), "margins")
  expect_error(asymcheck(rate.1s.ci, x = 2, exposure = -1), "positive")
  expect_error(asymcheck(stats::mantelhaen.test,
                         x = list(matrix(1, 2, 2))), "count array")
  expect_error(asymcheck(rate.1s.ci(2), method = "wald"), "post-checks")
  expect_error(indep.test(matrix(1, 2, 2), asymcheck = 1), "TRUE or FALSE")
})
