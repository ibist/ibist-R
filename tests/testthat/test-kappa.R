test_that("kappa inference matches the SAS PROC FREQ agreement example", {
  skin <- matrix(
    c(10, 4, 1, 0, 5, 10, 12, 2, 2, 4, 12, 5, 0, 2, 6, 13),
    nrow = 4,
    byrow = TRUE
  )

  result <- kappa.test(skin)

  expect_s3_class(result, "htest")
  expect_equal(unname(result$estimate), 0.3449, tolerance = 0.0001)
  expect_equal(result$std.error, 0.0724, tolerance = 0.0001)
  expect_equal(result$null.std.error, 0.0612, tolerance = 0.0005)
  expect_equal(unname(result$statistic), 5.6366, tolerance = 0.0001)
  expect_equal(
    as.numeric(result$conf.int),
    c(0.2030, 0.4868),
    tolerance = 0.0001
  )
  expect_lt(result$p.value, 0.0001)
  expect_equal(result$n, sum(skin))
})

test_that("linear weighted kappa inference matches the SAS example", {
  movies <- matrix(
    c(4, 0, 0, 0, 4, 0, 0, 4, 0),
    nrow = 3,
    byrow = TRUE
  )

  result <- kappa.test(movies, weights = "linear")

  expect_equal(unname(result$estimate), 0.5714, tolerance = 0.0001)
  expect_equal(result$std.error, 0.1323, tolerance = 0.0005)
  expect_equal(
    as.numeric(result$conf.int),
    c(0.3122, 0.8307),
    tolerance = 0.0001
  )
})

test_that("quadratic and custom weights use the supplied category order", {
  tab <- matrix(c(8, 2, 1, 5), nrow = 2, byrow = TRUE)
  quadratic <- kappa.test(tab, weights = "quadratic")
  custom <- kappa.test(tab, weights = diag(2))
  unweighted <- kappa.test(tab)

  expect_equal(quadratic$weights, matrix(c(1, 0, 0, 1), nrow = 2))
  expect_equal(unname(custom$estimate), unname(unweighted$estimate))
  expect_equal(custom$conf.int, unweighted$conf.int)
  expect_equal(custom$statistic, unweighted$statistic)

  original <- kappa.test(tab, weights = "linear")
  dimnames(tab) <- list(c("6", "10"), c("6", "10"))
  relabeled <- kappa.test(tab, weights = "linear")
  expect_equal(relabeled$estimate, original$estimate)
})

test_that("kappa test alternatives and input validation work", {
  tab <- matrix(c(9, 1, 2, 8), nrow = 2, byrow = TRUE)
  two.sided <- kappa.test(tab)
  greater <- kappa.test(tab, alternative = "greater")
  less <- kappa.test(tab, alternative = "less")

  expect_equal(greater$p.value, two.sided$p.value / 2)
  expect_equal(less$p.value, 1 - two.sided$p.value / 2)
  expect_error(kappa.test(matrix(1:6, nrow = 2)), "square table")
  expect_error(kappa.test(matrix(c(1, -1, 2, 3), nrow = 2)), "non-negative")
  expect_error(kappa.test(tab, weights = "other"), "'arg' should be one of")
  invalid.weights <- matrix(c(1, 1.2, 1.2, 1), nrow = 2)
  expect_error(
    kappa.test(tab, weights = invalid.weights),
    "custom 'weights'"
  )
  expect_error(kappa.test(matrix(0, 2, 2)), "positive total count")
})
