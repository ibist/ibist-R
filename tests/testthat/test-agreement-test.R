test_that("agreement.test preserves kappa inference without deprecation", {
  ratings <- matrix(c(9, 1, 2, 8), 2, byrow = TRUE)
  for (weights in c("unweighted", "linear", "quadratic")) {
    for (variance in c("fce", "delta")) {
      for (alternative in c("two.sided", "greater", "less")) {
        args <- list(x = ratings, weights = weights, conf.level = 0.9,
                     alternative = alternative, test.variance = variance)
        expect_warning(old <- do.call(cohen.kappa.test, args), NA)
        expect_warning(new <- do.call(agreement.test, args), NA)
        old$data.name <- new$data.name <- NULL
        expect_identical(new, old)
      }
    }
  }
  expect_identical(agreement.test(ratings)$data.name, "ratings")
  expect_identical(agreement.test(ratings, "linear")$weight.method, "linear")
})

test_that("agreement weights retain empty categories and category positions", {
  ratings <- matrix(c(4, 0, 1, 0, 0, 0, 2, 0, 5), 3, byrow = TRUE)
  result <- agreement.test(ratings, weights = "linear")
  expect_equal(result$weights, rbind(c(1, 0.5, 0), c(0.5, 1, 0.5),
                                   c(0, 0.5, 1)))
  expect_equal(result$n, 12)
  old <- cohen.kappa.test(ratings, weights = "linear")
  expect_equal(result$estimate, old$estimate)
  dimnames(ratings) <- list(c("1", "2", "10"), c("1", "2", "10"))
  expect_equal(agreement.test(ratings, weights = "linear")$estimate,
               result$estimate)
  custom <- diag(3)
  expect_identical(agreement.test(ratings, weights = custom)$weight.method,
                   "custom")
})

test_that("agreement.test validates precision without changing the old API", {
  ratings <- matrix(c(9, 1, 2, 8), 2, byrow = TRUE)
  expect_error(agreement.test(matrix(0, 2, 2)), "positive total")
  expect_error(agreement.test(matrix(1:6, 2)), "square")
  expect_error(agreement.test(matrix(c(1, 0, 0, 0), 2)), "undefined")
  expect_error(agreement.test(ratings, weights = matrix(1, 2, 2)), "undefined")
  expect_error(agreement.test(diag(4, 2), test.variance = "delta"),
                "selected test variance")
  expect_error(agreement.test(ratings, weights = "other"))
  expect_error(agreement.test(ratings, conf.level = 1))
  expect_error(agreement.test(ratings, test.variance = "other"))
  expect_error(agreement.test(matrix(c(2^53, 1, 1, 1), 2)), "2\\^53")
})
