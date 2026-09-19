test_that("Zelen test is vacuous for one stratum", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)

  result <- or.homogeneity.test(list(tab))

  expect_s3_class(result, "htest")
  expect_equal(result$p.value, 1)
})

test_that("Zelen test uses probability ordering conditional on the total", {
  tabs <- list(
    matrix(c(2, 3, 1, 4), nrow = 2, byrow = TRUE),
    matrix(c(1, 2, 3, 2), nrow = 2, byrow = TRUE)
  )
  result <- or.homogeneity.test(tabs)

  observed_weight <- dhyper(2, 5, 5, 3) * dhyper(1, 3, 5, 4)
  reference_weights <- c(
    dhyper(0, 5, 5, 3) * dhyper(3, 3, 5, 4),
    dhyper(1, 5, 5, 3) * dhyper(2, 3, 5, 4),
    observed_weight,
    dhyper(3, 5, 5, 3) * dhyper(0, 3, 5, 4)
  )
  expected <- sum(reference_weights[reference_weights <= observed_weight]) /
    sum(reference_weights)

  expect_equal(result$p.value, expected, tolerance = 1e-12)
})

test_that("Zelen test accepts a 2 x 2 x K array", {
  tabs <- list(
    matrix(c(2, 3, 1, 4), nrow = 2, byrow = TRUE),
    matrix(c(1, 2, 3, 2), nrow = 2, byrow = TRUE)
  )
  tables <- array(unlist(tabs), dim = c(2, 2, 2))

  expect_equal(
    or.homogeneity.test(tables)$p.value,
    or.homogeneity.test(tabs)$p.value,
    tolerance = 1e-12
  )
})

test_that("Zelen test reproduces a published StatXact example", {
  # Table 10.5 in Applied Nonparametric Statistical Methods, Third Edition.
  tabs <- list(
    matrix(c(8, 1, 11, 1), nrow = 2, byrow = TRUE),
    matrix(c(14, 2, 18, 0), nrow = 2, byrow = TRUE),
    matrix(c(25, 3, 42, 2), nrow = 2, byrow = TRUE),
    matrix(c(39, 3, 22, 6), nrow = 2, byrow = TRUE)
  )

  result <- or.homogeneity.test(tabs)

  expect_equal(round(result$p.value, 4), 0.0689)
})
