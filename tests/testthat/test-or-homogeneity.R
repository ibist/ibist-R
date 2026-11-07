test_that("Zelen test is vacuous for one stratum", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)

  result <- or.homogeneity.test(list(tab))

  expect_s3_class(result, "htest")
  expect_equal(result$p.value, 1)
  expect_identical(result, or.homogeneity.test(list(tab), method = "zelen"))
})

test_that("Breslow-Day handles a common odds ratio exactly equal to one", {
  tabs <- list(matrix(c(30, 20, 20, 30), 2, byrow = TRUE),
               matrix(c(20, 30, 30, 20), 2, byrow = TRUE))
  # MH odds ratio is one. Each fitted cell is 25, so each variance is 25/4.
  for (correct in c(FALSE, TRUE)) {
    result <- or.homogeneity.test(tabs, "breslow-day", correct = correct)
    expect_s3_class(result, "htest")
    expect_equal(unname(result$statistic), 8, tolerance = 1e-12)
    expect_equal(unname(result$parameter), 1)
    expect_equal(result$p.value, 2 * stats::pnorm(-sqrt(8)))
  }
  woolf <- or.homogeneity.test(tabs, "woolf")
  expect_equal(unname(woolf$statistic), 12 * log(2.25)^2, tolerance = 1e-12)
})

test_that("Woolf applies a zero-cell adjustment across all strata", {
  tabs <- list(matrix(c(0, 10, 20, 30), 2, byrow = TRUE),
               matrix(c(5, 15, 25, 35), 2, byrow = TRUE))
  # Once a zero is present, both tables receive 0.5 in every cell.
  l1 <- log(0.5 * 30.5 / (10.5 * 20.5))
  l2 <- log(5.5 * 35.5 / (15.5 * 25.5))
  v1 <- 1 / 0.5 + 1 / 10.5 + 1 / 20.5 + 1 / 30.5
  v2 <- 1 / 5.5 + 1 / 15.5 + 1 / 25.5 + 1 / 35.5
  expect_warning(result <- or.homogeneity.test(tabs, "woolf"),
                 "Chi-squared approximation")
  expect_equal(unname(result$statistic), (l1 - l2)^2 / (v1 + v2),
               tolerance = 1e-12)
})

test_that("homogeneity tests match the book's DescTools reference", {
  x <- array(c(88, 469, 48, 653, 27, 480, 34, 1324), c(2, 2, 2))
  # Independently compared with DescTools 0.99.60 default MH estimation.
  bd <- or.homogeneity.test(x, "breslow-day")
  tarone <- or.homogeneity.test(x, "breslow-day", correct = TRUE)
  woolf <- or.homogeneity.test(x, "woolf")
  expect_equal(unname(bd$statistic), 0.22363470327382157,
               tolerance = 1e-10)
  expect_equal(bd$p.value, 0.63628430370898648, tolerance = 1e-10)
  expect_equal(unname(tarone$statistic), 0.22280033792545345,
               tolerance = 1e-10)
  expect_equal(tarone$p.value, 0.63691443542242898, tolerance = 1e-10)
  expect_equal(unname(woolf$statistic), 0.22265002126181524,
               tolerance = 1e-10)
  expect_equal(woolf$p.value, 0.63702811131781434, tolerance = 1e-10)
  expect_lte(unname(tarone$statistic), unname(bd$statistic))
  expect_match(tarone$method, "Tarone")
})

test_that("asymptotic homogeneity tests preserve table symmetries", {
  x <- array(c(88, 469, 48, 653, 27, 480, 34, 1324), c(2, 2, 2))
  for (method in c("breslow-day", "woolf")) {
    result <- or.homogeneity.test(x, method)
    variants <- list(aperm(x, c(2, 1, 3)), x[2:1, , ], x[, , 2:1],
                     list(x[, , 1], x[, , 2]))
    for (variant in variants) {
      other <- or.homogeneity.test(variant, method)
      expect_equal(other$statistic, result$statistic, tolerance = 1e-10)
      expect_equal(other$p.value, result$p.value, tolerance = 1e-10)
    }
    expect_identical(result$data.name, "x")
  }
})

test_that("asymptotic methods reject undefined tests and invalid options", {
  tab <- matrix(c(10, 20, 30, 40), 2)
  for (method in c("breslow-day", "woolf")) {
    expect_error(or.homogeneity.test(list(tab), method), "at least two strata")
    empty.margin <- matrix(c(0, 0, 1, 2), 2)
    expect_error(or.homogeneity.test(list(tab, empty.margin), method),
                 "positive row and column totals")
  }
  separated <- list(matrix(c(0, 10, 20, 30), 2, byrow = TRUE),
                    matrix(c(0, 20, 30, 40), 2, byrow = TRUE))
  expect_error(or.homogeneity.test(separated, "breslow-day"),
               "finite positive Mantel-Haenszel")
  for (correct in list(NA, 1, logical(), c(TRUE, FALSE))) {
    expect_error(or.homogeneity.test(list(tab, tab), correct = correct),
                 "must be TRUE or FALSE")
  }
  expect_error(or.homogeneity.test(list(tab, tab), correct = TRUE),
               "only available for the Breslow-Day")
  expect_error(or.homogeneity.test(list(tab, tab), "woolf", correct = TRUE),
               "only available for the Breslow-Day")
  expect_error(or.homogeneity.test(list(tab, tab), "unknown"), "arg")
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
