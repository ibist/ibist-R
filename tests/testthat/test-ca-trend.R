test_that("Cochran-Armitage test matches its signed statistic", {
  tab <- matrix(c(8, 12, 12, 8, 16, 4), nrow = 3, byrow = TRUE)
  scores <- c(1, 2, 4)
  pbar <- sum(tab[, 1]) / sum(tab)
  xbar <- sum(rowSums(tab) * scores) / sum(tab)
  variance <- pbar * (1 - pbar) *
    sum(rowSums(tab) * (scores - xbar)^2)
  expected <- sum(scores * (tab[, 1] - rowSums(tab) * pbar)) /
    sqrt(variance)

  result <- ca.trend.test(tab, scores = scores)

  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), expected)
  expect_gt(unname(result$statistic), 0)
  expect_equal(result$p.value, 2 * stats::pnorm(-abs(expected)))
})

test_that("exact trend test agrees with direct conditional enumeration", {
  tab <- matrix(c(2, 0, 1, 1, 0, 2), nrow = 3, byrow = TRUE)
  scores <- c(0, 1, 3)
  observed <- sum(scores * tab[, 1])
  allocations <- expand.grid(rep(list(0:2), 3))
  allocations <- allocations[rowSums(allocations) == sum(tab[, 1]), ]
  row.totals <- rowSums(tab)
  event.total <- sum(tab[, 1])
  weights <- apply(allocations, 1, function(k) {
    prod(choose(row.totals, k)) / choose(sum(row.totals), event.total)
  })
  statistic <- as.matrix(allocations) %*% scores
  lower <- sum(weights[statistic <= observed + 1e-12])
  upper <- sum(weights[statistic >= observed - 1e-12])

  result <- ca.trend.test(tab, scores = scores, exact = TRUE)

  expect_s3_class(result, "htest")
  expect_equal(result$p.value, min(1, 2 * min(lower, upper)))
  exact.p <- function(alt) {
    ca.trend.test(tab, scores = scores, alternative = alt,
                  exact = TRUE)$p.value
  }
  expect_equal(exact.p("greater"), upper)
  expect_equal(exact.p("less"), lower)
})

test_that("trend test validates tables, scores, and exact flag", {
  tab <- matrix(c(4, 1, 2, 3), nrow = 2)

  expect_error(ca.trend.test(matrix(1:6, nrow = 2)), "R-by-2 table")
  expect_error(ca.trend.test(matrix(c(1, -1, 2, 3), nrow = 2)),
               "non-negative integer counts")
  expect_error(ca.trend.test(tab, scores = c(2, 1)), "strictly increasing")
  expect_error(ca.trend.test(tab, exact = 1), "'exact' must be TRUE or FALSE")
  expect_error(ca.trend.test(matrix(c(4, 0, 2, 0), nrow = 2, byrow = TRUE)),
               "Both outcome categories")
})
