test_that("power.p1s.test returns stable normal and exact results", {
  normal <- power.p1s.test(n = 50, p0 = 0.1, p1 = 0.25)
  exact <- power.p1s.test(n = 50, p0 = 0.1, p1 = 0.25, exact = TRUE)

  expect_s3_class(normal, "power.htest")
  expect_s3_class(exact, "power.htest")
  expect_equal(normal$power, 0.8625629, tolerance = 1e-7)
  expect_true(exact$achieved.sig.level <= exact$sig.level)
})

test_that("one-sample exact mid-p methods use their stated rejection rules", {
  n <- 35
  p0 <- 0.1
  p1 <- 0.25
  alpha <- 0.05
  x <- 0:n

  lower_midp <- pbinom(x - 1, n, p0) + 0.5 * dbinom(x, n, p0)
  upper_midp <- pbinom(x, n, p0, lower.tail = FALSE) +
    0.5 * dbinom(x, n, p0)
  reject <- 2 * pmin(lower_midp, upper_midp) <= alpha
  expected_midp <- sum(dbinom(x[reject], n, p1))

  midp <- power.p1s.test(
    n = n, p0 = p0, p1 = p1, power = NULL,
    alternative = "two.sided", exact = TRUE, exact.method = "midp"
  )
  expect_equal(midp$power, expected_midp)

  k_lo <- qbinom(alpha / 2, n, p0) - 1L
  k_hi <- qbinom(1 - alpha / 2, n, p0) + 1L
  expected_rand <-
    pbinom(k_lo - 1, n, p1) + 0.5 * dbinom(k_lo, n, p1) +
    pbinom(k_hi, n, p1, lower.tail = FALSE) +
    0.5 * dbinom(k_hi, n, p1)

  midp_rand <- power.p1s.test(
    n = n, p0 = p0, p1 = p1, power = NULL,
    alternative = "two.sided", exact = TRUE,
    exact.method = "midp-rand"
  )
  expect_equal(midp_rand$power, expected_rand)
})

test_that(
  "randomized mid-p tails treat upper and lower boundaries symmetrically",
  {
  upper <- power.p1s.test(
    n = 35, p0 = 0.1, p1 = 0.25, power = NULL,
    alternative = "greater", exact = TRUE,
    exact.method = "midp-rand"
  )
  lower <- power.p1s.test(
    n = 35, p0 = 0.9, p1 = 0.75, power = NULL,
    alternative = "less", exact = TRUE,
    exact.method = "midp-rand"
  )

  expect_equal(upper$power, lower$power)
  expect_equal(upper$achieved.sig.level, lower$achieved.sig.level)
  }
)

test_that("power.p2s.test returns stable unequal-allocation power", {
  result <- power.p2s.test(n = 100, p1 = 0.3, p2 = 0.5,
                           group.rate = 2)

  expect_s3_class(result, "power.htest")
  expect_equal(result$power, 0.8980878, tolerance = 1e-7)
  expect_equal(result$method, "Two-sample proportions power calculation")
  expect_equal(
    result$note,
    "n is number in the 1st group; continuity correction applied"
  )
})

test_that("power functions solve sample sizes", {
  one_sample <- power.p1s.test(p0 = 0.1, p1 = 0.25, power = 0.8)
  two_sample <- power.p2s.test(p1 = 0.3, p2 = 0.5, power = 0.8,
                               group.rate = 1.5)

  expect_equal(one_sample$n, 40.29131, tolerance = 1e-6)
  expect_equal(two_sample$n, 86.04591, tolerance = 1e-6)
})

test_that("power.p1s.test solves missing effect and significance inputs", {
  p1_result <- power.p1s.test(n = 50, p0 = 0.1, power = 0.8)
  p0_result <- power.p1s.test(n = 50, p1 = 0.25, power = 0.8)
  sig_result <- power.p1s.test(n = 50, p0 = 0.1, p1 = 0.25,
                               sig.level = NULL, power = 0.8)
  less_result <- power.p1s.test(n = 50, p0 = 0.25, power = 0.8,
                                alternative = "less")

  expect_equal(p1_result$p1, 0.2334746, tolerance = 1e-6)
  expect_equal(p0_result$p0, 0.1113088, tolerance = 1e-6)
  expect_equal(sig_result$sig.level, 0.02029299, tolerance = 1e-6)
  expect_equal(less_result$p1, 0.1117713, tolerance = 1e-6)
})

test_that("power.p2s.test agrees with stats::power.prop.test for equal groups", {
  result <- power.p2s.test(n = 100, p1 = 0.3, p2 = 0.5,
                           group.rate = 1, correct = FALSE)
  stats_result <- stats::power.prop.test(n = 100, p1 = 0.3, p2 = 0.5,
                                         sig.level = 0.05,
                                         alternative = "two.sided")

  expect_equal(result$power, stats_result$power)
  expect_equal(result$method, "Two-sample proportions power calculation")
  expect_equal(
    result$note,
    "n is number in the 1st group; no continuity correction"
  )
})

test_that("power functions validate inputs", {
  expect_error(power.p1s.test(n = 10, p0 = 0.1, p1 = 0.2,
                              power = 0.8),
               "Exactly one")
  expect_error(power.p1s.test(n = NULL, p0 = 0.1, p1 = 0.2,
                              power = 0.8, max_n = 0),
               "max_n")
  expect_error(power.p2s.test(n = 10, p1 = -0.1, p2 = 0.2),
               "p1")
  expect_error(power.p2s.test(n = 10, p1 = 0.1, p2 = 0.2,
                              group.rate = 0),
               "group.rate")
  expect_error(power.p2s.test(n = 10, p1 = 0.1, p2 = 1),
               "p2")
  expect_error(power.p1s.test(n = c(10, 20), p0 = 0.1, p1 = 0.2,
                              power = NULL),
               "n")
  expect_error(power.p1s.test(n = 10, p0 = 0.1, p1 = 0.2,
                              power = NULL, correct = NA),
               "correct")
  expect_error(power.p2s.test(n = 10, p1 = 0.1, p2 = 0.2,
                              correct = c(TRUE, FALSE)),
               "correct")
  expect_error(power.p2s.test(n = 10, p1 = 0.1, p2 = 0.2,
                              tol = c(1e-4, 1e-5)),
               "tol")
})
