test_that("rate.1s.ci returns stable exact intervals", {
  result <- rate.1s.ci(5, T = 10)

  expect_s3_class(result, "htest")
  expect_equal(unname(result$estimate), 0.5)
  expect_equal(unname(result$conf.int), c(0.1623486, 1.1668332),
               tolerance = 1e-7)
})

test_that("rate.test handles one- and two-sample tests", {
  one_sample <- rate.test(x = 411, T = 25800, r = 0.0119,
                          correct = FALSE)
  two_sample <- rate.test(x = c(12, 5), T = c(100, 80),
                          correct = FALSE)

  expect_s3_class(one_sample, "htest")
  expect_s3_class(two_sample, "htest")
  expect_equal(unname(one_sample$estimate), 411 / 25800)
  expect_equal(two_sample$p.value, 0.2122691, tolerance = 1e-6)
})

test_that("rate functions validate inputs", {
  expect_error(rate.1s.ci(-1), "non-negative integer")
  expect_error(rate.1s.ci(1, T = 0), "positive")
  expect_error(rate.test(x = c(1, 2), T = 1), "same length")
  expect_error(rate.test(x = c(0, 0), T = c(1, 1)),
               "at least one event")
})
