test_that("power.p1s.test returns stable normal and exact results", {
  normal <- power.p1s.test(n = 50, p0 = 0.1, p1 = 0.25)
  exact <- power.p1s.test(n = 50, p0 = 0.1, p1 = 0.25, exact = TRUE)

  expect_s3_class(normal, "power.htest")
  expect_s3_class(exact, "power.htest")
  expect_equal(normal$power, 0.8625629, tolerance = 1e-7)
  expect_true(exact$achieved.sig.level <= exact$sig.level)
})

test_that("power.p2s.test returns stable unequal-allocation power", {
  result <- power.p2s.test(n = 100, p1 = 0.3, p2 = 0.5,
                           group.rate = 2)

  expect_s3_class(result, "power.htest")
  expect_equal(result$power, 0.8980878, tolerance = 1e-7)
  expect_equal(result$note, "n is number in the 1st group")
})

test_that("power functions solve sample sizes", {
  one_sample <- power.p1s.test(p0 = 0.1, p1 = 0.25, power = 0.8)
  two_sample <- power.p2s.test(p1 = 0.3, p2 = 0.5, power = 0.8,
                               group.rate = 1.5)

  expect_equal(one_sample$n, 40.29131, tolerance = 1e-6)
  expect_equal(two_sample$n, 86.04591, tolerance = 1e-6)
})

test_that("power.p2s.test agrees with stats::power.prop.test for equal groups", {
  result <- power.p2s.test(n = 100, p1 = 0.3, p2 = 0.5,
                           group.rate = 1, correct = FALSE)
  stats_result <- stats::power.prop.test(n = 100, p1 = 0.3, p2 = 0.5,
                                         sig.level = 0.05,
                                         alternative = "two.sided")

  expect_equal(result$power, stats_result$power)
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
})
