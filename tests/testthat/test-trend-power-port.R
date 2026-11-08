test_that("trend.test preserves the existing interface and results", {
  tab <- matrix(c(2, 0, 1, 1, 0, 2), 3, byrow = TRUE)
  for (exact in c(FALSE, TRUE)) {
    expect_warning(old <- ca.trend.test(tab, exact = exact), NA)
    expect_warning(new <- trend.test(tab, exact = exact), NA)
    expect_identical(old, new)
  }
})

test_that("exact sample size is minimal and respects max_n", {
  result <- power.p1s.test(p0 = 0.01, p1 = 0.99, power = 0.8,
                            alternative = "greater", exact = TRUE,
                            size.rule = "minimal", max_n = 100)
  expect_equal(result$n, 1)
  expect_error(power.p1s.test(p0 = 0.5, p1 = 0.99, power = 0.5,
                              alternative = "greater", exact = TRUE,
                              size.rule = "minimal", max_n = 4), "max_n")
})
