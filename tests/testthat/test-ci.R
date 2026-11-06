test_that("printing changes presentation without rounding stored intervals", {
  result <- prop.1s.ci(5, 40)
  before <- result
  expect_output(print(result, digits = 3), "sample estimates")
  expect_output(print(result, digits = 3), "method: wilson")
  expect_output(print(result, digits = 3), "95% confidence interval")
  expect_identical(result, before)
  expect_identical(withVisible(print(result))$visible, FALSE)
  multiple <- prop.1s.ci(5, 40, method = c("wilson", "wald"))
  expect_output(print(multiple), "lower.*upper")
})
