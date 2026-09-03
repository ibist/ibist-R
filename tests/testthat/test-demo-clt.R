test_that("demo_clt returns a ggplot with expected simulated data", {
  set.seed(1)
  result <- demo_clt(
    runif,
    n = c(2, 4),
    nrep = 20,
    min = 0,
    max = 1,
    pmean = 0.5,
    psd = sqrt(1 / 12)
  )

  expect_s3_class(result, "ggplot")
  expect_equal(nrow(result$data), 40L)
  expect_equal(sort(unique(result$data$SampleSize)), c(2, 4))
  expect_equal(result$labels$x, "Standardized sample mean")
  expect_equal(length(result$layers), 2L)
})

test_that("demo_clt validates inputs", {
  expect_error(demo_clt("runif", n = 2), "must be a function")
  expect_error(demo_clt(runif, n = 0), "positive values")
  expect_error(demo_clt(runif, n = 2, nrep = 0), "positive integer")
})
