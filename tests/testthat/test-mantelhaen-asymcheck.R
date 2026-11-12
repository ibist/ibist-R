test_that("mantelhaen.asymcheck reproduces the book example", {
  strat <- array(
    c(88, 469, 48, 653, 27, 480, 34, 1324),
    dim = c(2, 2, 2)
  )

  result <- mantelhaen.asymcheck(strat)

  expect_s3_class(result, "ibist_mantelhaen_asymcheck")
  expect_equal(result$variance.sum, 41.6369241, tolerance = 1e-8)
  expect_true(result$condition.met)
  expect_output(print(result), "Condition \\(sum >= 5\\): met")
})

test_that("mantelhaen.asymcheck accepts lists and arrays", {
  tables <- list(
    matrix(c(12, 8, 5, 15), nrow = 2),
    matrix(c(4, 6, 7, 13), nrow = 2)
  )
  array <- simplify2array(tables)

  list_result <- mantelhaen.asymcheck(tables)
  array_result <- mantelhaen.asymcheck(array)

  expect_equal(list_result$variance, array_result$variance)
  expect_equal(list_result$condition.met, array_result$condition.met)
})

test_that("mantelhaen.asymcheck reports an inadequate approximation", {
  tables <- list(
    matrix(c(1, 1, 1, 1), nrow = 2),
    matrix(c(1, 0, 0, 1), nrow = 2)
  )

  result <- mantelhaen.asymcheck(tables)

  expect_false(result$condition.met)
  expect_lt(result$variance.sum, 5)
  expect_output(print(result), "not met")
})

test_that("mantelhaen.asymcheck validates its input", {
  expect_error(mantelhaen.asymcheck(matrix(1:4, nrow = 2)), "list")
  expect_error(
    mantelhaen.asymcheck(list(matrix(1:6, nrow = 2))),
    "2 x 2"
  )
  expect_error(
    mantelhaen.asymcheck(list(matrix(c(1, 2, 3, -1), nrow = 2))),
    "non-negative integer"
  )
})
