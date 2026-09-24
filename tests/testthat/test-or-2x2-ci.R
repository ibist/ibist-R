test_that("or.2x2.ci returns stable Wald intervals", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)

  wald <- or.2x2.ci(tab)
  gart <- or.2x2.ci(tab, method = "gart")

  expect_s3_class(wald, "ibist_ci")
  expect_equal(wald$estimate, c("odds ratio" = 8.555556), tolerance = 1e-6)
  expect_equal(as.numeric(wald$conf.int), c(0.9904903, 73.9003),
               tolerance = 1e-6)
  expect_equal(as.numeric(gart$conf.int), c(0.9827961, 37.7486),
               tolerance = 1e-6)
})

test_that("or.2x2.ci returns the four conditional intervals", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)

  exact <- or.2x2.ci(tab, method = c("cornfield", "bp"))
  midp <- or.2x2.ci(
    tab,
    method = c("cornfield", "bp"),
    midp = TRUE
  )

  expect_s3_class(exact, "ibist_ci")
  expect_s3_class(midp, "ibist_ci")
  expect_equal(exact$method, c("cornfield-exact", "bp-exact"))
  expect_equal(midp$method, c("cornfield-midp", "bp-midp"))
  expect_equal(exact$conf.int["bp-exact", ], c(1, 195.495),
               tolerance = 0.01)
  expect_equal(midp$conf.int["bp-midp", ], c(1.3277, 98.8359),
               tolerance = 1e-3)
})

test_that("or.2x2.ci supports multiple methods", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)
  result <- or.2x2.ci(tab, method = c("wald", "bp"))

  expect_s3_class(result, "ibist_ci")
  expect_equal(result$method, c("wald", "bp-exact"))
  expect_equal(result$estimate, c("odds ratio" = 8.555556), tolerance = 1e-6)
  expect_equal(rownames(result$conf.int), c("wald", "bp-exact"))
  expect_equal(colnames(result$conf.int), c("lower", "upper"))
})

test_that("or.2x2.ci validates inputs", {
  expect_error(or.2x2.ci(matrix(1:6, nrow = 2)), "2 x 2")
  expect_error(or.2x2.ci(matrix(c(1, 2, 3, -1), nrow = 2)), "non-negative")
  expect_error(or.2x2.ci(matrix(c(1, 2, 3, 4.5), nrow = 2)), "integer")
  expect_error(
    or.2x2.ci(matrix(c(1, 2, 3, 4), nrow = 2), conf.level = 1),
    "conf.level"
  )
  expect_error(
    or.2x2.ci(matrix(c(1, 2, 3, 4), nrow = 2), midp = NA),
    "midp"
  )
  expect_error(
    or.2x2.ci(matrix(c(1, 0, 3, 4), nrow = 2)),
    "positive"
  )
})
