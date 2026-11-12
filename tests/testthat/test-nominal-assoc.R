test_that("nominal measures agree with an independent 2 x 2 formula", {
  tab <- matrix(c(30, 10, 20, 40), nrow = 2, byrow = TRUE)
  result <- nominal.assoc(tab)
  # Signed binary correlation from the cross products, not chisq.test.
  phi <- (30 * 40 - 10 * 20) / sqrt(40 * 60 * 50 * 50)
  expected <- c(phi = abs(phi), contingency = abs(phi) / sqrt(1 + phi^2),
                `cramer-v` = abs(phi))
  expect_s3_class(result, "ibist_nominal_assoc")
  expect_equal(result$estimate, expected, tolerance = 1e-12)
  expect_equal(unname(result$statistic), sum(tab) * phi^2)
  expect_identical(result$observed, tab)
  expect_equal(result$n, 100)
  expect_identical(result$data.name, "tab")
  expect_equal(nominal.assoc(tab[2:1, ])$estimate, expected)
})

test_that("RxC formulas allow Phi above one and normalize Cramer's V", {
  tab <- diag(c(10, 20, 30))
  result <- nominal.assoc(tab)
  expect_equal(result$estimate,
               c(phi = sqrt(2), contingency = sqrt(2 / 3), `cramer-v` = 1))
  rectangular <- cbind(c(10, 0), c(10, 0), c(0, 20))
  expect_equal(nominal.assoc(rectangular)$estimate,
               c(phi = 1, contingency = sqrt(1 / 2), `cramer-v` = 1))
})

test_that("book tumor estimates match the DescTools reference", {
  tab <- matrix(c(109, 68, 14, 72, 79, 37, 9, 42, 139), nrow = 3,
                byrow = TRUE)
  # Independently compared with DescTools 0.99.60 default point estimates.
  reference <- c(phi = 0.63761827201204579,
                 contingency = 0.53762824441746093,
                 `cramer-v` = 0.45086420394816618)
  expect_equal(nominal.assoc(tab)$estimate, reference, tolerance = 1e-12)
})

test_that("nominal estimates handle independent and sparse tables", {
  tab <- matrix(c(10, 20, 30, 20, 40, 60), nrow = 2, byrow = TRUE)
  expect_equal(unname(nominal.assoc(tab)$estimate), rep(0, 3))
  sparse <- matrix(c(0, 2, 3, 0), nrow = 2)
  expect_no_warning(result <- nominal.assoc(sparse))
  expect_equal(result$estimate[["phi"]], 1)
  expect_warning(indep.test(sparse), "Chi-squared approximation")
})

test_that("nominal results preserve selection, labels, and symmetries", {
  tab <- matrix(c(109, 68, 14, 72, 79, 37, 9, 42, 139), nrow = 3,
                byrow = TRUE,
                dimnames = list(size = c("S", "M", "L"),
                                risk = c("low", "medium", "high")))
  result <- nominal.assoc(tab, measure = c("cramer-v", "phi", "cramer-v"))
  expect_identical(result$measure, c("cramer-v", "phi"))
  expect_identical(names(result$estimate), result$measure)
  expect_identical(dimnames(result$expected), dimnames(tab))
  expect_equal(nominal.assoc(as.table(tab))$estimate,
               nominal.assoc(tab)$estimate)
  expect_equal(nominal.assoc(t(tab))$estimate, nominal.assoc(tab)$estimate)
  expect_equal(nominal.assoc(tab[c(3, 1, 2), c(2, 3, 1)])$estimate,
               nominal.assoc(tab)$estimate)
  expect_equal(nominal.assoc(7 * tab)$estimate, nominal.assoc(tab)$estimate)
  original <- result$estimate
  expect_output(printed <- print(result, digits = 3), "Nominal association")
  expect_output(print(result), "Estimate")
  expect_identical(printed, result)
  expect_identical(result$estimate, original)
  expect_length(nominal.assoc(tab, "phi")$estimate, 1)
})

test_that("nominal measures reject invalid inputs and empty margins", {
  tab <- matrix(c(30, 10, 20, 40), nrow = 2)
  for (bad in list(1:4, matrix(1:4, nrow = 1), array(1, c(2, 2, 2)),
                   matrix("a", 2, 2), matrix(TRUE, 2, 2))) {
    expect_error(nominal.assoc(bad), "two-dimensional table")
  }
  for (value in c(-1, 0.5, NA, NaN, Inf)) {
    bad <- tab
    bad[1, 1] <- value
    expect_error(nominal.assoc(bad), "non-negative integer counts")
  }
  for (bad in list(matrix(0, 2, 2), rbind(c(0, 0), c(1, 2)),
                   cbind(c(0, 0), c(1, 2)), matrix(1e308, 2, 2))) {
    expect_error(nominal.assoc(bad), "finite total and positive row and column")
  }
  expect_error(nominal.assoc(tab, "signed-phi"), "arg")
  expect_error(nominal.assoc(tab, character()), "length")
})
