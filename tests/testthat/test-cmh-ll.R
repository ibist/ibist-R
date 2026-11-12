test_that("CMH test matches the Pearson-score correlation formula", {
  wdbc <- matrix(
    c(109, 68, 14, 72, 79, 37, 9, 42, 139),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("Small", "Medium", "Large"),
      c("Low", "Moderate", "High")
    )
  )
  result <- cmh.ll.test(wdbc)
  pearson <- ordinal.assoc(wdbc, measure = "pearson")$estimate

  expect_s3_class(result, "htest")
  expect_equal(unname(result$estimate), unname(pearson))
  expect_equal(
    unname(result$statistic),
    sqrt(sum(wdbc) - 1) * unname(pearson)
  )
  expect_equal(
    unname(result$chisq),
    (sum(wdbc) - 1) * unname(pearson)^2
  )
  expect_equal(unname(result$statistic), 14.093, tolerance = 0.001)
  expect_equal(unname(result$chisq), 198.6126, tolerance = 0.001)
  expect_equal(
    result$p.value,
    2 * stats::pnorm(-abs(unname(result$statistic)))
  )
})

test_that("CMH Spearman test uses marginal midranks", {
  tab <- matrix(c(4, 1, 2, 5, 1, 3), nrow = 2, byrow = TRUE)
  result <- cmh.ll.test(tab, measure = "spearman")
  cell <- which(tab > 0, arr.ind = TRUE)
  row.values <- unlist(Map(rep, cell[, "row"], tab[cell]))
  col.values <- unlist(Map(rep, cell[, "col"], tab[cell]))
  expected <- stats::cor(rank(row.values), rank(col.values))

  expect_equal(unname(result$estimate), expected)
  expect_equal(
    unname(result$statistic),
    sqrt(sum(tab) - 1) * expected
  )
  expect_equal(
    result$p.value,
    2 * stats::pnorm(-abs(sqrt(sum(tab) - 1) * expected))
  )
})

test_that("CMH test supports custom Pearson scores and one-sided alternatives", {
  tab <- matrix(c(10, 1, 2, 8), nrow = 2, byrow = TRUE)
  result <- cmh.ll.test(
    tab,
    row.scores = c(0, 10),
    col.scores = c(1, 3),
    alternative = "greater"
  )

  expect_gt(unname(result$statistic), 0)
  expect_equal(
    result$p.value,
    stats::pnorm(unname(result$statistic), lower.tail = FALSE)
  )
  expect_match(paste(capture.output(print(result)), collapse = " "),
               "greater than 0")
})

test_that("CMH test validates table and score inputs", {
  tab <- matrix(c(4, 1, 2, 5), nrow = 2)

  expect_error(cmh.ll.test(matrix(c(1, -1, 3, 4), nrow = 2)),
               "non-negative integer counts")
  expect_error(cmh.ll.test(tab, row.scores = 1:2),
               "must be supplied together")
  expect_error(cmh.ll.test(tab, measure = "spearman", row.scores = 1:2,
                           col.scores = 1:2),
               "cannot be supplied for Spearman")
})
