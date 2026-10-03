test_that("ordinal.assoc computes Pearson and Spearman inference", {
  tab <- matrix(
    c(12, 3, 1, 4, 9, 2, 1, 5, 13),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("low", "medium", "high"), c("1", "2", "3"))
  )

  result <- ordinal.assoc(
    tab,
    measure = c("pearson", "spearman")
  )

  expect_named(result$estimate, c("pearson", "spearman"))
  expect_named(result$statistic, c("pearson", "spearman"))
  expect_named(result$p.value, c("pearson", "spearman"))
  expect_true(all(is.finite(result$std.error)))
  expect_true(all(is.finite(result$conf.int)))
  expect_true(all(result$p.value >= 0 & result$p.value <= 1))

  cell <- which(tab > 0, arr.ind = TRUE)
  row.values <- unlist(Map(
    rep,
    cell[, "row"],
    tab[cell]
  ))
  col.values <- unlist(Map(
    rep,
    as.numeric(colnames(tab)[cell[, "col"]]),
    tab[cell]
  ))
  expected.pearson <- stats::cor(row.values, col.values)
  expect_equal(unname(result$estimate["pearson"]), expected.pearson)
})

test_that("ordinal.assoc reports shared tests for concordance measures", {
  tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), nrow = 3, byrow = TRUE)
  result <- ordinal.assoc(tab, measure = c("gamma", "tau-b", "pearson"))
  output <- capture.output(print(result))

  expect_true(any(grepl("Concordance-based measures share the same test", output)))
})

test_that("Pearson category scores can be specified", {
  tab <- matrix(c(10, 1, 2, 8), nrow = 2, byrow = TRUE)
  default <- ordinal.assoc(tab, measure = "pearson")
  custom <- ordinal.assoc(
    tab,
    measure = "pearson",
    row.scores = c(0, 10),
    col.scores = c(1, 3)
  )

  expect_equal(custom$estimate, default$estimate)
  expect_equal(custom$std.error, default$std.error)
  expect_error(
    ordinal.assoc(tab, measure = "pearson", row.scores = c(1, 1)),
    "at least two distinct values"
  )
})

test_that("numeric category labels are used as the default Pearson scores", {
  tab <- matrix(
    c(8, 1, 4, 2, 5, 9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("1", "2"), c("1", "2", "10"))
  )
  default <- ordinal.assoc(tab, measure = "pearson")
  explicit <- ordinal.assoc(
    tab,
    measure = "pearson",
    row.scores = c(1, 2),
    col.scores = c(1, 2, 10)
  )
  table.order <- ordinal.assoc(
    tab,
    measure = "pearson",
    row.scores = c(1, 2),
    col.scores = 1:3
  )

  expect_equal(default$estimate, explicit$estimate)
  expect_false(isTRUE(all.equal(default$estimate, table.order$estimate)))
})

test_that("Spearman correlation uses marginal midranks for tied categories", {
  tab <- matrix(c(4, 1, 2, 5), nrow = 2, byrow = TRUE)
  result <- ordinal.assoc(tab, measure = "spearman")
  cell <- which(tab > 0, arr.ind = TRUE)
  row.values <- unlist(Map(rep, cell[, "row"], tab[cell]))
  col.values <- unlist(Map(rep, cell[, "col"], tab[cell]))
  expected <- stats::cor(rank(row.values), rank(col.values))

  expect_equal(unname(result$estimate), expected)
})

test_that("Pearson and Spearman inference matches the SAS PROC FREQ example", {
  # SAS/STAT Example 3.8, Clinical Trial for Treatment of Pain:
  # https://support.sas.com/documentation/cdl/en/procstat/70116/HTML/default/procstat_freq_examples08.htm
  pain <- matrix(
    c(26, 26, 23, 18, 9, 6, 7, 9, 14, 23),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("No", "Yes"), as.character(0:4))
  )

  result <- ordinal.assoc(
    pain,
    measure = c("pearson", "spearman")
  )

  expect_equal(
    round(unname(result$estimate), 4),
    c(0.3776, 0.3771)
  )
  expect_equal(
    round(unname(result$std.error), 4),
    c(0.0714, 0.0718)
  )
  expect_equal(
    round(unname(result$conf.int[, "lower"]), 4),
    c(0.2378, 0.2363)
  )
  expect_equal(
    round(unname(result$conf.int[, "upper"]), 4),
    c(0.5175, 0.5178)
  )
})
