test_that("single and multiple measures preserve native result structures", {
  tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), 3, byrow = TRUE)
  single <- assoc.test(tab)
  expect_s3_class(single, "htest")
  expect_identical(single$measure, "pearson")
  expect_identical(single$test.method, "multinomial")
  expect_named(single$estimate, "pearson")
  expect_equal(unname(single$null.value), 0)
  expect_identical(single$data.name, "tab")
  names <- c("tau-b", "gamma", "tau-b", "spearman")
  multiple <- assoc.test(tab, measure = names)
  expect_s3_class(multiple, "ibist_assoc_tests")
  expect_identical(names(multiple), names)
  expect_true(all(vapply(multiple, inherits, logical(1), "htest")))
  expect_equal(multiple[[1]], multiple[[3]])
  expect_identical(attr(multiple, "shared.tests")[[1]], names[1:3])
  combined <- summary(multiple)
  expect_identical(combined$measure, names)
  expect_equal(as.data.frame(multiple), combined)
  before <- multiple
  printed <- capture.output(returned <- print(multiple, digits = 2))
  expect_identical(returned, before)
  expect_identical(multiple, before)
  expect_true(any(grepl("Shared test:", printed)))
  expect_true(any(grepl("p-values are unadjusted", printed)))
})

test_that("multinomial results preserve existing ordinal estimates and tails", {
  tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), 3, byrow = TRUE)
  measures <- c("gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c",
                "pearson", "spearman")
  for (alternative in c("two.sided", "greater", "less")) {
    old <- ordinal.assoc(tab, alternative = alternative, conf.level = 0.9)
    new <- assoc.test(tab, measures, alternative = alternative,
                      conf.level = 0.9)
    for (i in seq_along(measures)) {
      expect_equal(unname(new[[i]]$estimate), unname(old$estimate[i]))
      expect_equal(unname(new[[i]]$statistic), unname(old$statistic[i]))
      expect_equal(new[[i]]$p.value, unname(old$p.value[i]))
      expect_equal(new[[i]]$std.error, unname(old$std.error[i]))
      expect_equal(as.numeric(new[[i]]$conf.int), unname(old$conf.int[i, ]))
      expect_identical(attr(new[[i]]$conf.int, "conf.level"), 0.9)
      expect_identical(new[[i]]$alternative, alternative)
    }
  }
})

test_that("nominal coefficients share the actual independence test", {
  tab <- matrix(c(40, 10, 15, 35), 2, byrow = TRUE)
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  for (correct in c(FALSE, TRUE)) {
    result <- assoc.test(tab, c("phi", "contingency", "cramer-v"),
                         method = "pearson", correct = correct)
    difference <- abs(tab - expected)
    if (correct) difference <- pmax(0, difference - 0.5)
    statistic <- sum(difference^2 / expected)
    estimate <- sqrt(sum((tab - expected)^2 / expected) / sum(tab))
    for (fit in result) {
      expect_equal(unname(fit$statistic), statistic)
      expect_equal(fit$p.value, stats::pchisq(statistic, 1, lower.tail = FALSE))
      expect_equal(unname(fit$parameter), 1)
      expect_null(fit$conf.int)
    }
    expect_equal(unname(result[[1]]$estimate), estimate)
    expect_equal(unname(result[[3]]$estimate), estimate)
  }
})

test_that("CMH inference uses score correlation and has no invented interval", {
  tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), 3, byrow = TRUE)
  for (alternative in c("two.sided", "greater", "less")) {
    for (name in c("pearson", "spearman")) {
      fit <- assoc.test(tab, name, method = "cmh", alternative = alternative)
      rho <- ordinal.assoc(tab, measure = name)$estimate
      z <- unname(rho) * sqrt(sum(tab) - 1)
      expect_equal(unname(fit$statistic), z)
      expect_equal(unname(fit$chisq), z^2)
      expect_equal(fit$p.value, switch(alternative,
        two.sided = 2 * stats::pnorm(-abs(z)),
        greater = stats::pnorm(z, lower.tail = FALSE),
        less = stats::pnorm(z)))
      expect_null(fit$conf.int)
      expect_null(fit$std.error)
      expect_null(fit$parameter)
    }
  }
  scored <- assoc.test(tab, c("pearson", "spearman"), method = "cmh",
                       row.scores = c(1, 2, 10), col.scores = c(2, 3, 9))
  old <- cmh.ll.test(tab, row.scores = c(1, 2, 10),
                     col.scores = c(2, 3, 9))
  expect_equal(scored[[1]]$statistic, old$statistic)
  expect_equal(scored[[2]]$statistic,
               cmh.ll.test(tab, measure = "spearman")$statistic)
  expect_length(attr(scored, "shared.tests"), 0)
})

test_that("invalid combinations and inapplicable options fail clearly", {
  tab <- matrix(10, 3, 3)
  invalid <- list(
    list(measure = character()), list(measure = NULL),
    list(measure = c("gamma", NA_character_)),
    list(measure = "gamma", method = "cmh"),
    list(measure = "cramer-v"), list(measure = c("gamma", "phi")),
    list(measure = "spearman", method = "pearson"),
    list(method = "invalid"), list(method = NULL),
    list(alternative = "lower"), list(alternative = NULL),
    list(conf.level = TRUE), list(conf.level = 1), list(correct = NA),
    list(correct = TRUE),
    list(measure = "phi", method = "pearson", alternative = "greater"),
    list(measure = "phi", method = "pearson", row.scores = 1:3),
    list(method = "cmh", row.scores = 1:3),
    list(measure = "spearman", method = "cmh",
         row.scores = 1:3, col.scores = 1:3),
    list(measure = "phi", method = "pearson", correct = TRUE)
  )
  for (args in invalid) {
    expect_error(do.call(assoc.test, c(list(x = tab), args)))
  }
  empty <- rbind(c(10, 0, 2), c(0, 0, 0), c(2, 0, 10))
  expect_error(assoc.test(empty, "phi", method = "pearson"), "marginal|total")
  expect_s3_class(assoc.test(empty), "htest")
  expect_s3_class(assoc.test(empty, method = "cmh"), "htest")
})
