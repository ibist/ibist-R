test_that("one-proportion intervals reproduce SAS published limits", {
  ## SAS/STAT 14.2 FREQ, Output 41.4.1, p. 2857: 341 brown eyes / 762.
  expected <- list(wilson = c(0.4181, 0.4773),
                   `agresti-coull` = c(0.4181, 0.4773),
                   `clopper-pearson` = c(0.4174, 0.4779))
  for (method in names(expected)) {
    result <- prop.1s.ci(341, 762, conf.level = 0.9, method = method)
    expect_equal(round(as.numeric(result$conf.int), 4), expected[[method]])
    expect_equal(unname(result$estimate), 341 / 762)
  }
})

test_that("Wilson and Clopper-Pearson agree with base R at boundaries", {
  for (n in c(1, 2, 10, 40)) {
    for (x in unique(c(0, 1, floor(n / 2), n))) {
      for (alternative in c("two.sided", "less", "greater")) {
        for (correct in c(FALSE, TRUE)) {
          ours <- prop.1s.ci(x, n, correct = correct,
                             alternative = alternative)$conf.int
          ## prop.test caps its correction using the specified null p;
          ## choose a distant null so that the full half-count is used.
          null <- if (x / n < 0.5) 0.999 else 0.001
          reference <- suppressWarnings(prop.test(x, n, p = null,
                                                  correct = correct,
                                                  alternative = alternative))
          expect_equal(as.numeric(ours), as.numeric(reference$conf.int),
                       tolerance = 1e-10)
        }
        ours <- prop.1s.ci(x, n, method = "clopper-pearson",
                           alternative = alternative)$conf.int
        reference <- binom.test(x, n, alternative = alternative)
        expect_equal(as.numeric(ours), as.numeric(reference$conf.int),
                     tolerance = 1e-10)
      }
    }
  }
})

test_that("mid-p limits invert the inclusive half-mass tails", {
  result <- prop.1s.ci(5, 40, method = "clopper-pearson", midp = TRUE)
  lower <- result$conf.int[1L, "lower"]
  upper <- result$conf.int[1L, "upper"]
  expect_equal(pbinom(5, 40, lower, lower.tail = FALSE) +
                 0.5 * dbinom(5, 40, lower), 0.025, tolerance = 1e-8)
  expect_lt(abs(pbinom(4, 40, upper) + 0.5 * dbinom(5, 40, upper) - 0.025),
            1e-8)
  expect_equal(prop.1s.ci(0, 10, method = "clopper-pearson",
                         midp = TRUE)$conf.int[1L, "upper"],
               1 - 0.05^(1 / 10), tolerance = 1e-8)
})

test_that("Blaker bounds contain the exact acceptance set", {
  for (n in 1:8) {
    support <- 0:n
    result <- prop.1s.ci(support, n, method = "blaker", tol = 1e-10)
    intervals <- result$conf.int
    opposite <- prop.1s.ci(n - support, n, method = "blaker", tol = 1e-10)
    expect_equal(intervals[, "lower"],
                 1 - opposite$conf.int[, "upper"], ignore_attr = TRUE)
    for (p in seq(0.005, 0.995, length.out = 41)) {
      mass <- dbinom(support, n, p)
      ordering <- pmin(pbinom(support, n, p),
                       pbinom(support - 1, n, p, lower.tail = FALSE))
      accepted <- vapply(ordering, function(value) {
        sum(mass[ordering <= value + 1e-14]) > 0.05
      }, logical(1))
      included <- intervals[, "lower"] <= p & intervals[, "upper"] >= p
      expect_true(all(included[accepted]))
      expect_gte(sum(mass[included]), 0.95 - 1e-12)
    }
  }
})

test_that("vector inputs and multiple methods keep input and method order", {
  result <- prop.1s.ci(c(0, 5, 40), 40, method = c("wilson", "wald"))
  expect_s3_class(result, "ibist_ci")
  expect_equal(result$inputs$input, rep(1:3, each = 2))
  expect_equal(result$inputs$x, rep(c(0, 5, 40), each = 2))
  expect_equal(rownames(result$conf.int),
               c("wilson[1]", "wald[1]", "wilson[2]", "wald[2]",
                 "wilson[3]", "wald[3]"))
  for (i in 1:3) {
    single <- prop.1s.ci(c(0, 5, 40)[i], 40, method = c("wilson", "wald"))
    expect_equal(as.numeric(result$conf.int[(2 * i - 1):(2 * i), ]),
                 as.numeric(single$conf.int))
  }
})

test_that("independent-proportion interfaces accept counts and tables", {
  tab <- matrix(c(38, 62, 267, 717), nrow = 2, byrow = TRUE)
  methods <- c("newcombe", "wald", "miettinen-nurminen", "mee", "agresti-caffo")
  result <- prop.2s.ci(c(38, 267), c(100, 984), method = methods)
  expect_equal(result, prop.2s.ci(tab, method = methods))
  expect_equal(unname(result$estimate), 38 / 100 - 267 / 984)
  expect_equal(result$method, methods)
  for (correct in c(FALSE, TRUE)) {
    ours <- prop.2s.ci(c(38, 267), c(100, 984), method = "wald",
                       correct = correct)
    reference <- prop.test(c(38, 267), c(100, 984), correct = correct)
    expect_equal(as.numeric(ours$conf.int), as.numeric(reference$conf.int))
  }
})

test_that("Newcombe combines the separate Wilson distances", {
  x <- c(1366, 383)
  n <- c(15571, 5366)
  p <- x / n
  for (correct in c(FALSE, TRUE)) {
    one <- prop.test(x[1], n[1], correct = correct)$conf.int
    two <- prop.test(x[2], n[2], correct = correct)$conf.int
    expected <- p[1] - p[2] + c(-1, 1) *
      sqrt((as.numeric(one) - p[1])^2 + (rev(as.numeric(two)) - p[2])^2)
    ours <- prop.2s.ci(x, n, correct = correct)
    expect_equal(as.numeric(ours$conf.int), expected)
  }
})

test_that("Miettinen-Nurminen and Mee limits satisfy their score equations", {
  x <- c(12, 4)
  n <- c(40, 30)
  for (method in c("miettinen-nurminen", "mee")) {
    result <- prop.2s.ci(x, n, method = method)
    for (delta in as.numeric(result$conf.int)) {
      loglik <- function(q) {
        dbinom(x[1], n[1], q + delta, log = TRUE) +
          dbinom(x[2], n[2], q, log = TRUE)
      }
      fit <- optimize(loglik, c(max(0, -delta), min(1, 1 - delta)),
                      maximum = TRUE, tol = 1e-12)
      p <- c(fit$maximum + delta, fit$maximum)
      factor <- if (method == "miettinen-nurminen") sum(n) / (sum(n) - 1) else 1
      variance <- sum(p * (1 - p) / n) * factor
      statistic <- (x[1] / n[1] - x[2] / n[2] - delta)^2 / variance
      expect_equal(statistic, qchisq(0.95, 1), tolerance = 1e-6)
    }
  }
})

test_that("Mee and Miettinen-Nurminen have distinct zero-event limits", {
  n <- c(3, 5)
  critical <- qchisq(0.95, 1)
  for (method in c("mee", "miettinen-nurminen")) {
    factor <- if (method == "miettinen-nurminen") sum(n) / (sum(n) - 1) else 1
    expected <- c(-critical * factor / (n[2] + critical * factor),
                   critical * factor / (n[1] + critical * factor))
    result <- prop.2s.ci(c(0, 0), n, method = method)
    expect_equal(as.numeric(result$conf.int), expected, tolerance = 1e-8)
  }
  result <- prop.2s.ci(c(12, 4), c(40, 30),
                       method = c("miettinen-nurminen", "mee"))
  expect_lt(result$conf.int[1, "lower"], result$conf.int[2, "lower"])
  expect_gt(result$conf.int[1, "upper"], result$conf.int[2, "upper"])
})

test_that("independent intervals reverse sign and give one-sided bounds", {
  for (x in list(c(0, 0), c(0, 5), c(3, 0), c(3, 5), c(2, 4))) {
    n <- c(3, 5)
    for (method in c("newcombe", "wald", "miettinen-nurminen", "mee",
                     "agresti-caffo")) {
      result <- prop.2s.ci(x, n, conf.level = 0.9, method = method)
      reverse <- prop.2s.ci(rev(x), rev(n), conf.level = 0.9, method = method)
      expect_equal(as.numeric(result$conf.int),
                   -rev(as.numeric(reverse$conf.int)), tolerance = 1e-7)
      lower <- prop.2s.ci(x, n, method = method, alternative = "greater")
      upper <- prop.2s.ci(x, n, method = method, alternative = "less")
      expect_equal(as.numeric(lower$conf.int), c(result$conf.int[1], 1))
      expect_equal(as.numeric(upper$conf.int), c(-1, result$conf.int[2]))
    }
  }
})

test_that("binomial interval functions reject invalid and inapplicable inputs", {
  for (x in list(NA, Inf, -1, 1.5, "1", numeric())) {
    expect_error(prop.1s.ci(x, 10), "integer counts")
  }
  expect_error(prop.1s.ci(11, 10), "exceed")
  expect_error(prop.1s.ci(1:3, 1:2), "same length")
  expect_error(prop.1s.ci(1, 0), "integer counts")
  expect_error(prop.1s.ci(1, 10, conf.level = NA), "conf.level")
  expect_error(prop.1s.ci(1, 10, tol = 0), "tol")
  expect_error(prop.1s.ci(1, 10, correct = 1), "TRUE or FALSE")
  expect_error(prop.1s.ci(1, 10, midp = NA), "TRUE or FALSE")
  expect_error(prop.1s.ci(1, 10, midp = TRUE), "applies only")
  expect_error(prop.1s.ci(1, 10, method = "blaker", correct = TRUE),
               "applies only")
  expect_error(prop.2s.ci(1, 10), "two counts")
  expect_error(prop.2s.ci(c(1, 2), c(0, 3)), "integer counts")
  expect_error(prop.2s.ci(matrix(1:4, 2), n = c(3, 7)), "Omit")
  expect_error(prop.2s.ci(c(1, 2), c(3, 3), method = "agresti-caffo",
                         correct = TRUE), "applies only")
  expect_error(prop.2s.ci(c(1, 2), c(3, 3), method = "mee",
                         correct = TRUE), "applies only")
})
