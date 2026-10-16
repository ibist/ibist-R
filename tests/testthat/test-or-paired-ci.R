test_that("matched-pair Wald intervals reproduce the smoking example", {
  tab <- matrix(c(51, 30, 10, 41), 2, byrow = TRUE)
  result <- or.paired.ci(tab)
  expect_s3_class(result, "ibist_ci")
  expect_equal(result$estimate, c("matched-pair odds ratio" = 3))
  expect_equal(as.numeric(result$conf.int), c(1.466582, 6.136717),
               tolerance = 1e-6)
  expect_equal(result$method, "wald")
  expect_equal(colnames(result$conf.int), c("lower", "upper"))
  expect_equal(attr(result$conf.int, "conf.level"), 0.95)
})

test_that("paired Wald intervals equal pair-stratum Mantel-Haenszel limits", {
  patterns <- list(rbind(c(1, 1), c(0, 0)), diag(2)[2:1, ],
                   diag(2), rbind(c(0, 0), c(1, 1)))
  for (counts in list(c(51, 10, 30, 41), c(0, 7, 2, 0), c(5, 3, 3, 8))) {
    tab <- matrix(counts, 2)
    strata <- simplify2array(rep(patterns, counts))
    reference <- stats::mantelhaen.test(strata)
    result <- or.paired.ci(tab)
    expect_equal(unname(result$estimate), unname(reference$estimate))
    expect_equal(as.numeric(result$conf.int), as.numeric(reference$conf.int))
    expect_equal(stats::mcnemar.test(tab)$p.value, reference$p.value)
  }
})

test_that("transformed Wilson and exact intervals match base R", {
  for (conf.level in c(0.8, 0.95, 0.99)) {
    for (b in c(0, 1, 7, 30, 40)) {
      tab <- matrix(c(2, b, 40 - b, 3), 2, byrow = TRUE)
      result <- or.paired.ci(tab, conf.level,
                             method = c("wilson", "clopper-pearson"))
      wilson <- suppressWarnings(
        stats::prop.test(b, 40, conf.level = conf.level,
                        correct = FALSE)$conf.int
      )
      exact <- stats::binom.test(b, 40, conf.level = conf.level)$conf.int
      expect_equal(as.numeric(result$conf.int["wilson", ]),
                   as.numeric(wilson / (1 - wilson)), tolerance = 1e-10)
      expect_equal(as.numeric(result$conf.int["clopper-pearson", ]),
                   as.numeric(exact / (1 - exact)), tolerance = 1e-10)
    }
  }
})

test_that("mid-p intervals solve their defining tail equations", {
  for (b in c(0, 1, 7, 30, 40)) {
    probability <- paired_or_binomial_ci(b, 40 - b, 0.95, "midp")
    if (b > 0) {
      tail <- stats::pbinom(b - 1, 40, probability[1], lower.tail = FALSE)
      expect_equal(tail - stats::dbinom(b, 40, probability[1]) / 2,
                   0.025, tolerance = 1e-9)
    }
    if (b < 40) {
      expect_equal(stats::pbinom(b, 40, probability[2]) -
                     stats::dbinom(b, 40, probability[2]) / 2, 0.025,
                   tolerance = 1e-9)
    }
  }
})

test_that("Blaker intervals match independently calculated limits", {
  ## exactci::binom.exact(30, 40, tsmethod = "blaker"), with tol = 1e-10.
  expected <- c(0.589598027709465, 0.871025929500122)
  expect_equal(paired_or_binomial_ci(30, 10, 0.95, "blaker"), expected,
               tolerance = 1e-9)

  ## Check exact coverage directly, without requiring a reference package.
  probabilities <- seq(0.01, 0.99, by = 0.01)
  for (n in 1:12) {
    limits <- vapply(0:n, function(b) {
      paired_or_binomial_ci(b, n - b, 0.95, "blaker")
    }, numeric(2))
    coverage <- vapply(probabilities, function(p) {
      covered <- p >= limits[1, ] & p <= limits[2, ]
      sum(stats::dbinom(0:n, n, p)[covered])
    }, numeric(1))
    expect_gte(min(coverage), 0.95 - 1e-10)
  }
})

test_that("matched-pair intervals ignore concordance and respect direction", {
  methods <- c("wald", "wilson", "clopper-pearson", "midp", "blaker")
  tab <- matrix(c(0, 30, 10, 0), 2, byrow = TRUE)
  result <- or.paired.ci(tab, method = methods)
  tab[cbind(1:2, 1:2)] <- c(100, 200)
  expect_equal(result$conf.int, or.paired.ci(tab, method = methods)$conf.int)
  reversed <- or.paired.ci(t(tab), method = methods)
  expect_equal(reversed$conf.int[, 1], 1 / result$conf.int[, 2])
  expect_equal(reversed$conf.int[, 2], 1 / result$conf.int[, 1])
  expect_equal(unname(reversed$estimate), 1 / unname(result$estimate))
})

test_that("zero discordant counts are handled explicitly", {
  methods <- c("wilson", "clopper-pearson", "midp", "blaker")
  tab <- matrix(c(5, 4, 0, 6), 2, byrow = TRUE)
  result <- or.paired.ci(tab, method = methods)
  expect_equal(unname(result$estimate), Inf)
  expect_true(all(is.infinite(result$conf.int[, 2])))
  expect_true(all(is.finite(result$conf.int[, 1])))
  expect_error(or.paired.ci(tab), "discordant counts")
  tab[1, 2] <- 0
  expect_warning(result <- or.paired.ci(tab, method = methods), "unidentified")
  expect_true(is.na(result$estimate))
  expect_equal(as.numeric(result$conf.int[, 1]), rep(0, length(methods)))
  expect_true(all(is.infinite(result$conf.int[, 2])))
})

test_that("matched-pair intervals validate their inputs", {
  tab <- matrix(c(5, 4, 2, 6), 2, byrow = TRUE)
  expect_error(or.paired.ci(matrix(1:6, 2)), "2 x 2")
  expect_error(or.paired.ci(tab + 0.5), "integer")
  expect_error(or.paired.ci(-tab), "non-negative")
  tab[1, 1] <- NA
  expect_error(or.paired.ci(tab), "integer")
  tab[1, 1] <- 5
  for (level in list(NA, Inf, 0, 1, c(0.9, 0.95))) {
    expect_error(or.paired.ci(tab, conf.level = level), "conf.level")
  }
  expect_error(or.paired.ci(tab, method = "fisher"), "arg")
})
