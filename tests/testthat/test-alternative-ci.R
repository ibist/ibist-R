test_that("kappa interval directions use the estimation variance", {
  tab <- matrix(c(9, 1, 2, 8), nrow = 2, byrow = TRUE)
  for (level in c(0.4, 0.8, 0.95, 0.99)) {
    for (variance in c("fce", "delta")) {
      for (weights in c("unweighted", "linear", "quadratic")) {
        two <- cohen.kappa.test(tab, weights = weights, conf.level = level,
                                test.variance = variance)
        greater <- cohen.kappa.test(tab, weights = weights, conf.level = level,
                                    alternative = "greater",
                                    test.variance = variance)
        less <- cohen.kappa.test(tab, weights = weights, conf.level = level,
                                 alternative = "less", test.variance = variance)
        estimate <- unname(two$estimate)
        half <- qnorm(level) * two$std.error
        expect_equal(as.numeric(greater$conf.int),
                     c(max(-1, min(1, estimate - half)), 1))
        expect_equal(as.numeric(less$conf.int),
                     c(-1, min(1, max(-1, estimate + half))))
        expect_equal(attr(greater$conf.int, "conf.level"), level)
        expect_equal(attr(less$conf.int, "conf.level"), level)
        expect_equal(greater$p.value + less$p.value, 1)
        expect_equal(greater$std.error, two$std.error)
        expect_equal(greater$statistic, two$statistic)
      }
    }
  }
})

test_that("all ordinal measures have direction-sensitive intervals", {
  tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), nrow = 3, byrow = TRUE)
  for (level in c(0.4, 0.8, 0.95, 0.99)) {
    two <- ordinal.assoc(tab, conf.level = level)
    greater <- ordinal.assoc(tab, conf.level = level, alternative = "greater")
    less <- ordinal.assoc(tab, conf.level = level, alternative = "less")
    half <- qnorm(level) * two$std.error
    expect_equal(unname(greater$conf.int[, "lower"]),
                 pmax(-1, pmin(1, unname(two$estimate - half))))
    expect_equal(unname(less$conf.int[, "upper"]),
                 pmin(1, pmax(-1, unname(two$estimate + half))))
    expect_true(all(greater$conf.int[, "upper"] == 1))
    expect_true(all(less$conf.int[, "lower"] == -1))
    expect_equal(rownames(greater$conf.int), two$measure)
    expect_equal(greater$p.value + less$p.value, rep(1, 7), ignore_attr = TRUE)
    expect_equal(greater$statistic, two$statistic)
  }
})

test_that("one-sided rate score limits invert the directional score test", {
  for (count in c(0, 1, 5, 411)) {
    for (correct in c(FALSE, TRUE)) {
      for (level in c(0.4, 0.8, 0.95, 0.99)) {
        greater <- rate.test(count, exposure = 10, alternative = "greater",
                              conf.level = level, correct = correct)
        less <- rate.test(count, exposure = 10, alternative = "less",
                           conf.level = level, correct = correct)
        z <- qnorm(level)
        lower.count <- max(0, count - if (correct) 0.5 else 0)
        upper.count <- count + if (correct) 0.5 else 0
        lower <- ((sqrt(z ^ 2 + 4 * lower.count) - z) / 2) ^ 2 / 10
        upper <- ((sqrt(z ^ 2 + 4 * upper.count) + z) / 2) ^ 2 / 10
        expect_equal(as.numeric(greater$conf.int), c(lower, Inf))
        expect_equal(as.numeric(less$conf.int), c(0, upper))
        expect_equal(attr(greater$conf.int, "conf.level"), level)
        expect_named(greater$conf.int, c("lower", "upper"))
        if (!correct) expect_equal(greater$p.value + less$p.value, 1)
      }
    }
  }
})

test_that("two-sample rate bounds use directional log-Wald limits", {
  counts <- c(151, 55)
  exposure <- c(57518.1, 74573.5)
  estimate <- (counts[1] / exposure[1]) / (counts[2] / exposure[2])
  se <- sqrt(sum(1 / counts))
  for (level in c(0.4, 0.8, 0.95, 0.99)) {
    greater <- rate.test(counts, exposure, alternative = "greater",
                          conf.level = level, correct = FALSE)
    less <- rate.test(counts, exposure, alternative = "less",
                       conf.level = level, correct = FALSE)
    expect_equal(as.numeric(greater$conf.int),
                 c(exp(log(estimate) - qnorm(level) * se), Inf))
    expect_equal(as.numeric(less$conf.int),
                 c(0, exp(log(estimate) + qnorm(level) * se)))
    expect_equal(greater$p.value + less$p.value, 1)
    expect_equal(attr(greater$conf.int, "conf.level"), level)
  }
})
