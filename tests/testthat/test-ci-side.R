ci_side_calls <- function() {
  tab <- matrix(c(7, 27, 1, 33), 2, byrow = TRUE)
  list(
    function(...) prop.1s.ci(5, 40, ...),
    function(...) prop.2s.ci(c(12, 4), c(40, 40), ...),
    function(...) rr.2x2.ci(tab, ...),
    function(...) rate.1s.ci(5, 10, ...),
    function(...) rate.2s.ci(c(12, 5), c(100, 80), ...),
    function(...) or.2x2.ci(tab, ...),
    function(...) or.paired.ci(tab, ...),
    function(...) prop.paired.ci(7, 3, 20, ...),
    function(...) or.common.ci(list(tab, tab), ...)
  )
}

test_that("every public CI entry point uses side and checks aliases", {
  for (call in ci_side_calls()) {
    expect_identical(call()$side, "two.sided")
    for (side in c("two.sided", "lower", "upper")) {
      alias <- switch(side, two.sided = "two.sided",
                      lower = "greater", upper = "less")
      result <- call(side = side)
      expect_identical(result$side, side)
      expect_equal(result$conf.int, call(alternative = alias)$conf.int)
      expect_equal(result$conf.int,
                   call(side = side, alternative = alias)$conf.int)
    }
    expect_error(call(side = "two.sided", alternative = "greater"),
                 "conflicting")
    expect_error(call(side = "lower", alternative = "less"), "conflicting")
    for (side in list("greater", "less", "two-sided", NA, TRUE,
                      c("lower", "upper"))) {
      expect_error(call(side = side), "side")
    }
  }
})

test_that("one-sided Poisson limits invert the full Poisson tail", {
  for (level in c(0.2, 0.8, 0.95, 0.99)) {
    for (x in c(0, 1, 5)) {
      lower <- rate.1s.ci(x, 10, conf.level = level, side = "lower")
      upper <- rate.1s.ci(x, 10, conf.level = level, side = "upper")
      if (x > 0) {
        expect_equal(ppois(x - 1, lower$conf.int[1] * 10,
                           lower.tail = FALSE), 1 - level)
      } else expect_equal(lower$conf.int[1], 0)
      expect_equal(ppois(x, upper$conf.int[2] * 10), 1 - level)
      expect_equal(lower$conf.int[2], Inf)
      expect_equal(upper$conf.int[1], 0)
    }
  }
})

test_that("rate ratio exact bounds match conditional Poisson inference", {
  for (level in c(0.2, 0.8, 0.95)) {
    for (side in c("two.sided", "lower", "upper")) {
      alternative <- switch(side, two.sided = "two.sided",
                            lower = "greater", upper = "less")
      result <- rate.2s.ci(c(9, 12), c(2, 7), method = "exact",
                           conf.level = level, side = side)
      reference <- poisson.test(c(9, 12), c(2, 7), conf.level = level,
                                 alternative = alternative)
      expect_equal(as.numeric(result$conf.int),
                   as.numeric(reference$conf.int), tolerance = 1e-10)
    }
  }
})

test_that("new one-sided methods use one tail rather than two-sided limits", {
  tab <- matrix(c(7, 27, 1, 33), 2, byrow = TRUE)
  calls <- list(
    function(...) rate.1s.ci(5, 10, method =
                              c("exact", "score", "wh", "wald", "log"), ...),
    function(...) rate.2s.ci(c(12, 5), c(100, 80),
                              method = c("log", "score", "exact"), ...),
    function(...) or.2x2.ci(tab, method = c("wald", "gart", "cornfield"), ...),
    function(...) or.paired.ci(tab, method =
                                c("wald", "wilson", "clopper-pearson",
                                  "midp"), ...),
    function(...) prop.paired.ci(7, 3, 20, method =
                                  c("wald", "waldcc", "agresti-min"), ...),
    function(...) or.common.ci(list(tab, tab), ...)
  )
  for (call in calls) {
    two <- call(conf.level = 0.9)$conf.int
    lower <- call(side = "lower")$conf.int
    upper <- call(side = "upper")$conf.int
    expect_equal(unname(lower[, 1]), unname(two[, 1]), tolerance = 1e-7)
    expect_equal(unname(upper[, 2]), unname(two[, 2]), tolerance = 1e-7)
  }
  for (level in c(0.2, 0.8, 0.95)) {
    for (side in c("lower", "upper")) {
      bp <- or.2x2.ci(tab, method = "bp", side = side, conf.level = level)
      cf <- or.2x2.ci(tab, method = "cornfield", side = side,
                       conf.level = level)
      expect_equal(as.numeric(bp$conf.int), as.numeric(cf$conf.int))
      blaker <- or.paired.ci(tab, method = "blaker", side = side,
                             conf.level = level)
      cp <- or.paired.ci(tab, method = "clopper-pearson", side = side,
                         conf.level = level)
      expect_equal(as.numeric(blaker$conf.int), as.numeric(cp$conf.int))
    }
  }
})

test_that("paired score one-sided limits invert the signed score", {
  estimate <- (7 - 3) / 20
  for (level in c(0.2, 0.8, 0.95)) {
    for (side in c("lower", "upper")) {
      result <- prop.paired.ci(7, 3, 20, side = side, conf.level = level)
      bound <- result$conf.int[if (side == "lower") 1 else 2]
      signed <- sign(estimate - bound) *
        ibist:::paired_score_stat(7, 3, 20, bound)
      expect_equal(signed, if (side == "lower") qnorm(level) else -qnorm(level),
                   tolerance = 1e-6)
    }
  }
})

test_that("low-confidence upper Wald rate bounds stay in the rate domain", {
  result <- rate.1s.ci(0, 10, method = "wald", conf.level = 0.2,
                       side = "upper")
  expect_equal(as.numeric(result$conf.int), c(0, 0))
})
