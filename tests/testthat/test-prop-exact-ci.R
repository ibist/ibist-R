test_that("Chan-Zhang reproduces published comparison limits", {
  ## WUSS 2016, Table 2, p. 9: search-extracted values, still awaiting
  ## confirmation against the original PDF. Independent checks follow below.
  cases <- list(list(x = c(56, 48), n = c(70, 80), ci = c(0.0382, 0.3453)),
                list(x = c(9, 3), n = c(10, 10), ci = c(0.1608, 0.8795)),
                list(x = c(10, 0), n = c(10, 20), ci = c(0.6915, 1)))
  for (case in cases) {
    result <- prop.2s.ci(case$x, case$n, method = "chan-zhang")
    expect_equal(round(as.numeric(result$conf.int), 4), case$ci)
    expect_equal(unname(result$estimate), diff(rev(case$x / case$n)))
  }
})

test_that("exact intervals preserve reversal, alternatives, and method order", {
  x <- c(4, 1); n <- c(6, 5)
  result <- prop.2s.ci(x, n, method = c("mee", "chan-zhang"))
  reversed <- prop.2s.ci(rev(x), rev(n), method = "chan-zhang")
  expect_equal(as.numeric(result$conf.int[2L, ]),
               -rev(as.numeric(reversed$conf.int)), tolerance = 1e-8)
  expect_equal(rownames(result$conf.int), c("mee", "chan-zhang"))
  for (alternative in c("less", "greater")) {
    one <- prop.2s.ci(x, n, conf.level = 0.975, method = "chan-zhang",
                      alternative = alternative)$conf.int
    if (alternative == "less") {
      expect_equal(one[1L, "lower"], -1)
      expect_equal(one[1L, "upper"], result$conf.int[2L, "upper"])
    } else {
      expect_equal(one[1L, "lower"], result$conf.int[2L, "lower"])
      expect_equal(one[1L, "upper"], 1)
    }
  }
  expect_equal(prop.2s.ci(cbind(x, n - x), method = "chan-zhang"),
               prop.2s.ci(x, n, method = "chan-zhang"))
})

test_that("exact bounds cover all outcomes in small binomial experiments", {
  ## Enumerate outcomes and their actual probabilities; no other CI routine
  ## or copy of the production search is used for this coverage check.
  for (n in list(c(1, 1), c(2, 3), c(3, 2))) {
    counts <- as.matrix(expand.grid(0:n[1L], 0:n[2L]))
    for (level in c(0.8, 0.95)) {
      bounds <- t(apply(counts, 1L, function(x) {
        as.numeric(prop.2s.ci(x, n, method = "chan-zhang",
                              conf.level = level, tol = 1e-6)$conf.int)
      }))
      for (p1 in seq(0, 1, by = 0.1)) {
        for (p2 in seq(0, 1, by = 0.1)) {
          mass <- dbinom(counts[, 1L], n[1L], p1) *
            dbinom(counts[, 2L], n[2L], p2)
          expect_gte(sum(mass[bounds[, 1L] <= p1 - p2 &
                               bounds[, 2L] >= p1 - p2]), level - 1e-12)
        }
      }
    }
  }
})

test_that("constant responses and tolerance controls are handled explicitly", {
  for (x in list(c(0, 0), c(3, 5), c(3, 0), c(0, 5))) {
    result <- prop.2s.ci(x, c(3, 5), method = "chan-zhang")
    expect_true(all(is.finite(result$conf.int)))
    expect_lte(result$conf.int[1L, "lower"], unname(result$estimate))
    expect_gte(result$conf.int[1L, "upper"], unname(result$estimate))
  }
  for (tol in list(0, -1, NA, Inf, c(1e-8, 1e-6))) {
    expect_error(prop.2s.ci(c(1, 2), c(3, 4), nuisance.tol = tol),
                 "nuisance.tol")
  }
  expect_error(prop.2s.ci(c(1, 2), c(1000, 1000), method = "chan-zhang"),
               "1,000,000")
  expect_error(prop.2s.ci(c(1, 2), c(3, 4), method = "chan-zhang",
                          correct = TRUE), "only")
})


test_that("the outer search retains a disconnected accepted component", {
  ## Hawila and Berg (2023), Figure 1, uses x=(5,10), n=(9,19).
  ## Direct likelihood optimization and full table summation independently
  ## establish acceptance at -0.095 and rejection at the larger -0.09.
  x <- c(5, 10); n <- c(9, 19)
  support <- as.matrix(expand.grid(0:n[1L], 0:n[2L]))
  tail.probability <- function(delta, q) {
    score <- apply(support, 1L, function(counts) {
      loglik <- function(p) sum(dbinom(counts, n, c(p + delta, p), log = TRUE))
      range <- c(max(0, -delta), min(1, 1 - delta))
      fit <- optimize(loglik, range, maximum = TRUE, tol = 1e-12)
      candidates <- c(range, fit$maximum)
      fitted <- candidates[which.max(vapply(candidates, loglik, 0.0))]
      p <- c(fitted + delta, fitted)
      variance <- sum(p * (1 - p) / n)
      (counts[1L] / n[1L] - counts[2L] / n[2L] - delta) / sqrt(variance)
    })
    observed <- score[x[1L] + 1L + x[2L] * (n[1L] + 1L)]
    mass <- dbinom(support[, 1L], n[1L], q + delta) *
      dbinom(support[, 2L], n[2L], q)
    sum(mass[score >= observed - 1e-9])
  }
  expect_gt(tail.probability(-0.095, 0.095), 0.44)
  ## At -0.09, even a dense independent nuisance grid is below 0.44.
  ## This is a diagnostic, not the production global maximization algorithm.
  expect_lt(max(vapply(seq(0.09, 1, length.out = 31L), function(q) {
    tail.probability(-0.09, q)
  }, 0.0)), 0.44)
  result <- prop.2s.ci(x, n, conf.level = 0.12, method = "chan-zhang")
  expect_lt(result$conf.int[1L, "lower"], -0.095)
})


test_that("complete separation inverts the analytic worst-case probability", {
  for (n in list(c(3, 5), c(10, 20))) {
    result <- prop.2s.ci(c(n[1L], 0), n, method = "chan-zhang")
    lower <- result$conf.int[1L, "lower"]
    q <- max(0, min(1 - lower, (n[1L] - n[2L] * lower) / sum(n)))
    probability <- (q + lower)^n[1L] * (1 - q)^n[2L]
    expect_lt(abs(probability - 0.025), 1e-8)
    expect_equal(result$conf.int[1L, "upper"], 1)
  }
})

test_that("one-sided low-confidence bounds search beyond the estimate", {
  result <- prop.2s.ci(c(4, 1), c(6, 5), method = "chan-zhang",
                       conf.level = 0.2, alternative = "greater")
  expect_gt(result$conf.int[1L, "lower"], unname(result$estimate))
  expect_equal(result$conf.int[1L, "upper"], 1)
  reversed <- prop.2s.ci(c(1, 4), c(5, 6), method = "chan-zhang",
                         conf.level = 0.2, alternative = "less")
  expect_equal(as.numeric(result$conf.int),
               -rev(as.numeric(reversed$conf.int)))
})
