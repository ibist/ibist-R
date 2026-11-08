test_that("rate power reproduces published Minitab examples", {
  # Fixed external values and matching settings: RATE_POWER_VALIDATION.md.
  for (i in seq_len(2L)) {
    one <- power.r1s.test(n = c(25, 30)[i], r0 = 15, r1 = 13,
                          alternative = "less")
    expect_lte(abs(one$power - c(0.842947, 0.898200)[i]), 5e-7)
    rate <- c(32, 38)[i]
    size <- c(79, 86)[i]
    two <- power.r2s.test(n = size, r1 = rate, r2 = 35, strict = TRUE)
    expect_lte(abs(two$power - c(0.902793, 0.902550)[i]), 5e-7)
    design <- power.r2s.test(r1 = rate, r2 = 35, power = 0.9, strict = TRUE)
    expect_equal(ceiling(design$n), size)
    expect_lt(power.r2s.test(n = size - 1, r1 = rate, r2 = 35,
                             strict = TRUE)$power, 0.9)
  }
})

test_that("rate designs agree with independent closed-form sample sizes", {
  target <- 0.8
  expected <- (qnorm(0.95) * sqrt(0.1) +
    qnorm(target) * sqrt(0.2))^2 / (2 * (0.2 - 0.1)^2)
  one <- power.r1s.test(r0 = 0.1, r1 = 0.2, power = target,
                        exposure = 2, alternative = "greater")
  expect_s3_class(one, "power.htest")
  expect_equal(one$n, expected, tolerance = 1e-8)
  pooled <- (0.1 * 2 + 0.2 * 6) / 8
  expected <- (qnorm(0.975) * sqrt(pooled * (1 / 2 + 1 / 6)) +
    qnorm(target) * sqrt(0.1 / 2 + 0.2 / 6))^2 / (0.2 - 0.1)^2
  two <- power.r2s.test(r1 = 0.1, r2 = 0.2, power = target,
                        group.rate = 2, exposure = c(2, 3))
  expect_equal(two$n, expected, tolerance = 1e-8)
  expect_equal(two$exposure, c(2, 3))
  expect_equal(two$group.rate, 2)
})

test_that("each missing rate design quantity can be recovered", {
  designs <- list(
    list(fun = power.r1s.test,
         args = list(n = 100, r0 = 0.1, r1 = 0.2, sig.level = 0.05)),
    list(fun = power.r2s.test,
         args = list(n = 100, r1 = 0.1, r2 = 0.2, sig.level = 0.05)))
  for (design in designs) {
    target <- do.call(design$fun, design$args)$power
    for (name in names(design$args)) {
      args <- design$args
      args[name] <- list(NULL)
      args$power <- target
      answer <- do.call(design$fun, args)
      expect_equal(answer[[name]], design$args[[name]], tolerance = 1e-7)
    }
  }
  expect_equal(power.r1s.test(n = 20, r0 = 0.2, r1 = 0.2,
                             strict = TRUE)$power, 0.05)
  expect_equal(power.r2s.test(n = 20, r1 = 0.2, r2 = 0.2,
                             strict = TRUE)$power, 0.05)
})

test_that("rate power rejects ambiguous and unattainable designs", {
  expect_error(power.r1s.test(n = 10, r0 = 0.1), "exactly one")
  expect_error(power.r1s.test(n = 10, r0 = 0, r1 = 0.2), "positive")
  expect_error(power.r2s.test(n = 10, r1 = 0.1, r2 = 0.2,
                             exposure = c(1, 2, 3)), "exposure")
  expect_error(power.r2s.test(r1 = 0.2, r2 = 0.2, power = 0.8), "equal")
  expect_error(power.r1s.test(r0 = 0.2, r1 = 0.1, power = 0.8,
                             alternative = "greater"), "could not solve")
})
