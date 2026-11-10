test_that("grouped binomial checks use variance rather than cell positivity", {
  d <- data.frame(x = 0:3, y = rep(6, 4), n = rep(24, 4))
  fit <- glm(cbind(y, n-y) ~ x, d, family = binomial())
  pre <- asymcheck(glm.gof, fit, method = c("pearson", "deviance"))
  expect_false(pre$condition.met)
  expect_equal(pre$checks$value, c(4.5, 4.5))
  expect_equal(pre$checks$threshold, c(5, 5))
  details <- pre$details$methods[[1L]]
  expect_true(all(details$groups$expected > 5))
  expect_true(all(details$groups$n-details$groups$expected > 5))
  expect_equal(details$binomial.variance, rep(4.5, 4))
  expect_true(details$computability$ok)
  expect_equal(details$pattern.fraction, 4/96)
  fit.gof <- suppressWarnings(glm.gof(fit, c("pearson", "deviance"),
                                      asymcheck = TRUE))
  expect_equal(pre, asymcheck(fit.gof))
  expect_equal(pre$details$methods[[1L]], fit.gof[[1L]]$diagnostics$details)
  # At p=.5 the exact book threshold is reached by n=20.
  d$y <- 10
  d$n <- 20
  boundary <- glm(cbind(y, n-y) ~ x, d, family = binomial())
  expect_true(asymcheck(glm.gof, boundary)$condition.met)
})

test_that("HL checks actual group count after ties rather than requested g", {
  d <- data.frame(x = 0:4, y = seq(20, 60, 10), n = 100)
  fit <- glm(cbind(y, n-y) ~ x, d, family = binomial())
  pre <- asymcheck(glm.gof, fit, method = "hosmer-lemeshow", g = 10)
  expect_equal(pre$details$requested.groups, 10)
  expect_equal(pre$details$actual.groups, 5)
  expect_equal(pre$checks$value[1L], 5)
  expect_equal(pre$checks$threshold[1L], 6)
  expect_true(all(pre$details$binomial.variance > 5))
  expect_false(pre$condition.met)
  expect_warning(result <- glm.gof(fit, "hosmer-lemeshow", asymcheck = TRUE),
                   "approximation criterion")
  expect_equal(pre, result$diagnostics)
  d <- rbind(d, data.frame(x = 5, y = 70, n = 100))
  fit <- glm(cbind(y, n-y) ~ x, d, family = binomial())
  expect_true(asymcheck(glm.gof, fit,
                        method = "hosmer-lemeshow")$condition.met)
})

test_that("Poisson checks report small counts without inventing many cutoff", {
  d <- data.frame(x = 0:3, y = c(0, 1, 2, 3))
  fit <- glm(y ~ x, d, family = poisson())
  check <- asymcheck(glm.gof, fit, method = "deviance")
  expected <- fitted(fit)
  expect_equal(check$checks$value,
    c(min(expected), sum(expected < 2), mean(expected < 5)))
  expect_equal(check$checks$threshold, c(5, NA, NA))
  expect_false(check$condition.met)
  expect_equal(check$details$expected.below.two, sum(expected < 2))
  expect_equal(check$details$fraction.expected.below.five, mean(expected < 5))
  d <- data.frame(x = 0:3, y = rep(6, 4))
  fit <- glm(y ~ x, d, family = poisson())
  expect_true(asymcheck(glm.gof, fit)$condition.met)
})

test_that("computability never assesses Osius-Rojek or Stukel adequacy", {
  d <- data.frame(x = 0:5, y = c(20, 25, 35, 50, 65, 75), n = 100)
  fit <- glm(cbind(y, n-y) ~ x, d, family = binomial())
  for (method in c("osius-rojek", "stukel")) {
    check <- asymcheck(glm.gof, fit, method = method)
    expect_true(check$details$computability$ok)
    expect_true(check$applicable)
    expect_true(is.na(check$condition.met))
    expect_true(all(is.na(check$checks$threshold)))
    result <- glm.gof(fit, method, asymcheck = TRUE)
    expect_equal(check, result$diagnostics)
    expect_equal(check, asymcheck(result))
  }
  mixed <- asymcheck(glm.gof, fit, method = c("pearson", "osius-rojek"))
  expect_true(is.na(mixed$condition.met))
  d$n <- 4
  d$y <- c(1, 1, 1, 2, 3, 3)
  sparse <- glm(cbind(y, n-y) ~ x, d, family = binomial())
  expect_false(asymcheck(glm.gof, sparse,
    method = c("pearson", "osius-rojek"))$condition.met)
  local_mocked_bindings(.glm_gof_stat = function(...) stop("test ran"))
  expect_s3_class(asymcheck(glm.gof, fit, method = "stukel"),
                   "ibist_asymcheck")
})
