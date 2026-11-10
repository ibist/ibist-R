test_that("GLM tests reproduce published SAS examples", {
  d <- utils::read.csv(test_path("fixtures", "glm-gof-sas-mroz.csv"))
  fit <- glm(observed ~ ., d, family = binomial(),
              control = glm.control(epsilon = 1e-10))
  methods <- c("pearson", "deviance", "hosmer-lemeshow", "osius-rojek",
                "stukel")
  result <- suppressWarnings(glm.gof(fit, methods))
  expected <- c(538.4877, 606.5035, 13.2967, 0.8412, 2.2581)
  values <- vapply(result, function(z) unname(z$statistic), numeric(1))
  expect_true(all(abs(values - expected) <= 5.1e-5))
  expect_equal(vapply(result, function(z) unname(z$parameter), numeric(1)),
                 setNames(c(544, 544, 8, 1, 2), methods))
  expect_equal(result[[3]]$groups$n, c(rep(55, 9), 57))
  expect_equal(nrow(as.data.frame(result)), 5L)
})

test_that("Pearson and deviance match independent native statistics", {
  d <- data.frame(y = c(42, 37, 1, 101, 73, 14),
    car = factor(rep(c("small", "medium", "large"), 2)),
    age = factor(rep(1:2, each = 3)),
    exposure = c(500, 1200, 100, 400, 500, 300))
  fit <- glm(y ~ car + age + offset(log(exposure)), d, family = poisson())
  result <- suppressWarnings(glm.gof(fit, c("pearson", "deviance")))
  expect_equal(unname(result[[1]]$statistic), 2.8416, tolerance = 2e-5)
  expect_equal(unname(result[[2]]$statistic), 2.8207, tolerance = 2e-5)
  expect_equal(unname(result[[1]]$statistic),
                 sum(residuals(fit, type = "pearson")^2))
  expect_equal(unname(result[[2]]$statistic), deviance(fit))
})

test_that("pre and post diagnostics agree and pre-checks do not run tests", {
  d <- data.frame(x = 1:6, y = c(1, 3, 4, 2, 5, 7), n = rep(50, 6))
  fit <- glm(cbind(y, n-y) ~ x, d, family = binomial())
  result <- glm.gof(fit, "stukel", asymcheck = TRUE)
  expect_equal(unname(result$parameter), 1)
  expect_equal(result$diagnostics, asymcheck(result))
  expect_equal(result$diagnostics, asymcheck(glm.gof, fit, method = "stukel"))
  repeated <- glm.gof(fit, c("stukel", "stukel"))
  expect_equal(repeated[[1]], repeated[[2]])
  expect_equal(asymcheck(repeated),
                 asymcheck(glm.gof, fit, method = c("stukel", "stukel")))
  local_mocked_bindings(.glm_gof_stat = function(...) stop("test ran"))
  expect_s3_class(asymcheck(glm.gof, fit, method = "stukel"),
                   "ibist_asymcheck")
})

test_that("offsets, aliasing, and response encodings are preserved", {
  d <- data.frame(x = rep(0:2, each = 2), y = c(3, 5, 9, 10, 13, 16),
    n = rep(20, 6), off = rep(c(0, 0.5), 3))
  fit <- glm(cbind(y, n-y) ~ x + offset(off), d, family = binomial())
  alias <- glm(I(y/n) ~ x + I(2*x) + offset(off), d, weights = n,
                family = binomial())
  first <- suppressWarnings(glm.gof(fit))
  second <- suppressWarnings(glm.gof(alias))
  expect_equal(first$groups, second$groups)
  expect_equal(unname(first$parameter), 4)
  expect_equal(first$statistic, second$statistic)
})

test_that("unsupported families and singular tests fail clearly", {
  d <- data.frame(x = 0:7, y = c(0, 1, 1, 0, 0, 1, 1, 0))
  fit <- glm(y ~ x, d, family = binomial())
  expect_error(glm.gof(fit, "osius-rojek"), "variance")
  expect_error(glm.gof(fit, "stukel"), "degrees of freedom")
  expect_error(glm.gof(fit, "hosmer-lemeshow"), "degrees of freedom")
  check <- asymcheck(glm.gof, fit, method = "osius-rojek")
  expect_true(is.na(check$condition.met))
  expect_false(check$details$computability$ok)
  expect_error(glm.gof(glm(y ~ x, data = d)), "dispersion one")
  expect_error(glm.gof(fit, "unknown"), "method")
  expect_error(glm.gof(fit, g = 4), "Hosmer-Lemeshow")
  fit$converged <- FALSE
  expect_error(glm.gof(fit), "converged")
})

test_that("zero-weight boundary predictions do not affect GOF", {
  d <- data.frame(x = c(0:5, 1000), y = c(2, 4, 6, 7, 9, 12, 0)/20,
                  w = c(rep(20, 6), 0))
  fit <- suppressWarnings(glm(y ~ x, data = d, weights = w,
                               family = binomial()))
  base <- glm(y ~ x, data = d[1:6, ], weights = w, family = binomial())
  expect_equal(suppressWarnings(glm.gof(fit))$statistic,
                 suppressWarnings(glm.gof(base))$statistic)
})
