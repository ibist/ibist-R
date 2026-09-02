test_that("logistic_gof handles binary responses with weights", {
  data(rds, package = "ibist")
  fit <- glm(death ~ surf + bwt,
    data = rds, weights = count,
    family = "binomial"
  )

  result <- logistic_gof(fit, pool = FALSE)

  expect_s3_class(result, "htest")
  expect_named(result$statistic, "X2")
  expect_true(result$parameter > 0)
  expect_equal(sum(result$observed), sum(rds$death * rds$count))
  expect_equal(sum(result$group_n), sum(rds$count))
})

test_that("logistic_gof handles grouped two-column responses", {
  data(rds, package = "ibist")
  grouped <- reshape(
    rds,
    idvar = c("bwt", "surf"),
    timevar = "death",
    direction = "wide"
  )
  names(grouped) <- sub("count\\.1", "death_count", names(grouped))
  names(grouped) <- sub("count\\.0", "alive_count", names(grouped))

  fit <- glm(cbind(death_count, alive_count) ~ surf + bwt,
    data = grouped, family = "binomial"
  )
  expect_warning(
    result <- logistic_gof(fit, pool = FALSE),
    "No replication detected"
  )

  expect_s3_class(result, "htest")
  expect_equal(sum(result$observed), sum(grouped$death_count))
  expect_equal(
    sum(result$group_n),
    sum(grouped$death_count + grouped$alive_count)
  )
})

test_that("logistic_gof handles proportion responses with weights", {
  data(rds, package = "ibist")
  grouped <- reshape(
    rds,
    idvar = c("bwt", "surf"),
    timevar = "death",
    direction = "wide"
  )
  names(grouped) <- sub("count\\.1", "death_count", names(grouped))
  names(grouped) <- sub("count\\.0", "alive_count", names(grouped))
  grouped$total <- grouped$death_count + grouped$alive_count
  grouped$death_prop <- grouped$death_count / grouped$total

  fit <- glm(death_prop ~ surf + bwt,
    data = grouped,
    weights = total, family = "binomial"
  )
  expect_warning(
    result <- logistic_gof(fit, pool = FALSE),
    "No replication detected"
  )

  expect_s3_class(result, "htest")
  expect_equal(sum(result$observed), sum(grouped$death_count))
  expect_equal(sum(result$group_n), sum(grouped$total))
})

test_that("logistic_gof validates model type", {
  expect_error(
    logistic_gof(lm(mpg ~ wt, data = mtcars)),
    "glm"
  )
})

test_that("logistic_gof validates control arguments", {
  data(rds, package = "ibist")
  fit <- glm(death ~ surf + bwt,
    data = rds, weights = count,
    family = "binomial"
  )

  expect_error(logistic_gof(fit, min_n = 0), "min_n")
  expect_error(
    logistic_gof(fit, min_expected = -1),
    "min_expected"
  )
  expect_error(logistic_gof(fit, pool = NA), "pool")
})

test_that("logistic_gof validates binomial response encodings", {
  prop_data <- data.frame(y = c(0.25, 0.75), x = c(0, 1))
  fit_prop <- suppressWarnings(
    glm(y ~ x, data = prop_data, family = "binomial")
  )

  expect_error(
    logistic_gof(fit_prop),
    "Proportion responses require binomial weights"
  )
})

test_that("logistic_gof uses fitted model rank for degrees of freedom", {
  data <- data.frame(
    x = rep(1:6, each = 2),
    y = rep(c(0, 1), 6)
  )
  data$x_dup <- data$x
  fit <- glm(y ~ x + x_dup, data = data, family = "binomial")

  result <- logistic_gof(fit, pool = FALSE)

  expect_true(anyNA(coef(fit)))
  expect_equal(unname(result$parameter), 4)
})
