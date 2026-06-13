test_that("gof_logistic_pearson handles binary responses with weights", {
  data(rds, package = "ibist")
  fit <- glm(death ~ surf + bwt, data = rds, weights = count,
             family = "binomial")

  result <- gof_logistic_pearson(fit, pool = FALSE)

  expect_s3_class(result, "htest")
  expect_named(result$statistic, "X2")
  expect_true(result$parameter > 0)
  expect_equal(sum(result$observed), sum(rds$death * rds$count))
  expect_equal(sum(result$group_n), sum(rds$count))
})

test_that("gof_logistic_pearson handles grouped two-column responses", {
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
             data = grouped, family = "binomial")
  expect_warning(
    result <- gof_logistic_pearson(fit, pool = FALSE),
    "No replication detected"
  )

  expect_s3_class(result, "htest")
  expect_equal(sum(result$observed), sum(grouped$death_count))
  expect_equal(sum(result$group_n),
               sum(grouped$death_count + grouped$alive_count))
})

test_that("gof_logistic_pearson handles proportion responses with weights", {
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

  fit <- glm(death_prop ~ surf + bwt, data = grouped,
             weights = total, family = "binomial")
  expect_warning(
    result <- gof_logistic_pearson(fit, pool = FALSE),
    "No replication detected"
  )

  expect_s3_class(result, "htest")
  expect_equal(sum(result$observed), sum(grouped$death_count))
  expect_equal(sum(result$group_n), sum(grouped$total))
})

test_that("gof_logistic_pearson validates model type", {
  expect_error(gof_logistic_pearson(lm(mpg ~ wt, data = mtcars)),
               "glm")
})
