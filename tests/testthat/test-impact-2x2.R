test_that("book cohort attributable measures use the stated definitions", {
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  p <- tab[, 1]/rowSums(tab)
  e <- rowSums(tab)[1]/sum(tab)
  exposed <- impact.2x2(tab, measure = "fraction-exposed")
  expect_s3_class(exposed, "ibist_impact")
  expect_equal(unname(exposed$estimate), 1-p[2]/p[1])
  # epiR 2.0.93 AFRisk.strata.wald, using two identical strata.
  expect_equal(as.numeric(exposed$conf.int),
    c(0.06416131558155, 0.45516607837169), tolerance = 1e-12)
  risk <- impact.2x2(tab, measure = "risk-exposed")
  expect_equal(unname(risk$estimate), p[1]-p[2])
  expect_equal(as.numeric(risk$conf.int),
    as.numeric(prop.2s.ci(tab, method = "wald")$conf.int))
  population <- impact.2x2(tab, measure = "risk-population")
  expect_equal(unname(population$estimate), e*(p[1]-p[2]))
  # Independent parameterization includes uncertainty in exposure prevalence.
  v <- p*(1-p)/rowSums(tab)
  var.e <- e*(1-e)/sum(tab)
  variance <- e^2*sum(v) + (p[1]-p[2])^2*var.e
  expect_equal(as.numeric(population$conf.int),
    unname(population$estimate) + c(-1, 1)*qnorm(.975)*sqrt(variance))
  known <- impact.2x2(tab, measure = "risk-population", exposure = e)
  expect_false(known$exposure.estimated)
  expect_equal(unname(known$estimate), unname(population$estimate))
  expect_lt(diff(known$conf.int[1, ]), diff(population$conf.int[1, ]))
})

test_that("population fraction variance includes the shared data covariance", {
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  p <- tab[, 1]/rowSums(tab)
  e <- rowSums(tab)[1]/sum(tab)
  risk <- e*p[1] + (1-e)*p[2]
  gradient <- c(-e/risk, 1/p[2]-(1-e)/risk, -(p[1]-p[2])/risk)
  variances <- c(p*(1-p)/rowSums(tab), e*(1-e)/sum(tab))
  expected <- -expm1(log(p[2]/risk) +
    c(1, -1)*qnorm(.975)*sqrt(sum(gradient^2*variances)))
  result <- impact.2x2(tab, measure = "fraction-population")
  expect_equal(unname(result$estimate), 1-p[2]/risk)
  expect_equal(as.numeric(result$conf.int), expected)
  expect_equal(as.numeric(result$conf.int), c(.001924318874, .068189204257),
                 tolerance = 1e-10)
  known <- impact.2x2(tab, measure = "fraction-population",
                               exposure = e)
  conditional <- -expm1(log(p[2]/risk) + c(1, -1)*qnorm(.975)*
    sqrt(sum(gradient[1:2]^2*variances[1:2])))
  expect_equal(as.numeric(known$conf.int), conditional)
})

test_that("case-control fractions match the published Jewell example", {
  # Jewell (2004), p. 90: coffee drinking and pancreatic cancer.
  tab <- matrix(c(347, 555, 20, 88), 2, byrow = TRUE)
  result <- impact.2x2(tab, measure = "fraction-population",
                                design = "case-control", rare = TRUE)
  estimate <- (347*88-555*20)/(88*367)
  variance <- 347/(20*367) + 555/(88*643)
  expected <- 1-exp(log(1-estimate) + c(1, -1)*qnorm(.975)*sqrt(variance))
  expect_equal(unname(result$estimate), estimate)
  expect_equal(as.numeric(result$conf.int), expected)
  expect_equal(round(as.numeric(result$conf.int), 2), c(.36, .75))
  expect_match(result$approximation, "Rare-outcome")
  coffee <- matrix(c(265, 483, 112, 280), 2, byrow = TRUE)
  population <- impact.2x2(coffee, measure = "fraction-population",
                                    design = "case-control", rare = TRUE)
  expect_equal(as.numeric(population$conf.int),
    c(0.02973564892762, 0.32454521426735), tolerance = 1e-12)
  exposed <- impact.2x2(coffee, measure = "fraction-exposed",
                         design = "case-control", rare = TRUE)
  odds <- or.2x2.ci(coffee, method = "wald")
  expect_equal(unname(exposed$estimate), unname(1-1/odds$estimate))
  expect_equal(as.numeric(exposed$conf.int), as.numeric(1-1/odds$conf.int))
})

test_that("sides and zero counts preserve protective associations", {
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  protective <- impact.2x2(tab[2:1, ], measure = "fraction-exposed")
  expect_lt(unname(protective$estimate), 0)
  lower <- impact.2x2(tab, side = "lower")
  expect_equal(lower$conf.int[1, "upper"], 1)
  expect_equal(lower$conf.int,
    impact.2x2(tab, alternative = "greater")$conf.int)
  upper <- impact.2x2(tab, measure = "fraction-population",
                               side = "upper")
  expect_equal(upper$conf.int[1, "lower"], -Inf)
  zero <- matrix(c(0, 20, 4, 16), 2, byrow = TRUE)
  expect_warning(result <- impact.2x2(zero), "unavailable")
  expect_true(all(is.na(result$conf.int["fraction-exposed", ])))
  expect_true(all(is.finite(result$conf.int["risk-exposed", ])))
  expect_match(result$notes["fraction-exposed"], "positive event counts")
  score <- impact.2x2(zero, measure = "fraction-exposed",
                      rr.method = "koopman")
  expect_true(is.na(score$estimate))
  expect_equal(score$conf.int[1, "lower"], -Inf)
  no.baseline <- impact.2x2(zero[2:1, ], measure = "fraction-exposed",
                            rr.method = "koopman")
  expect_equal(unname(no.baseline$estimate), 1)
  expect_equal(no.baseline$conf.int[1, "upper"], 1)
  zero[2, ] <- c(0, 20)
  score <- impact.2x2(zero, measure = "fraction-exposed",
                      rr.method = "koopman")
  expect_equal(as.numeric(score$conf.int), c(-Inf, 1))
  expect_warning(result <- impact.2x2(zero), "unavailable")
  expect_true(is.na(result$estimate["fraction-population"]))

})

test_that("invalid designs and unsupported inputs fail before inference", {
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  expect_error(impact.2x2(tab, design = "case-control"), "rare")
  expect_error(impact.2x2(tab, measure = "risk-exposed",
    design = "case-control", rare = TRUE), "absolute risks")
  expect_error(impact.2x2(tab, design = "case-control", rare = TRUE,
    rr.method = "koopman"), "only to the cohort")
  expect_error(impact.2x2(tab, exposure = NA_real_), "known proportion")
  expect_error(impact.2x2(tab, side = "upper", alternative = "greater"),
                 "conflicting")
  expect_error(impact.2x2(tab + .1), "integer counts")
  expect_error(impact.2x2(tab, prevalence = .1), "only for case-control")
  expect_error(impact.2x2(tab, design = "case-control", exposure = .1),
                 "only for cohort")
  expect_error(impact.2x2(tab, design = "case-control", prevalence = 0),
                 "strictly between")
  expect_error(impact.2x2(tab, measure = character()), "empty")
  expect_error(impact.2x2(tab, conf.level = 1), "conf.level")
  expect_error(impact.2x2(tab, tol = -1), "tol")
  expect_error(impact.2x2(tab, rare = NA), "rare")
  expect_error(impact.2x2(matrix(c(0, 10, 0, 10), 2, byrow = TRUE),
    design = "case-control", rare = TRUE), "totals must be positive")
  repeated <- impact.2x2(tab, measure = rep("risk-population", 2))
  expect_equal(repeated$conf.int[1, ], repeated$conf.int[2, ])
  expect_equal(nrow(as.data.frame(repeated)), 2L)
})

test_that("summary defaults retain every identified measure", {
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  result <- impact.2x2(tab)
  expect_length(result$estimate, 4L)
  expect_named(result$estimate, result$measure)
  expect_named(as.data.frame(result),
    c("measure", "estimate", "lower", "upper", "method"))
  expect_true(result$exposure.estimated)
  expect_true(result$prevalence.estimated)
  expect_output(print(result), "Epidemiological impact")
  case <- impact.2x2(tab, design = "case-control", rare = TRUE)
  expect_identical(case$measure,
    c("fraction-exposed", "fraction-population"))
  expect_output(print(case), "Rare-outcome")
})

test_that("known case-control prevalence reconstructs population risks", {
  # Population exposure .25, risks .3/.1; sampled case fraction is 1/3.
  tab <- matrix(c(85, 70, 85, 270), 2, byrow = TRUE)
  result <- impact.2x2(tab, design = "case-control", prevalence = .15)
  expect_equal(unname(result$estimate), c(.2, .05, 2/3, 1/3))
  expect_equal(result$exposure, .25)
  expect_equal(result$prevalence, .15)
  expect_true(result$exposure.estimated)
  expect_false(result$prevalence.estimated)
  expect_null(result$approximation)
  # Independent numerical derivatives on the two exposure probabilities.
  q <- tab[1, ]/colSums(tab)
  centers <- function(q) {
    e <- .15*q[1] + .85*q[2]
    p1 <- .15*q[1]/e
    p0 <- .15*(1-q[1])/(1-e)
    c(p1-p0, .15-p0, log(p0/p1), log(p0/.15))
  }
  h <- 1e-6
  gradients <- vapply(1:2, function(i) {
    delta <- numeric(2)
    delta[i] <- h
    (centers(q+delta)-centers(q-delta))/(2*h)
  }, numeric(4))
  variances <- as.numeric(gradients^2 %*% (q*(1-q)/colSums(tab)))
  raw <- cbind(centers(q)-qnorm(.975)*sqrt(variances),
               centers(q)+qnorm(.975)*sqrt(variances))
  raw[3:4, ] <- -expm1(raw[3:4, 2:1])
  expect_equal(as.numeric(result$conf.int), as.numeric(raw),
                 tolerance = 1e-9)
  expect_equal(impact.2x2(tab, design = "case-control", prevalence = .15,
    rare = TRUE)$conf.int, result$conf.int)
  boundary <- matrix(c(10, 5, 0, 15), 2, byrow = TRUE)
  expect_warning(zero <- impact.2x2(boundary, design = "case-control",
    prevalence = .1), "unavailable")
  expect_equal(unname(zero$estimate[3:4]), c(1, 1))
  expect_true(all(is.finite(zero$conf.int[1:2, ])))
  expect_true(all(is.na(zero$conf.int[3:4, ])))
})
