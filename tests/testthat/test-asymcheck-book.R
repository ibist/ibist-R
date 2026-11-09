test_that("book binomial criteria distinguish null tests and estimation", {
  # Ch. 11: 74 observations, null probability .10 gives variance 6.66.
  one <- asymcheck(stats::prop.test, x = 16, n = 74, p = 0.1)
  expect_equal(one$checks$value, 6.66)
  expect_true(one$condition.met)
  # Positive expected counts do not suffice: the variance must reach five.
  expect_false(asymcheck(stats::prop.test, x = 8, n = 16)$condition.met)
  two <- asymcheck(stats::prop.test, x = c(0, 10), n = c(10, 10))
  expect_equal(two$checks$value, rep(5, 4))
  expect_true(two$condition.met)
  expect_equal(two, asymcheck(stats::prop.test,
    x = matrix(c(0, 10, 10, 0), 2, byrow = TRUE)))
  exact <- asymcheck(stats::binom.test, x = c(0, 5))
  expect_false(exact$applicable)
  expect_true(is.na(exact$condition.met))
  expect_error(asymcheck(stats::prop.test(10, 20)), "pre-check")
})

test_that("original binomial variance drives all approximate CI screens", {
  # Ch. 11 OSA: 16/40 gives 9.6; boundary 0/40 fails without pseudo-counts.
  methods <- c("wald", "wilson", "agresti-coull", "clopper-pearson", "wald")
  args <- list(x = c(0, 16), n = 40, method = methods, correct = TRUE,
               midp = TRUE, side = "lower")
  pre <- do.call(asymcheck, c(list(object = prop.1s.ci), args))
  fit <- do.call(prop.1s.ci, c(args, list(asymcheck = TRUE)))
  plain <- do.call(prop.1s.ci, args)
  expect_equal(pre, fit$diagnostics)
  expect_equal(pre, asymcheck(fit))
  expect_equal(pre, asymcheck(plain))
  expect_equal(fit$conf.int, plain$conf.int)
  expect_equal(pre$checks$value, c(rep(0, 5), rep(9.6, 5)))
  expect_false(pre$condition.met)
  expect_identical(pre$method, methods)
  expect_equal(pre$checks$threshold, rep(c(5, 5, 5, NA, 5), 2))
  # All cells positive is insufficient: n=4, p.hat=.5 gives only 1.
  small <- matrix(2, 2, 2)
  expect_false(asymcheck(rr.2x2.ci, x = small, method = "katz")$condition.met)
  expect_false(asymcheck(or.2x2.ci, x = small, method = "gart")$condition.met)
  expect_false(asymcheck(prop.2s.ci, x = small,
                         method = "agresti-caffo")$condition.met)
  boundary <- matrix(c(0, 40, 16, 24), 2, byrow = TRUE)
  expect_false(asymcheck(or.2x2.ci, x = boundary,
                         method = "gart")$condition.met)
  expect_true(asymcheck(prop.1s.ci, x = 10, n = 20)$condition.met)
})

test_that("book VF risk-ratio and odds-ratio examples agree", {
  # Ch. 15 VF: original row variances are 23.56 and 194.5518292683.
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  expected <- c(23.56, 267*(717/984))
  requests <- list(
    list(fun = rr.2x2.ci, method = c("katz", "koopman", "katz")),
    list(fun = or.2x2.ci, method = c("wald", "gart", "cornfield", "bp")),
    list(fun = prop.2s.ci, method = c("newcombe", "miettinen-nurminen",
                                     "mee", "agresti-caffo", "brown-li")))
  for (request in requests) {
    args <- list(x = tab, method = request$method, side = "upper")
    pre <- do.call(asymcheck, c(list(object = request$fun), args))
    fit <- do.call(request$fun, c(args, list(asymcheck = TRUE)))
    plain <- do.call(request$fun, args)
    expect_equal(pre$checks$value, rep(expected, length(request$method)))
    expect_true(pre$condition.met)
    expect_equal(pre, fit$diagnostics)
    expect_equal(pre, asymcheck(plain))
    expect_equal(fit$conf.int, plain$conf.int)
  }
  exact <- asymcheck(or.2x2.ci, x = tab, method = c("cornfield", "bp"),
                     midp = TRUE, tol = 1e-7)
  expect_false(exact$applicable)
  expect_true(is.na(exact$condition.met))
  expect_true(all(is.na(exact$checks$threshold)))
})

test_that("rate tests use null quantities and CI checks use observed events", {
  # Ch. 21 AF: testing uses .0119*25800, not the 411 observed events.
  one <- asymcheck(rate.test, x = 411, exposure = 25800, r = 0.0119)
  expect_equal(one$checks$value, 307.02)
  expect_true(one$condition.met)
  expect_false(asymcheck(rate.test, x = 411, exposure = 1, r = 9)$condition.met)
  # Ch. 22 low-HDL: 21*.1951479*(1-.1951479) < 5.
  x <- c(9, 12)
  exposure <- c(1817.6, 7496.3)
  two <- asymcheck(rate.test, x = x, exposure = exposure)
  pi0 <- 1817.6/(1817.6+7496.3)
  expect_equal(two$checks$value, 21*pi0*(1-pi0))
  expect_false(two$condition.met)
  fit <- suppressWarnings(rate.test(x, exposure, asymcheck = TRUE))
  expect_equal(two, fit$diagnostics)
  expect_equal(two, asymcheck(fit))
  # Ch. 22 MI: both 151 and 55 exceed ten.
  args <- list(x = c(151, 55), exposure = c(57518.1, 74573.5),
               method = c("log", "score", "exact"))
  pre <- do.call(asymcheck, c(list(object = rate.2s.ci), args))
  fit <- do.call(rate.2s.ci, c(args, list(asymcheck = TRUE)))
  expect_equal(pre, fit$diagnostics)
  expect_equal(pre, asymcheck(fit))
  expect_true(pre$condition.met)
  expect_false(asymcheck(rate.2s.ci, x = c(10, 11))$condition.met)
  expect_true(asymcheck(rate.2s.ci, x = c(11, 11))$condition.met)
})

test_that("one-rate book guidelines distinguish interval constructions", {
  expect_false(asymcheck(rate.1s.ci, x = 49, method = "wald")$condition.met)
  expect_true(asymcheck(rate.1s.ci, x = 50, method = "wald")$condition.met)
  expect_false(asymcheck(rate.1s.ci, x = 9, method = "log")$condition.met)
  expect_true(asymcheck(rate.1s.ci, x = 10, method = "log")$condition.met)
  expect_true(asymcheck(rate.1s.ci, x = 10, method = "score",
                        correct = FALSE)$condition.met)
  # The book recommends WH for small counts; don't invent a cutoff for it.
  expect_true(is.na(asymcheck(rate.1s.ci, x = 3, method = "wh")$condition.met))
  expect_false(asymcheck(rate.1s.ci, x = 0, method = "log")$applicable)
  # Ch. 12 uses all four expected cells >=5 for a 2x2 test.
  expect_true(asymcheck(matrix(c(1, 9, 9, 1), 2),
                        method = "pearson")$condition.met)
  tab <- matrix(c(1, 1, 5, 5), 2, byrow = TRUE)
  expect_false(asymcheck(tab, method = "pearson")$condition.met)
})

test_that("book pre-checks do not run interval or test backends", {
  fail <- function(...) stop("inference ran")
  local_mocked_bindings(binomial_ci_limits = fail,
    independent_exact_limits = fail, koopman_ratio_lower = fail,
    or_ci_cornfield = fail, or_ci_baptista_pike = fail,
    ci_rate_ratio_log = fail)
  tab <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
  expect_s3_class(asymcheck(prop.1s.ci, x = 2, n = 10,
                            method = "blaker"), "ibist_asymcheck")
  expect_false(asymcheck(prop.2s.ci, x = tab,
                         method = "chan-zhang")$applicable)
  expect_true(asymcheck(rr.2x2.ci, x = tab)$condition.met)
  expect_false(asymcheck(or.2x2.ci, x = tab, method = "bp")$applicable)
  expect_true(asymcheck(rate.2s.ci, x = c(151, 55))$condition.met)
  expect_error(asymcheck(prop.1s.ci, x = 2, n = 10, midp = TRUE), "midp")
  expect_error(asymcheck(or.2x2.ci, x = tab, unused = 1), "additional")
  expect_error(rate.test(1, asymcheck = 1), "asymcheck")
})
