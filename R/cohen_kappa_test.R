#' Asymptotic Inference for Cohen's Kappa
#'
#' Computes Cohen's kappa or weighted kappa from a square table of paired
#' ratings, with a Wald confidence interval and an asymptotic test of
#' \eqn{H_0: \kappa = 0}.
#'
#' @param x A square matrix or table of non-negative integer counts. Rows and
#'   columns represent the same ordered rating categories.
#' @param weights The agreement weights: \code{"unweighted"}, \code{"linear"},
#'   \code{"quadratic"}, or a numeric matrix of user-specified weights. The
#'   named choices use category positions in \code{x}; the quadratic weights
#'   are Fleiss--Cohen weights. A supplied matrix must have the same dimensions
#'   as \code{x}, be symmetric, have ones on the diagonal, and contain values
#'   between zero and one.
#' @param conf.level Confidence level for the two-sided Wald confidence
#'   interval.
#' @param alternative Alternative hypothesis: \code{"two.sided"},
#'   \code{"greater"}, or \code{"less"}.
#' @param test.variance Variance used for the test: \code{"fce"} for the
#'   Fleiss--Cohen--Everitt null variance, or \code{"delta"} for the
#'   multinomial delta-method variance evaluated at the observed proportions.
#'
#' @details
#' The confidence interval uses the Fleiss--Cohen--Everitt (FCE) large-sample
#' variance of the kappa estimate. This variance is algebraically equivalent
#' to the multinomial delta-method variance used by
#' \code{psych::cohen.kappa()}. By default, the test uses the distinct FCE
#' variance under the null hypothesis \eqn{\kappa = 0}. Setting
#' \code{test.variance = "delta"} instead uses the delta-method variance
#' evaluated at the observed cell proportions. Evaluating the delta-method
#' variance under the null gives the FCE null variance, so it is not a separate
#' test option. The interval is an untruncated Wald interval.
#'
#' This function provides asymptotic inference only. SAS \code{PROC FREQ}
#' also provides exact tests for simple and weighted kappa; exact inference is
#' not implemented here.
#'
#' @return An object of class \code{"htest"} containing the estimate,
#'   confidence interval, Z statistic, and p-value. It also includes the
#'   standard error used for the confidence interval, the standard error under
#'   the null hypothesis, the selected test variance method, the weight
#'   matrix, and the sample size.
#'
#' @references
#' Fleiss, J. L., Cohen, J., and Everitt, B. S. (1969). Large sample standard
#' errors of kappa and weighted kappa. \emph{Psychological Bulletin}, 72(5),
#' 323--327. \doi{10.1037/h0028106}
#'
#' SAS Institute Inc. (2015). Tests and measures of agreement. In
#' \emph{Base SAS 9.4 Procedures Guide: Statistical Procedures, Third Edition}.
#' \url{https://support.sas.com/documentation/cdl/en/procstat/67528/HTML/default/procstat_freq_details76.htm}
#'
#' @examples
#' peff <- matrix(
#'   c(581, 134, 12, 6, 0, 2,
#'     91, 73, 18, 31, 3, 1,
#'     4, 18, 15, 19, 3, 3,
#'     1, 8, 7, 24, 11, 6,
#'     0, 0, 2, 4, 3, 5,
#'     0, 0, 0, 9, 11, 13),
#'   nrow = 6,
#'   byrow = TRUE
#' )
#' cohen.kappa.test(peff)
#' cohen.kappa.test(peff, weights = "quadratic")
#' cohen.kappa.test(peff, weights = "quadratic", test.variance = "delta")
#'
#' @export
cohen.kappa.test <- function(
  x,
  weights = "unweighted",
  conf.level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  test.variance = c("fce", "delta")
) {
  alternative <- match.arg(alternative)
  test.variance <- match.arg(test.variance)

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || !is.finite(conf.level) || conf.level <= 0 ||
      conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }

  tab <- as.matrix(x)
  if (length(dim(tab)) != 2L || nrow(tab) != ncol(tab) || nrow(tab) < 2L ||
      !is.numeric(tab) || anyNA(tab) || any(!is.finite(tab)) ||
      any(tab < 0) || any(tab != floor(tab))) {
    stop("'x' must be a square table of non-negative integer counts.")
  }

  n <- sum(tab)
  if (!is.finite(n) || n <= 0) {
    stop("'x' must have a positive total count.")
  }

  weight.method <- NULL
  if (is.character(weights) && length(weights) == 1L && !is.na(weights)) {
    weight.method <- match.arg(
      weights,
      choices = c("unweighted", "linear", "quadratic")
    )
    category <- seq_len(nrow(tab))
    distance <- abs(outer(category, category, "-")) / (nrow(tab) - 1)
    weight.matrix <- switch(
      weight.method,
      unweighted = diag(1, nrow(tab)),
      linear = 1 - distance,
      quadratic = 1 - distance^2
    )
  } else if (is.matrix(weights) && is.numeric(weights)) {
    weight.matrix <- weights
    if (!identical(dim(weight.matrix), dim(tab)) || anyNA(weight.matrix) ||
        any(!is.finite(weight.matrix)) || any(weight.matrix < 0) ||
        any(weight.matrix > 1) ||
        !isTRUE(all.equal(weight.matrix, t(weight.matrix))) ||
        any(diag(weight.matrix) != 1)) {
      stop(paste0(
        "A custom 'weights' matrix must match 'x', be symmetric, have ones ",
        "on the diagonal, and contain finite values between 0 and 1."
      ))
    }
    weight.method <- "custom"
  } else {
    stop("'weights' must be a supported name or a numeric matrix.")
  }

  probability <- tab / n
  row.probability <- rowSums(probability)
  col.probability <- colSums(probability)
  independent <- outer(row.probability, col.probability)
  observed.agreement <- sum(probability * weight.matrix)
  chance.agreement <- sum(independent * weight.matrix)
  denominator <- 1 - chance.agreement

  if (denominator <= 0) {
    stop("Kappa is undefined because chance agreement equals one.")
  }

  estimate <- (observed.agreement - chance.agreement) / denominator
  row.weight <- drop(weight.matrix %*% col.probability)
  col.weight <- drop(crossprod(row.probability, weight.matrix))
  row.col.weight <- outer(row.weight, rep(1, ncol(tab))) +
    outer(rep(1, nrow(tab)), col.weight)

  variance <- (
    sum(probability * (
      weight.matrix - row.col.weight * (1 - estimate)
    )^2) - (estimate - chance.agreement * (1 - estimate))^2
  ) / (denominator^2 * n)

  null.variance <- (
    sum(independent * (weight.matrix - row.col.weight)^2) -
      chance.agreement^2
  ) / (denominator^2 * n)

  tolerance <- 100 * .Machine$double.eps
  if (variance < 0 && variance > -tolerance) variance <- 0
  if (null.variance < 0 && null.variance > -tolerance) null.variance <- 0
  if (!is.finite(variance) || variance < 0 ||
      !is.finite(null.variance) || null.variance <= 0) {
    stop("A valid kappa variance could not be computed for 'x'.")
  }

  test.variance.value <- if (test.variance == "fce") {
    null.variance
  } else {
    variance
  }
  if (test.variance.value <= 0) {
    stop("The selected test variance must be positive for 'x'.")
  }
  std.error <- sqrt(variance)
  null.std.error <- sqrt(test.variance.value)
  z <- estimate / null.std.error
  p.value <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(-abs(z)),
    greater = stats::pnorm(z, lower.tail = FALSE),
    less = stats::pnorm(z)
  )
  z.critical <- stats::qnorm(1 - (1 - conf.level) / 2)
  conf.int <- estimate + c(-1, 1) * z.critical * std.error
  attr(conf.int, "conf.level") <- conf.level

  estimate.name <- if (weight.method == "unweighted") {
    "kappa"
  } else {
    "weighted kappa"
  }
  structure(
    list(
      statistic = c(Z = z),
      p.value = p.value,
      conf.int = conf.int,
      estimate = stats::setNames(estimate, estimate.name),
      null.value = stats::setNames(0, estimate.name),
      alternative = alternative,
      method = if (weight.method == "unweighted") {
        "Asymptotic test for Cohen's kappa"
      } else {
        "Asymptotic test for weighted Cohen's kappa"
      },
      data.name = deparse(substitute(x)),
      std.error = std.error,
      null.std.error = null.std.error,
      test.variance = test.variance,
      weights = weight.matrix,
      weight.method = weight.method,
      n = n
    ),
    class = "htest"
  )
}
