#' Confidence Intervals for a Difference of Independent Proportions
#'
#' Computes intervals for the first population proportion minus the second.
#'
#' @param x Numeric vector of two success counts, or a 2 x 2 table whose
#'   rows are independent samples and columns are successes and failures.
#' @param n Numeric vector of two positive sample sizes. Omit for a table.
#' @param conf.level Confidence level, between 0 and 1.
#' @param method One or more of \code{"newcombe"} (default), \code{"wald"},
#'   \code{"miettinen-nurminen"}, \code{"mee"}, \code{"agresti-caffo"},
#'   or \code{"brown-li"}.
#' @param correct Logical. Apply continuity correction to Newcombe and Wald
#'   intervals. Other requested methods are unchanged.
#' @param alternative One of \code{"two.sided"}, \code{"less"}, or
#'   \code{"greater"}, referring to the proportion difference.
#' @param tol Positive numerical tolerance for score inversion.
#'
#' @details
#' Newcombe's hybrid interval combines separate Wilson bounds using the
#' square root of the sum of squared distances from the sample proportions.
#' With \code{correct = TRUE}, it combines continuity-corrected Wilson bounds.
#' This is not inversion of the pooled score test for equal proportions.
#'
#' Wald intervals use the unpooled binomial variance, adding
#' \eqn{1/(2n_1)+1/(2n_2)} to the margin when corrected. Agresti--Caffo adds
#' one success and one failure to each sample before computing a Wald interval.
#' Miettinen--Nurminen inverts the score statistic using constrained binomial
#' maximum likelihood estimates and the variance factor
#' \eqn{(n_1+n_2)/(n_1+n_2-1)}. Mee's method uses the same constrained
#' estimates and score inversion but omits this variance factor. The factor
#' is always included for \code{"miettinen-nurminen"} and never included for
#' \code{"mee"}; it is not controlled by the continuity-correction argument
#' \code{correct}. All limits are restricted to \eqn{[-1,1]}.
#'
#' Brown--Li's Jeffreys-based normal approximation uses
#' \eqn{\bar p_i=(x_i+1/2)/(n_i+1)} and the limits
#' \deqn{\bar p_1-\bar p_2 \mathbin{\pm} z
#'   \sqrt{\bar p_1(1-\bar p_1)/n_1+\bar p_2(1-\bar p_2)/n_2}.}
#' Here \eqn{z} is the standard-normal quantile corresponding to the requested
#' confidence level and alternative. The variance denominators are the
#' original sample sizes, not the adjusted sample sizes or posterior-variance
#' denominators. This approximate interval is not obtained from quantiles of
#' the posterior difference. As for Agresti--Caffo, the reported point estimate
#' remains the observed proportion difference, even though the interval uses
#' adjusted proportions internally.
#'
#' @return An \code{ibist_ci} object with the observed named \code{estimate},
#'   a \code{conf.int} matrix with one row per method, \code{conf.level},
#'   \code{method}, \code{alternative}, and the two \code{sample.proportions}.
#' @references
#' Agresti, A., and Caffo, B. (2000). Simple and effective confidence intervals
#' for proportions and differences of proportions result from adding two
#' successes and two failures. \emph{The American Statistician}, 54, 280--288.
#'
#' Brown, L., and Li, X. (2005). Confidence intervals for two sample binomial
#' distribution. \emph{Journal of Statistical Planning and Inference}, 130,
#' 359--375. \doi{10.1016/j.jspi.2003.09.039}
#'
#' Mee, R. W. (1984). Confidence bounds for the difference between two
#' probabilities. \emph{Biometrics}, 40, 1175--1176.
#'
#' Miettinen, O., and Nurminen, M. (1985). Comparative analysis of two rates.
#' \emph{Statistics in Medicine}, 4, 213--226.
#'
#' Newcombe, R. G. (1998). Interval estimation for the difference between
#' independent proportions: comparison of eleven methods.
#' \emph{Statistics in Medicine}, 17, 873--890.
#' @examples
#' prop.2s.ci(c(1366, 383), c(15571, 5366))
#' prop.2s.ci(c(12, 4), c(40, 40), method = c("newcombe", "wald"))
#' prop.2s.ci(c(12, 4), c(40, 30),
#'             method = c("miettinen-nurminen", "mee"))
#' prop.2s.ci(c(9, 3), c(10, 10), method = c("newcombe", "brown-li"))
#' prop.2s.ci(c(12, 4), c(40, 40), alternative = "greater")
#' tab <- matrix(c(38, 62, 267, 717), nrow = 2, byrow = TRUE)
#' prop.2s.ci(tab)
#' @export
prop.2s.ci <- function(
  x, n = NULL, conf.level = 0.95,
  method = c("newcombe", "wald", "miettinen-nurminen", "mee",
             "agresti-caffo", "brown-li"),
  correct = FALSE, alternative = c("two.sided", "less", "greater"),
  tol = 1e-8
) {
  counts <- independent_binomial_counts(x, n)
  settings <- validate_binomial_ci_settings(conf.level, alternative, tol)
  methods <- if (missing(method)) method[1L] else {
    match.arg(method, several.ok = TRUE)
  }
  validate_ci_flag(correct, "correct")
  if (correct && !any(methods %in% c("wald", "newcombe"))) {
    stop("'correct' applies only to Wald and Newcombe intervals.")
  }
  p <- counts$x / counts$n
  z <- stats::qnorm(settings$tail, lower.tail = FALSE)
  intervals <- lapply(methods, function(method) {
    if (method == "newcombe") {
      one <- binomial_ci_limits(counts$x[1L], counts$n[1L], "wilson",
                                settings$tail, correct, FALSE, tol)
      two <- binomial_ci_limits(counts$x[2L], counts$n[2L], "wilson",
                                settings$tail, correct, FALSE, tol)
      limits <- p[1L] - p[2L] + c(-1, 1) * sign(z) *
        sqrt((one - p[1L])^2 + (rev(two) - p[2L])^2)
    } else if (method %in% c("miettinen-nurminen", "mee")) {
      variance.factor <- if (method == "miettinen-nurminen") {
        sum(counts$n) / (sum(counts$n) - 1)
      } else 1
      limits <- independent_score_limits(counts$x, counts$n, z, tol,
                                          variance.factor)
    } else if (method == "brown-li") {
      adjusted <- (counts$x + 0.5) / (counts$n + 1)
      margin <- z * sqrt(sum(adjusted * (1 - adjusted) / counts$n))
      limits <- adjusted[1L] - adjusted[2L] + c(-1, 1) * margin
    } else {
      size <- counts$n + if (method == "agresti-caffo") 2 else 0
      center <- (counts$x + if (method == "agresti-caffo") 1 else 0) / size
      margin <- z * sqrt(sum(center * (1 - center) / size))
      if (correct && method == "wald") {
        margin <- margin + sum(0.5 / counts$n)
      }
      limits <- center[1L] - center[2L] + c(-1, 1) * margin
    }
    limits <- pmax(-1, pmin(1, limits))
    if (settings$alternative == "greater") limits[2L] <- 1
    if (settings$alternative == "less") limits[1L] <- -1
    limits
  })
  labels <- methods
  corrected <- correct & methods %in% c("wald", "newcombe")
  labels[corrected] <- paste0(labels[corrected], "-cc")
  result <- new_ibist_ci(
    method = labels, estimate = c("proportion difference" = p[1L] - p[2L]),
    intervals = intervals, conf.level = conf.level,
    parameter = "difference of independent proportions",
    data.name = paste0("x = (", paste(counts$x, collapse = ", "),
                       "), n = (", paste(counts$n, collapse = ", "), ")"),
    null.value = c("proportion difference" = 0)
  )
  result$alternative <- settings$alternative
  result$correct <- correct
  result$sample.proportions <- stats::setNames(p, c("prop 1", "prop 2"))
  result
}

independent_binomial_counts <- function(x, n = NULL) {
  if (is.matrix(x) || inherits(x, "table")) {
    if (!is.null(n)) stop("Omit 'n' when 'x' is a 2 x 2 table.")
    tab <- as.matrix(x)
    if (!is.numeric(tab) || !identical(dim(tab), c(2L, 2L)) ||
          any(!is.finite(tab)) || any(tab < 0) || any(tab != floor(tab))) {
      stop("'x' must be a 2 x 2 table of non-negative integer counts.")
    }
    storage.mode(tab) <- "double"
    n <- rowSums(tab)
    x <- tab[, 1L]
  }
  if (length(x) != 2L || length(n) != 2L) {
    stop("'x' and 'n' must each contain two counts.")
  }
  validate_binomial_counts(x, n)
}

independent_score_limits <- function(x, n, z, tol, variance.factor) {
  estimate <- x[1L] / n[1L] - x[2L] / n[2L]
  score <- function(delta) {
    if (delta == -1) return(if (estimate == -1) 0 else Inf)
    if (delta == 1) return(if (estimate == 1) 0 else -Inf)
    range <- c(max(0, -delta), min(1, 1 - delta))
    derivative <- function(q) {
      p <- pmax(0, pmin(1, c(q + delta, q)))
      success <- x > 0
      failure <- x < n
      sum(x[success] / p[success]) -
        sum((n - x)[failure] / (1 - p[failure]))
    }
    q <- if (derivative(range[1L]) <= 0) range[1L] else {
      if (derivative(range[2L]) >= 0) range[2L] else {
        stats::uniroot(derivative, range, tol = tol / 100)$root
      }
    }
    p <- c(q + delta, q)
    variance <- sum(p * (1 - p) / n) * variance.factor
    if (variance <= 0) {
      return(if (estimate == delta) 0 else sign(estimate - delta) * Inf)
    }
    (estimate - delta) / sqrt(variance)
  }
  limit <- function(target) {
    if (score(-1) <= target) return(-1)
    if (score(1) >= target) return(1)
    stats::uniroot(function(d) score(d) - target, c(-1, 1), tol = tol)$root
  }
  c(limit(z), limit(-z))
}
