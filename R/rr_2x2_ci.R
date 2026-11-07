#' Confidence Intervals for a Risk Ratio in a 2 x 2 Table
#'
#' Computes intervals for the first row's event probability divided by the
#' second row's event probability in two independent binomial samples.
#'
#' @param x A 2 x 2 table of non-negative integer counts. Rows are samples;
#'   the first column contains events and the second contains non-events.
#'   Both row totals must be positive.
#' @param side Interval form: \code{"two.sided"} (default),
#'   \code{"lower"}, or \code{"upper"}. An omitted (\code{NULL}) value
#'   selects two-sided limits. One-sided bounds use the full error probability.
#' @param alternative Compatibility alias: \code{"two.sided"},
#'   \code{"greater"} (lower bound), or \code{"less"} (upper bound).
#'   Conflicting \code{side} and \code{alternative} values are errors.
#' @param conf.level Confidence level, between 0 and 1.
#' @param method One or more of \code{"koopman"} (default) or \code{"katz"}.
#' @param tol Positive numerical tolerance for score inversion.
#'
#' @details
#' Katz's method constructs a Wald interval for the log risk ratio, using
#' variance \eqn{1/x_1-1/n_1+1/x_2-1/n_2}, and exponentiates its limits.
#' It requires positive event counts in both samples and does not add any
#' pseudo-counts.
#'
#' Koopman's method inverts the constrained Pearson score statistic. For a
#' candidate ratio \eqn{r}, the binomial likelihood is maximized under
#' \eqn{p_1=r p_2}. The squared statistic is the sum of the two constrained
#' Pearson contributions. No continuity correction or finite-sample variance
#' factor is applied. Inversion is performed on the log-ratio scale, using
#' row reversal for reciprocal symmetry and numerical stability.
#'
#' If only one event count is zero, the estimate is zero or infinity and
#' Koopman's interval has the corresponding boundary limit. If both event
#' counts are zero, the estimate is \code{NA} and the interval is
#' \eqn{[0,\infty)}: the score construction does not identify a ratio.
#' One-sided intervals have upper limit infinity or lower limit zero.
#'
#' @return An \code{ibist_ci} object containing the observed named
#'   \code{estimate}, a \code{conf.int} matrix with one row per method,
#'   \code{conf.level}, \code{method}, \code{side}, and the two
#'   \code{sample.proportions}.
#' @references
#' Katz, D., Baptista, J., Azen, S. P., and Pike, M. C. (1978). Obtaining
#' confidence intervals for the risk ratio in cohort studies.
#' \emph{Biometrics}, 34, 469--474.
#'
#' Koopman, P. A. R. (1984). Confidence intervals for the ratio of two binomial
#' proportions. \emph{Biometrics}, 40, 513--517.
#' @examples
#' vf <- matrix(c(38, 62, 267, 717), nrow = 2, byrow = TRUE)
#' rr.2x2.ci(vf)
#' rr.2x2.ci(vf, method = c("katz", "koopman"))
#' rr.2x2.ci(vf, side = "lower")
#' @export
rr.2x2.ci <- function(
  x, conf.level = 0.95, method = c("koopman", "katz"),
  alternative = NULL, tol = 1e-8, side = NULL
) {
  if (!is.matrix(x) && !inherits(x, "table")) {
    stop("'x' must be a 2 x 2 table of non-negative integer counts.")
  }
  counts <- independent_binomial_counts(x)
  orientation <- resolve_ci_side(side, alternative)
  alternative <- orientation$alternative
  settings <- validate_binomial_ci_settings(conf.level, alternative, tol)
  methods <- if (missing(method)) method[1L] else {
    match.arg(method, several.ok = TRUE)
  }
  p <- counts$x / counts$n
  estimate <- if (all(counts$x == 0)) NA_real_ else p[1L] / p[2L]
  z <- stats::qnorm(settings$tail, lower.tail = FALSE)
  intervals <- lapply(methods, function(method) {
    if (method == "katz") {
      if (any(counts$x == 0)) {
        stop("Katz intervals require positive event counts in both samples.")
      }
      se <- sqrt(sum(1 / counts$x - 1 / counts$n))
      limits <- exp(log(estimate) + c(-1, 1) * z * se)
    } else {
      limits <- c(koopman_ratio_lower(counts$x, counts$n, z, tol),
                  1 / koopman_ratio_lower(rev(counts$x), rev(counts$n), z, tol))
    }
    if (settings$alternative == "greater") limits[2L] <- Inf
    if (settings$alternative == "less") limits[1L] <- 0
    limits
  })
  result <- new_ibist_ci(
    method = methods, estimate = c("risk ratio" = estimate),
    intervals = intervals, conf.level = conf.level, parameter = "risk ratio",
    data.name = deparse(substitute(x)), null.value = c("risk ratio" = 1)
  )
  result$side <- orientation$side
  result$alternative <- settings$alternative
  result$sample.proportions <- stats::setNames(p, c("prop 1", "prop 2"))
  result
}

koopman_ratio_score <- function(log.ratio, x, n) {
  if (log.ratio > 0) {
    return(-koopman_ratio_score(-log.ratio, rev(x), rev(n)))
  }
  ratio <- exp(log.ratio)
  total <- sum(x)
  if (total == 0) return(0)
  b <- ratio * (n[1L] + x[2L]) + n[2L] + x[1L]
  discriminant <- max(0, b^2 - 4 * ratio * sum(n) * total)
  q <- min(1, 2 * total / (b + sqrt(discriminant)))
  p <- c(ratio * q, q)
  residual <- x - n * p
  variance <- n * p * (1 - p)
  terms <- numeric(2L)
  positive <- variance > 0
  terms[positive] <- residual[positive]^2 / variance[positive]
  terms[!positive & abs(residual) > 1e-10] <- Inf
  difference <- x[1L] / n[1L] - ratio * x[2L] / n[2L]
  if (difference == 0) return(0)
  sign(difference) * sqrt(sum(terms))
}

koopman_ratio_lower <- function(x, n, z, tol) {
  if (all(x == 0) || (x[1L] == 0 && z >= 0)) return(0)
  objective <- function(log.ratio) koopman_ratio_score(log.ratio, x, n) - z
  lower <- -1
  upper <- 1
  while (objective(lower) < 0 && lower > -700) lower <- max(-700, lower * 2)
  while (objective(upper) > 0 && upper < 700) upper <- min(700, upper * 2)
  if (objective(lower) <= 0) return(0)
  if (objective(upper) >= 0) return(Inf)
  exp(stats::uniroot(objective, c(lower, upper), tol = tol)$root)
}
