#' Confidence Intervals for a One-Sample Poisson Rate
#'
#' Computes confidence intervals for a Poisson rate parameter
#' \eqn{\lambda = x / T} using several methods.
#'
#' @param x Non-negative integer. Observed number of events.
#' @param exposure Positive numeric. Exposure time (or total time at risk).
#' @param side Interval form: \code{"two.sided"} (default),
#'   \code{"lower"}, or \code{"upper"}. An omitted (\code{NULL}) value
#'   selects two-sided limits. One-sided bounds use the full error probability.
#' @param alternative Compatibility alias: \code{"two.sided"},
#'   \code{"greater"} (lower bound), or \code{"less"} (upper bound).
#'   Conflicting \code{side} and \code{alternative} values are errors.
#' @param conf.level Confidence level. Default is 0.95.
#' @param method Method or methods for confidence interval. One or more of:
#'   \itemize{
#'     \item \code{"exact"}: Exact (Garwood) interval
#'     \item \code{"score"}: Score interval (inversion of score test)
#'     \item \code{"wh"}: Wilson-Hilferty interval
#'     \item \code{"wald"}: Wald interval (normal approximation)
#'     \item \code{"log"}: Log-transformed interval
#'   }
#' @param correct Logical. Apply continuity correction for
#'   \code{"score"}, \code{"wh"}, and \code{"wald"} methods. Default is TRUE.
#' @param ... Reserved for future extensions.
#'
#' @details
#' This function provides several confidence intervals for the Poisson
#' rate \eqn{\lambda = x / T}.
#'
#' \strong{Exact (Garwood):}
#' Based on inversion of the Poisson test using chi-square quantiles:
#' \deqn{
#'   \left[
#'     \frac{1}{2T} \chi^2_{2x, \alpha/2},
#'     \frac{1}{2T} \chi^2_{2(x+1), 1-\alpha/2}
#'   \right].
#' }
#'
#' \strong{Score:}
#' Obtained by inverting the score test. Computed from solving quadratic
#' equations for boundaries.
#'
#' \strong{Wilson-Hilferty:}
#' The Wilson–Hilferty (WH) interval is based on a cube-root transformation
#' of a chi-square approximation, which reduces skewness and yields an
#' approximately normal pivot, leading to the cubic form of the limits.
#' The lower limit is
#' \deqn{
#' \frac{x}{T}
#' \left(1 - \frac{1}{9x} + \frac{z_{\alpha/2}}{3\sqrt{x}}\right)^3
#' }
#' and the upper limit is
#' \deqn{
#' \frac{x + 1}{T}
#' \left(1 - \frac{1}{9(x + 1)} +
#' \frac{z_{1-\alpha/2}}{3\sqrt{x + 1}}\right)^3.
#' }
#'
#' \strong{Wald:}
#' For \eqn{x > 0},
#' \deqn{\hat{\lambda} \pm z_{\alpha/2} \sqrt{\hat{\lambda}/T}.}
#' For \eqn{x = 0}, lower limit 0, upper limit \eqn{-\log(\alpha/2)/T}.
#'
#' \strong{Log:}
#' For \eqn{x > 0}, applies normal approximation to \eqn{\log(\lambda)}
#' and transforms back. For \eqn{x = 0}, lower limit 0, upper limit
#' \eqn{-\log(\alpha/2)/T}.
#'
#' When \code{correct = TRUE}, continuity correction is applied on the count
#' scale for methods that support it.
#' 
#' @param asymcheck Logical; attach shared approximation diagnostics in
#'   \code{diagnostics}, default \code{FALSE}. See \code{\link{asymcheck}}.
#' @return An object of class \code{"ibist_ci"} containing the common estimate
#'   and a confidence-limit matrix with one row per method. The retained
#'   \code{inputs} data frame and \code{correct} permit later post-checks.
#'   Checks follow the book's Chapter 21 observed-count guidance: Wald
#'   \eqn{x\ge50}, log-Wald and uncorrected score \eqn{x\ge10}. WH is
#'   recommended even for small counts; no binary cutoff is assigned to WH
#'   or corrected score. These guidelines are not coverage guarantees.
#'
#' @examples
#' rate.1s.ci(5, 10)
#' rate.1s.ci(5, 10, method = "score")
#' rate.1s.ci(5, 10, method = c("exact", "score", "log"))
#' rate.1s.ci(0, 10, method = "exact")
#'
#' @importFrom stats qchisq

#' @export
rate.1s.ci <- function(
  x, exposure = 1.0,
  conf.level = 0.95,
  method = c("exact", "score", "wh", "wald", "log"),
  correct = TRUE,
  side = NULL, alternative = NULL, ..., asymcheck = FALSE
) {
  .asymcheck_flag(asymcheck, "asymcheck")
  orientation <- resolve_ci_side(side, alternative)
  side <- orientation$side

  if (length(x) != 1 || x < 0 || x != floor(x)) {
    stop("x must be a non-negative integer.")
  }

  if (length(exposure) != 1 || exposure <= 0) {
    stop("exposure must be positive.")
  }

  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be in (0, 1).")
  }

  methods <- if (missing(method)) method[1L] else {
    match.arg(method, several.ok = TRUE)
  }

  ci_methods <- list(
    exact = ci_exact,
    score = ci_score,
    wald  = ci_wald,
    log   = ci_log,
    wh    = ci_wh
  )

  intervals <- lapply(
    methods,
    function(method) {
      ci_methods[[method]](x, exposure, conf.level, correct,
                           alternative = orientation$alternative, ...)
    }
  )

  result <- new_ibist_ci(
    method = methods,
    estimate = c(rate = x / exposure),
    intervals = intervals,
    conf.level = conf.level,
    side = side,
    parameter = "Poisson rate",
    data.name = paste0("x = ", x, ", exposure = ", exposure)
  )
  result$inputs <- data.frame(x = x, exposure = exposure)
  result$correct <- correct
  if (asymcheck) result$diagnostics <- .asymcheck_rate(list(
    x = x, exposure = exposure, method = methods, conf.level = conf.level,
    correct = correct, side = side, alternative = NULL))
  result
}

ci_exact <- function(x, exposure, conf.level, correct = TRUE,
                     alternative = "two.sided", ...) {
  side <- resolve_ci_side(alternative = alternative)$side
  tail <- ci_tail(conf.level, side)
  lower <- if (x == 0) 0 else stats::qchisq(tail, 2 * x) /
    (2 * exposure)
  upper <- qchisq(1 - tail, 2 * (x + 1)) / (2 * exposure)
  ci_bounds(c(lower, upper), side)
}

ci_score <- function(x, exposure, conf.level, correct = TRUE,
                     alternative = "two.sided", ...) {
  alpha <- 1 - conf.level
  z <- qnorm(if (alternative == "two.sided") 1 - alpha / 2 else conf.level)

  if (correct) {
    ## optional: adjust x → x ± 0.5
    x_lower <- max(0, x - 0.5)
    x_upper <- x + 0.5
  } else {
    x_lower <- x_upper <- x
  }

  lower <- (2 * x_lower + z^2 -
            z * sqrt(z^2 + 4 * x_lower)) / (2 * exposure)

  upper <- (2 * x_upper + z^2 +
            z * sqrt(z^2 + 4 * x_upper)) / (2 * exposure)

  switch(alternative,
         two.sided = c(lower, upper),
         greater = c(max(0, lower), Inf),
         less = c(0, upper))
}

ci_wh <- function(x, exposure, conf.level, correct = TRUE,
                   alternative = "two.sided", ...) {
  side <- resolve_ci_side(alternative = alternative)$side
  tail <- ci_tail(conf.level, side)
  zL <- qnorm(tail)
  zU <- qnorm(1 - tail)

  if (correct) {
    lower_x <- upper_x <- x + 0.5
  } else {
    lower_x <- x
    upper_x <- x + 1
  }

  lower <- if (lower_x <= 0) {
    0
  } else {
    lower_x * (1 - 1 / (9 * lower_x) + zL / (3 * sqrt(lower_x)))^3
  }

  upper <- upper_x * (1 - 1 / (9 * upper_x) + zU / (3 * sqrt(upper_x)))^3

  ci_bounds(c(lower, upper) / exposure, side)
}

ci_wald <- function(x, exposure, conf.level, correct = TRUE,
                   alternative = "two.sided", ...) {
  side <- resolve_ci_side(alternative = alternative)$side
  tail <- ci_tail(conf.level, side)
  z <- qnorm(1 - tail)

  if (x == 0 && !correct) {
    return(ci_bounds(c(0, -log(tail) / exposure), side))
  }

  # continuity correction offset
  d <- if (correct) 0.5 else 0

  lower_x <- max(0, x - d)
  upper_x <- x + d

  lower <- max(0, (lower_x - z * sqrt(lower_x)) / exposure)
  upper <- max(0, (upper_x + z * sqrt(upper_x)) / exposure)

  ci_bounds(c(lower, upper), side)
}

ci_log <- function(x, exposure, conf.level, correct = TRUE,
                   alternative = "two.sided", ...) {
  side <- resolve_ci_side(alternative = alternative)$side
  tail <- ci_tail(conf.level, side)
  z <- qnorm(1 - tail)

  if (x == 0) return(ci_bounds(c(0, -log(tail) / exposure), side))

  rate <- x / exposure
  se <- 1 / sqrt(x)

  lower <- rate * exp(-z * se)
  upper <- rate * exp(z * se)

  ci_bounds(c(lower, upper), side)
}
