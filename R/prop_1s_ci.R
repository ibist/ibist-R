#' Confidence Intervals for a Binomial Proportion
#'
#' Computes intervals for a binomial probability, optionally comparing methods
#' or applying the same methods to several observed counts.
#'
#' @param x Non-negative integer count or vector of counts of successes.
#' @param n Positive integer number of trials, either a scalar or a vector
#'   of the same length as \code{x}.
#' @param side Interval form: \code{"two.sided"} (default),
#'   \code{"lower"}, or \code{"upper"}. An omitted (\code{NULL}) value
#'   selects two-sided limits. One-sided bounds use the full error probability.
#' @param alternative Compatibility alias: \code{"two.sided"},
#'   \code{"greater"} (lower bound), or \code{"less"} (upper bound).
#'   Conflicting \code{side} and \code{alternative} values are errors.
#' @param conf.level Confidence level, between 0 and 1.
#' @param method One or more of \code{"wilson"} (default), \code{"wald"},
#'   \code{"agresti-coull"}, \code{"clopper-pearson"}, or \code{"blaker"}.
#' @param correct Logical. Apply continuity correction to Wilson and Wald
#'   intervals. Other requested methods are unchanged.
#' @param midp Logical. Replace Clopper--Pearson tail probabilities with
#'   mid-p probabilities. Other requested methods are unchanged.
#' @param tol Positive numerical tolerance for interval inversion.
#' @param asymcheck Logical. Attach the book's approximation checks under
#'   \code{diagnostics}; default \code{FALSE}. See \code{\link{asymcheck}}.
#'
#' @details
#' Wilson intervals invert the binomial score statistic. Continuity correction
#' replaces the observed proportion by \eqn{(x-1/2)/n} for the lower limit and
#' \eqn{(x+1/2)/n} for the upper limit. Wald intervals use the observed binomial
#' standard error, adding \eqn{1/(2n)} to the margin when corrected.
#' The full half-count correction is used, independently of any null
#' probability; unlike \code{stats::prop.test()}, no null-dependent cap is
#' applied to the correction.
#' Agresti--Coull intervals use \eqn{n+z^2} trials and \eqn{x+z^2/2} successes.
#' All limits are restricted to the parameter space \eqn{[0,1]}.
#'
#' Clopper--Pearson intervals invert binomial tails; \code{midp = TRUE}
#' subtracts half the observed-count probability from each tail. Mid-p
#' intervals do not guarantee at least nominal coverage. Blaker's two-sided
#' interval is the smallest interval containing the acceptance set obtained
#' by ordering counts by their smaller inclusive tail probability. The
#' algorithm partitions the lower-tail search at changes in this ordering,
#' rather than assuming a continuous, monotone two-sided p-value. One-sided
#' Blaker bounds use the corresponding exact binomial tail.
#'
#' A correction flag that applies to none of the requested methods is an
#' error. The point estimate is always \eqn{x/n}, including when a method
#' uses an adjusted estimate internally.
#'
#' @return An \code{ibist_ci} object with named \code{estimate}, a
#'   \code{conf.int} matrix with columns \code{lower} and \code{upper},
#'   \code{conf.level}, \code{method}, and \code{side}. For vector
#'   inputs, rows are ordered by input and then method; \code{inputs} records
#'   the input index, count, and sample size corresponding to each row.
#'   Original counts and settings are retained for \code{asymcheck(result)}.
#'   Diagnostics do not change intervals or select a different method.
#'
#' @references
#' Agresti, A., and Coull, B. A. (1998). Approximate is better than
#' ``exact'' for interval estimation of binomial proportions.
#' \emph{The American Statistician}, 52, 119--126.
#'
#' Blaker, H. (2000). Confidence curves and improved exact confidence intervals
#' for discrete distributions. \emph{Canadian Journal of Statistics}, 28,
#' 783--798. See also the correction (2001), 29, 681.
#'
#' Brown, L. D., Cai, T. T., and DasGupta, A. (2001). Interval estimation for
#' a binomial proportion. \emph{Statistical Science}, 16, 101--133.
#'
#' @examples
#' prop.1s.ci(5, 40)
#' prop.1s.ci(5, 40, method = c("wilson", "clopper-pearson", "blaker"))
#' prop.1s.ci(5, 40, method = "clopper-pearson", midp = TRUE)
#' prop.1s.ci(5, 40, side = "lower")
#' prop.1s.ci(0:5, 40)
#' @export
prop.1s.ci <- function(
  x, n, conf.level = 0.95,
  method = c("wilson", "wald", "agresti-coull", "clopper-pearson", "blaker"),
  correct = FALSE, midp = FALSE,
  alternative = NULL, tol = 1e-8, side = NULL, asymcheck = FALSE
) {
  validate_ci_flag(asymcheck, "asymcheck")
  counts <- validate_binomial_counts(x, n)
  orientation <- resolve_ci_side(side, alternative)
  alternative <- orientation$alternative
  settings <- validate_binomial_ci_settings(conf.level, alternative, tol)
  alternative <- settings$alternative
  methods <- if (missing(method)) method[1L] else {
    match.arg(method, several.ok = TRUE)
  }
  validate_ci_flag(correct, "correct")
  validate_ci_flag(midp, "midp")
  if (correct && !any(methods %in% c("wald", "wilson"))) {
    stop("'correct' applies only to Wald and Wilson intervals.")
  }
  if (midp && !"clopper-pearson" %in% methods) {
    stop("'midp' applies only to Clopper-Pearson intervals.")
  }

  rows <- expand.grid(method = methods, input = seq_along(counts$x),
                      stringsAsFactors = FALSE)
  intervals <- lapply(seq_len(nrow(rows)), function(i) {
    j <- rows$input[i]
    binomial_ci_limits(counts$x[j], counts$n[j], rows$method[i],
                       settings$tail, correct, midp, tol, alternative)
  })
  labels <- rows$method
  labels[correct & labels %in% c("wald", "wilson")] <-
    paste0(labels[correct & labels %in% c("wald", "wilson")], "-cc")
  labels[midp & rows$method == "clopper-pearson"] <- "clopper-pearson-midp"
  if (length(counts$x) > 1L) {
    labels <- paste0(labels, "[", rows$input, "]")
  }
  estimate <- counts$x / counts$n
  names(estimate) <- if (length(estimate) == 1L) "proportion" else {
    paste0("proportion[", seq_along(estimate), "]")
  }
  result <- new_ibist_ci(
    method = labels, estimate = estimate, intervals = intervals,
    conf.level = conf.level, parameter = "binomial proportion",
    data.name = if (length(counts$x) == 1L) {
      paste0("x = ", counts$x, ", n = ", counts$n)
    } else paste(length(counts$x), "binomial samples")
  )
  result$side <- orientation$side
  result$alternative <- alternative
  result$correct <- correct
  result$midp <- midp
  result$inputs <- data.frame(input = rows$input, x = counts$x[rows$input],
                              n = counts$n[rows$input])
  context <- list(family = "p1", x = counts$x, n = counts$n, method = methods,
    conf.level = conf.level, correct = correct, midp = midp, tol = tol,
    side = orientation$side, alternative = NULL)
  attr(result, "asymcheck.context") <- context
  if (asymcheck) result$diagnostics <- .asymcheck_book_ci(context)
  result
}

validate_binomial_counts <- function(x, n) {
  valid <- function(v) {
    is.numeric(v) && length(v) > 0L && all(is.finite(v)) &&
      all(v >= 0) && all(v == floor(v))
  }
  if (!valid(x) || !valid(n) || any(n == 0)) {
    stop("'x' and 'n' must be integer counts, with 'x' >= 0 and 'n' > 0.")
  }
  if (!length(n) %in% c(1L, length(x))) {
    stop("'n' must be a scalar or have the same length as 'x'.")
  }
  n <- rep_len(as.double(n), length(x))
  x <- as.double(x)
  if (any(x > n)) stop("'x' must not exceed 'n'.")
  list(x = x, n = n)
}

validate_ci_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", name, "' must be TRUE or FALSE.")
  }
}

validate_binomial_ci_settings <- function(conf.level, alternative, tol) {
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
        tol <= 0 || tol >= 1) {
    stop("'tol' must be a single number between 0 and 1.")
  }
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  list(alternative = alternative,
       tail = (1 - conf.level) / if (alternative == "two.sided") 2 else 1)
}

binomial_ci_limits <- function(x, n, method, tail, correct, midp, tol,
                                alternative = "two.sided") {
  p <- x / n
  z <- stats::qnorm(tail, lower.tail = FALSE)
  if (method == "wilson") {
    offset <- if (correct) 0.5 / n else 0
    lower.p <- max(0, p - offset)
    upper.p <- min(1, p + offset)
    endpoint <- function(q, sign) {
      (q + z^2 / (2 * n) + sign * z *
         sqrt(q * (1 - q) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)
    }
    limits <- c(if (x == 0) 0 else endpoint(lower.p, -1),
                if (x == n) 1 else endpoint(upper.p, 1))
  } else if (method %in% c("wald", "agresti-coull")) {
    size <- if (method == "agresti-coull") n + z^2 else n
    center <- if (method == "agresti-coull") (x + z^2 / 2) / size else p
    margin <- z * sqrt(center * (1 - center) / size)
    if (correct && method == "wald") margin <- margin + 0.5 / n
    limits <- center + c(-1, 1) * margin
  } else if (method == "blaker" && alternative == "two.sided") {
    limits <- c(blaker_binomial_lower(x, n, 2 * tail, tol),
                1 - blaker_binomial_lower(n - x, n, 2 * tail, tol))
  } else if (midp && method == "clopper-pearson") {
    lower.tail <- function(p) {
      stats::pbinom(x, n, p, lower.tail = FALSE) +
        0.5 * stats::dbinom(x, n, p) - tail
    }
    upper.tail <- function(p) {
      stats::pbinom(x - 1, n, p) + 0.5 * stats::dbinom(x, n, p) - tail
    }
    limits <- c(if (x == 0) 0 else if (lower.tail(1) <= 0) 1 else {
      stats::uniroot(lower.tail, c(0, 1), tol = tol)$root
    }, if (x == n) 1 else if (upper.tail(0) <= 0) 0 else {
      stats::uniroot(upper.tail, c(0, 1), tol = tol)$root
    })
  } else {
    limits <- c(if (x == 0) 0 else stats::qbeta(tail, x, n - x + 1),
                if (x == n) 1 else {
                  stats::qbeta(tail, x + 1, n - x, lower.tail = FALSE)
                })
  }
  # A signed Wilson bound beyond a constant response hits the domain edge.
  # Snap it exactly before odds transformations amplify roundoff near 0/1.
  if (method == "wilson" && !correct && z < 0) {
    if (x == n) limits[1L] <- 1
    if (x == 0) limits[2L] <- 0
  }
  limits <- pmax(0, pmin(1, limits))
  if (alternative == "greater") limits[2L] <- 1
  if (alternative == "less") limits[1L] <- 0
  limits
}

blaker_binomial_lower <- function(x, n, alpha, tol) {
  if (x == 0) return(0)
  lower <- stats::qbeta(alpha / 2, x, n - x + 1)
  upper <- min(x / n, stats::qbeta(alpha, x, n - x + 1))
  right <- function(p) stats::pbinom(x - 1, n, p, lower.tail = FALSE)
  ## Ordering changes only where an opposite tail equals the observed tail.
  first <- max(0, stats::qbinom(right(lower), n, lower) - 1)
  last <- min(x - 1, stats::qbinom(right(upper), n, upper))
  breaks <- c(lower, upper)
  if (first <= last) {
    for (k in seq.int(first, last)) {
      crossing <- function(p) stats::pbinom(k, n, p) - right(p)
      if (crossing(lower) > 0 && crossing(upper) < 0) {
        breaks <- c(breaks, stats::uniroot(crossing, c(lower, upper),
                                          tol = tol / 10)$root)
      }
    }
  }
  breaks <- sort(unique(breaks))
  for (i in seq_len(length(breaks) - 1L)) {
    interval <- breaks[c(i, i + 1L)]
    middle <- mean(interval)
    k <- min(x - 1, stats::qbinom(right(middle), n, middle))
    if (stats::pbinom(k, n, middle) > right(middle)) k <- k - 1
    acceptability <- function(p) right(p) + stats::pbinom(k, n, p) - alpha
    if (acceptability(interval[1L]) >= -1e-14) return(interval[1L])
    if (acceptability(interval[2L]) >= 0) {
      return(stats::uniroot(acceptability, interval, tol = tol)$root)
    }
  }
  upper
}
