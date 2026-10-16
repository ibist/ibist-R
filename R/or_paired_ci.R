#' Confidence Intervals for a Matched-Pair Odds Ratio
#'
#' Computes confidence intervals for the conditional odds ratio from a
#' 2 x 2 table of matched pairs.
#'
#' @param x A 2 x 2 table of non-negative integer pair counts. Rows describe
#'   the first member of each pair and columns the second, with the positive
#'   category first in both dimensions. Thus \code{x[1, 2]} is the number
#'   positive only for the first member and \code{x[2, 1]} the number positive
#'   only for the second. In a matched case-control study, use cases in rows
#'   and controls in columns, with exposure present first.
#' @param conf.level Confidence level for the intervals.
#' @param method One or more of \code{"wald"} (the default), \code{"wilson"},
#'   \code{"clopper-pearson"}, \code{"midp"}, or \code{"blaker"}.
#'
#' @details
#' Write \eqn{b = x_{12}} and \eqn{c = x_{21}}. The matched-pair odds ratio
#' is estimated by \eqn{b/c}, not by the cross-product ratio of the pair-count
#' table. Concordant pairs do not contribute information about this
#' conditional odds ratio. Transposing the table or reversing both binary
#' category orders reciprocates the odds ratio. Reversing only one category
#' order changes which pairs are discordant and is not the same analysis.
#'
#' The \code{"wald"} interval is
#' \eqn{\exp\{\log(b/c) \mathbin{\pm} z_{1-\alpha/2}\sqrt{1/b+1/c}\}}.
#' It equals the asymptotic Mantel--Haenszel interval when each pair is
#' represented as an individual stratum. It requires both discordant counts
#' to be positive; zero concordant counts are allowed.
#'
#' Conditional on \eqn{b+c}, \eqn{B} is binomial with success probability
#' \eqn{p = \theta/(1+\theta)}. The remaining methods transform binomial
#' probability limits \eqn{(L,U)} to \eqn{(L/(1-L),U/(1-U))}.
#' The Wilson method inverts the score test, Clopper--Pearson inverts
#' equal-tailed exact tests, and mid-p subtracts half the observed mass from
#' each inclusive tail. Blaker's method takes the smallest interval
#' containing the acceptance set of its exact two-sided test. Its
#' acceptability statistic is the smaller inclusive binomial tail, and its
#' p-value sums probabilities of outcomes with acceptability no larger
#' than the observed value.
#'
#' Clopper--Pearson and Blaker intervals have at least nominal conditional
#' coverage; Wald, Wilson, and mid-p intervals have no such guarantee.
#' With only one nonzero discordant count the estimate is zero or infinity,
#' and transformed intervals remain available. With no discordant pairs,
#' the odds ratio is unidentified: the estimate is \code{NA}, transformed
#' intervals are \code{c(0, Inf)}, and a warning is issued.
#'
#' This function estimates a conditional matched-pair odds ratio, not a
#' marginal odds ratio or a difference of paired proportions. For the latter,
#' use \code{prop.paired.ci()}.
#'
#' @return An object of class \code{"ibist_ci"} with a matched-pair odds-ratio
#'   estimate and one confidence-limit row per requested method.
#'
#' @references
#' Fagerland, M. W., Lydersen, S., and Laake, P. (2014). Recommended tests
#' and confidence intervals for paired binomial proportions.
#' \emph{Statistics in Medicine}, 33, 2850--2875.
#' \doi{10.1002/sim.6148}
#'
#' Blaker, H. (2000). Confidence curves and improved exact confidence
#' intervals for discrete distributions. \emph{Canadian Journal of
#' Statistics}, 28, 783--798. \doi{10.2307/3315916}
#'
#' @examples
#' twins <- matrix(c(51, 30, 10, 41), nrow = 2, byrow = TRUE)
#' or.paired.ci(twins)
#' or.paired.ci(twins,
#'              method = c("wilson", "clopper-pearson", "midp", "blaker"))
#' or.paired.ci(matrix(c(8, 4, 0, 6), 2, byrow = TRUE), method = "wilson")
#'
#' @export
or.paired.ci <- function(x, conf.level = 0.95,
                         method = c("wald", "wilson", "clopper-pearson",
                                    "midp", "blaker")) {
  tab <- validate_2x2_table(x)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }
  methods <- if (missing(method)) method[1L] else {
    match.arg(method, several.ok = TRUE)
  }
  b <- as.numeric(tab[1L, 2L])
  c <- as.numeric(tab[2L, 1L])
  if ("wald" %in% methods && (b == 0 || c == 0)) {
    stop("Wald interval requires both discordant counts to be positive.")
  }
  if (b + c == 0) {
    warning("No discordant pairs: the matched-pair odds ratio is unidentified.")
  }
  estimate <- if (b + c == 0) NA_real_ else b / c
  intervals <- lapply(methods, function(method) {
    if (b + c == 0) return(c(0, Inf))
    if (method == "wald") {
      z <- stats::qnorm((1 + conf.level) / 2)
      return(exp(log(b / c) + c(-1, 1) * z * sqrt(1 / b + 1 / c)))
    }
    probability <- paired_or_binomial_ci(b, c, conf.level, method)
    probability / (1 - probability)
  })
  new_ibist_ci(method = methods,
               estimate = c("matched-pair odds ratio" = estimate),
               intervals = intervals, conf.level = conf.level,
               parameter = "matched-pair odds ratio",
               data.name = deparse(substitute(x)),
               null.value = c("matched-pair odds ratio" = 1))
}

paired_or_binomial_ci <- function(b, c, conf.level, method) {
  n <- b + c
  alpha <- 1 - conf.level
  if (method == "wilson") {
    z <- stats::qnorm(1 - alpha / 2)
    center <- (b + z ^ 2 / 2) / (n + z ^ 2)
    half <- z * sqrt(b * c / n + z ^ 2 / 4) / (n + z ^ 2)
    return(c(if (b == 0) 0 else center - half,
             if (c == 0) 1 else center + half))
  }
  if (method == "clopper-pearson") {
    interval <- stats::binom.test(b, n, conf.level = conf.level)$conf.int
    return(as.numeric(interval))
  }
  lower <- function(x) {
    if (x == 0) return(0)
    if (method == "blaker") return(paired_or_blaker_lower(x, n, alpha))
    stats::uniroot(function(p) {
      stats::pbinom(x - 1, n, p, lower.tail = FALSE) -
        stats::dbinom(x, n, p) / 2 - alpha / 2
    }, c(0, 1), tol = 1e-12)$root
  }
  c(lower(b), 1 - lower(c))
}

paired_or_blaker_lower <- function(x, n, alpha) {
  ## The lower limit is between equal-tail and one-tail exact limits.
  low <- stats::qbeta(alpha / 2, x, n - x + 1)
  high <- stats::qbeta(alpha, x, n - x + 1)
  right <- function(p) stats::pbinom(x - 1, n, p, lower.tail = FALSE)
  best <- high
  for (j in seq.int(0, x - 1)) {
    left <- function(p) stats::pbinom(j, n, p)
    ## The opposite tail contributes only when it is no larger than right.
    if (left(high) > right(high)) next
    if (left(high) + right(high) <= alpha) next
    start <- if (left(low) <= right(low)) low else {
      stats::uniroot(function(p) left(p) - right(p), c(low, high),
                     tol = 1e-12)$root
    }
    ## Check each tail-ordering boundary, including discontinuous crossings.
    candidate <- if (left(start) + right(start) >= alpha) start else {
      stats::uniroot(function(p) left(p) + right(p) - alpha,
                     c(start, high), tol = 1e-12)$root
    }
    best <- min(best, candidate)
  }
  best
}
