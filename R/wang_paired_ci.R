#' @useDynLib ibist, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Wang Exact Confidence Interval for Paired Proportions
#'
#' Computes Wang's exact inductive confidence interval for the paired risk
#' difference \eqn{p_{10} - p_{01}} from a paired binary table.
#'
#' @param n10 Number of pairs with success under treatment and failure under
#'   control.
#' @param t Number of concordant pairs, \eqn{n_{11} + n_{00}}.
#' @param n01 Number of pairs with failure under treatment and success under
#'   control.
#' @param conf.level Confidence level.
#' @param CItype Type of interval: \code{"Two.sided"}, \code{"Lower"}, or
#'   \code{"Upper"}.
#' @param precision Numerical precision for confidence limits.
#' @param grid.one Number of grid points in the first nuisance-parameter search.
#' @param grid.two Number of grid points in the second nuisance-parameter search.
#'
#' @return A list with elements \code{conf.level}, \code{CItype},
#'   \code{estimate}, and \code{ExactCI}.
#'
#' @references
#' Wang, W. (2012). An inductive order construction for the difference of two
#' dependent proportions. \emph{Statistics & Probability Letters}, 82,
#' 1623--1628.
#'
#' @examples
#' wang.paired.ci(3, 1, 0, CItype = "Lower")
#' wang.paired.ci(3, 1, 0, conf.level = 0.9)
#'
#' @noRd
wang.paired.ci <- function(n10, t, n01, conf.level = 0.95,
                           CItype = c("Two.sided", "Lower", "Upper"),
                           precision = 0.00001,
                           grid.one = 30, grid.two = 20) {
  if (length(n10) != 1L || length(t) != 1L || length(n01) != 1L) {
    stop("'n10', 't', and 'n01' must be single counts", call. = FALSE)
  }
  if (any(!is.finite(c(n10, t, n01))) || any(c(n10, t, n01) < 0) ||
      any(abs(c(n10, t, n01) - round(c(n10, t, n01))) > 0)) {
    stop("'n10', 't', and 'n01' must be nonnegative integers", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number in (0, 1)", call. = FALSE)
  }
  if (!is.numeric(precision) || length(precision) != 1L ||
      !is.finite(precision) || precision <= 0) {
    stop("'precision' must be a positive number", call. = FALSE)
  }
  if (any(!is.finite(c(grid.one, grid.two))) ||
      any(c(grid.one, grid.two) < 2) ||
      any(abs(c(grid.one, grid.two) - round(c(grid.one, grid.two))) > 0)) {
    stop("'grid.one' and 'grid.two' must be integers at least 2", call. = FALSE)
  }

  CItype <- match.arg(CItype)
  ci <- wang_paired_ci_cpp(
    as.integer(n10),
    as.integer(t),
    as.integer(n01),
    conf.level,
    CItype,
    precision,
    as.integer(grid.one),
    as.integer(grid.two)
  )

  list(
    conf.level = conf.level,
    CItype = CItype,
    estimate = unname(ci[["estimate"]]),
    ExactCI = unname(c(ci[["lower"]], ci[["upper"]]))
  )
}

#' Wang Interval-Inversion Test for Enumerated Paired Tables
#'
#' Evaluates whether Wang's exact paired confidence interval excludes zero for
#' each row of an enumerated paired-table data frame.
#'
#' @param tabs A data frame with columns \code{n11}, \code{n10}, \code{n01},
#'   and \code{n00}.
#' @param alpha Test size; the confidence level is \eqn{1 - \alpha}.
#' @inheritParams wang.paired.ci
#'
#' @return A logical vector indicating whether each row rejects
#'   \eqn{H_0: p_{10} - p_{01} = 0}.
#'
#' @examples
#' tabs <- data.frame(n11 = c(0, 1), n10 = c(3, 1), n01 = c(0, 1), n00 = c(1, 1))
#' wang.paired.reject(tabs, alpha = 0.05, precision = 0.0001)
#'
#' @noRd
wang.paired.reject <- function(tabs, alpha = 0.05, precision = 0.00001,
                               grid.one = 30, grid.two = 20) {
  if (!is.data.frame(tabs) ||
      !all(c("n11", "n10", "n01", "n00") %in% names(tabs))) {
    stop("'tabs' must contain n11, n10, n01, and n00 columns", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1)", call. = FALSE)
  }

  counts <- as.matrix(tabs[, c("n11", "n10", "n01", "n00")])
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(abs(counts - round(counts)) > 0)) {
    stop("table counts must be nonnegative integers", call. = FALSE)
  }

  N <- rowSums(counts)
  out <- logical(nrow(tabs))
  one_sided_level <- 1 - alpha / 2

  for (n in unique(N)) {
    rows <- which(N == n)
    if (n == 0) next

    lower_by_code <- wang_paired_lower_all_cpp(
      as.integer(n),
      one_sided_level,
      precision,
      as.integer(grid.one),
      as.integer(grid.two)
    )

    t <- tabs$n11[rows] + tabs$n00[rows]
    lower_code <- tabs$n10[rows] * (n + 2) + t + 1
    upper_code <- tabs$n01[rows] * (n + 2) + t + 1
    lower <- lower_by_code[lower_code]
    upper <- -lower_by_code[upper_code]
    out[rows] <- lower > 0 | upper < 0
  }

  out
}
