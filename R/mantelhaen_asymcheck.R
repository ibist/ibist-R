#' Check the Mantel-Haenszel Asymptotic Condition
#'
#' Checks whether the large-sample condition for the Mantel-Haenszel test is
#' satisfied for a collection of stratified 2 x 2 tables.
#'
#' @param x A list of 2 x 2 tables, or a 2 x 2 x K array of non-negative
#'   integer counts. Rows are the two groups and columns are the two outcome
#'   categories.
#'
#' @details
#' Conditional on the margins of stratum \eqn{s}, the upper-left cell count
#' \eqn{N_{s11}} has hypergeometric variance
#' \deqn{
#' V(N_{s11}) =
#' \frac{n_{s1.}n_{s2.}n_{s.1}n_{s.2}}{n_s^2(n_s-1)}.
#' }
#' The chi-squared approximation used by the Mantel-Haenszel test is considered
#' adequate when
#' \deqn{\sum_s V(N_{s11}) \ge 5.}
#'
#' This diagnostic applies to the asymptotic test. It is not required for an
#' exact conditional test such as \code{\link{or.common.test}}.
#'
#' @return An object of class \code{"ibist_mantelhaen_asymcheck"} containing:
#' \item{variance}{the hypergeometric variance for each stratum.}
#' \item{variance.sum}{the sum of the stratum-specific variances.}
#' \item{condition.met}{a logical value indicating whether the sum is at least
#'   5.}
#' \item{data.name}{a character string describing the data.}
#'
#' @examples
#' strat <- array(
#'   c(88, 469, 48, 653, 27, 480, 34, 1324),
#'   dim = c(2, 2, 2)
#' )
#' mantelhaen.asymcheck(strat)
#'
#' @seealso \code{\link[stats]{mantelhaen.test}},
#'   \code{\link{or.common.test}}
#'
#' @export
mantelhaen.asymcheck <- function(x) {
  tables <- validate_stratified_tables(x)
  variance <- vapply(tables, mantelhaen_stratum_variance, numeric(1))
  names(variance) <- paste("Stratum", seq_along(variance))
  variance_sum <- sum(variance)

  structure(
    list(
      variance = variance,
      variance.sum = variance_sum,
      condition.met = variance_sum >= 5,
      data.name = deparse(substitute(x))
    ),
    class = "ibist_mantelhaen_asymcheck"
  )
}

mantelhaen_stratum_variance <- function(tab) {
  n <- sum(tab)
  if (n <= 1) {
    return(0)
  }
  prod(c(rowSums(tab), colSums(tab))) / ((n - 1) * n^2)
}

#' @export
print.ibist_mantelhaen_asymcheck <- function(
    x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nMantel-Haenszel asymptotic condition check\n\n")
  cat("Data: ", x$data.name, "\n\n", sep = "")
  cat("Hypergeometric variance by stratum:\n")
  print(x$variance, digits = digits)
  cat(
    "\nSum of variances: ",
    format(x$variance.sum, digits = digits),
    "\nCondition (sum >= 5): ",
    if (x$condition.met) "met" else "not met",
    "\n",
    sep = ""
  )
  invisible(x)
}
