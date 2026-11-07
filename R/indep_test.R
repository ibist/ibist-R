#' Tests of Independence in a Contingency Table
#'
#' Tests independence of the row and column variables using Pearson's
#' chi-squared statistic or the likelihood-ratio G-squared statistic.
#'
#' @param x A two-dimensional matrix or table of non-negative integer counts,
#'   with at least two rows and two columns and positive marginal totals.
#' @param method The test statistic: \code{"pearson"} (default) or
#'   \code{"lr"} (likelihood ratio).
#' @param correct Logical. Apply Yates' continuity correction to Pearson's
#'   test for a 2 x 2 table. Default is \code{FALSE}. Setting this to
#'   \code{TRUE} for another method or table size is an error.
#'
#' @details
#' Under independence, the expected count in cell (i,j) is
#' \eqn{E_{ij}=n_{i+}n_{+j}/n}. Pearson's statistic is
#' \deqn{X^2=\sum_{ij}(n_{ij}-E_{ij})^2/E_{ij}.}
#' The likelihood-ratio statistic is
#' \deqn{G^2=2\sum_{ij}n_{ij}\log(n_{ij}/E_{ij}),}
#' where cells with zero observed counts contribute zero. Both tests use
#' the upper tail of a chi-squared distribution with
#' \eqn{(R-1)(C-1)} degrees of freedom.
#'
#' The Pearson method wraps \code{\link[stats]{chisq.test}}. Unlike that
#' function, \code{indep.test()} defaults to no continuity correction and
#' accepts only contingency tables with integer counts and positive row
#' and column totals. Empty margins are rejected rather than dropped.
#' With \code{correct = TRUE}, the Pearson method uses the same capped
#' Yates correction as \code{chisq.test()}.
#'
#' The likelihood-ratio statistic is computed directly. Expected counts and
#' Pearson residual diagnostics for both methods come from
#' \code{chisq.test()}. Both methods retain its warning when any expected
#' count is below 5; small expected counts can make the asymptotic p-value
#' unreliable. Exact tests, simulated p-values, goodness-of-fit tests, and
#' continuity corrections for the likelihood-ratio test are not implemented.
#'
#' @return An object of class \code{"htest"} containing:
#' \describe{
#'   \item{statistic}{The named test statistic, \code{X-squared} or
#'     \code{G-squared}.}
#'   \item{parameter}{Degrees of freedom, named \code{df}.}
#'   \item{p.value}{The asymptotic upper-tail p-value.}
#'   \item{method, data.name}{Descriptions of the test and input.}
#'   \item{observed, expected}{Observed and expected count matrices,
#'     retaining the input's dimension names.}
#'   \item{residuals, stdres}{Pearson and standardized Pearson residuals,
#'     respectively, for either test method. These are not deviance residuals.}
#' }
#'
#' @references
#' SAS Institute Inc. (2016). \emph{SAS/STAT 14.2 User's Guide}.
#' The FREQ Procedure, pp. 2759--2763: Chi-Square Tests and Statistics.
#' Cary, NC: SAS Institute Inc.
#' \url{https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf}
#'
#' @seealso \code{\link[stats]{chisq.test}},
#'   \code{\link[stats]{fisher.test}}
#'
#' @examples
#' screening <- matrix(c(1366, 14205, 383, 4983), nrow = 2, byrow = TRUE)
#' indep.test(screening)
#' indep.test(screening, method = "lr")
#' indep.test(screening, correct = TRUE)
#' indep.test(screening)$expected
#'
#' @export
indep.test <- function(x, method = c("pearson", "lr"), correct = FALSE) {
  method <- match.arg(method)
  tab <- as.matrix(x)
  if (length(dim(x)) != 2L || any(dim(tab) < 2L) || !is.numeric(tab) ||
      anyNA(tab) || any(!is.finite(tab)) || any(tab < 0) ||
      any(tab != floor(tab))) {
    stop("'x' must be a two-dimensional table of non-negative integer counts,",
         " with at least two rows and two columns.", call. = FALSE)
  }
  if (!is.finite(sum(tab)) || any(rowSums(tab) == 0) ||
      any(colSums(tab) == 0)) {
    stop("'x' must have a finite total and positive row and column totals.",
         call. = FALSE)
  }
  if (!is.logical(correct) || length(correct) != 1L || is.na(correct)) {
    stop("'correct' must be TRUE or FALSE.", call. = FALSE)
  }
  if (correct && (method != "pearson" || any(dim(tab) != 2L))) {
    stop("'correct = TRUE' is only available for Pearson's test ",
         "of a 2 x 2 table.",
         call. = FALSE)
  }

  result <- stats::chisq.test(tab, correct = correct)
  result$data.name <- paste(deparse(substitute(x)), collapse = "\n")
  if (method == "lr") {
    positive <- tab > 0
    statistic <- 2 * sum(tab[positive] *
                          log(tab[positive] / result$expected[positive]))
    result$statistic <- c("G-squared" = max(0, statistic))
    result$p.value <- stats::pchisq(result$statistic, result$parameter,
                                  lower.tail = FALSE)
    result$method <- "Likelihood-ratio G-squared test of independence"
  }
  result
}
