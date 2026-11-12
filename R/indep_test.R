#' Tests of Independence in a Contingency Table
#'
#' Tests independence of the row and column variables using Pearson's
#' chi-squared statistic, the likelihood-ratio G-squared statistic, or
#' Fisher's conditional test.
#'
#' @param x A two-dimensional matrix or table of non-negative integer counts,
#'   with at least two rows and two columns and positive marginal totals.
#' @param method The test: \code{"pearson"} (default), \code{"lr"}
#'   (likelihood ratio), or \code{"fisher"}.
#' @param correct Logical. Apply Yates' continuity correction to Pearson's
#'   test for a 2 x 2 table. Default is \code{FALSE}. Setting this to
#'   \code{TRUE} for another method or table size is an error.
#' @param asymcheck Logical; attach shared approximation diagnostics in
#'   \code{diagnostics}, default \code{FALSE}. See \code{\link{asymcheck}}.
#' @param ... Named options for Fisher's method, passed to
#'   \code{\link[stats]{fisher.test}}: \code{alternative}, \code{conf.int},
#'   \code{conf.level}, \code{workspace}, \code{simulate.p.value}, and
#'   \code{B}. Other names, unnamed options, and options supplied to another
#'   method are errors. Fisher's null odds ratio is fixed at one.
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
#' unreliable.
#'
#' Fisher's method wraps \code{fisher.test()} for 2 x 2 and larger tables,
#' conditioning on the row and column margins. It defaults to a two-sided
#' exact test of independence. For 2 x 2 tables, \code{alternative = "less"}
#' or \code{"greater"} tests an odds ratio below or above one, respectively.
#' The returned estimate is the conditional maximum likelihood odds ratio;
#' its confidence interval is for that odds ratio, not a proportion
#' difference. Fisher's method does not use a continuity correction and does
#' not issue the small-expected-count approximation warning.
#'
#' For larger tables, only a two-sided test is available. The exact network
#' algorithm may require a larger \code{workspace}; its errors are passed
#' through. Alternatively, \code{simulate.p.value = TRUE} requests a Monte
#' Carlo p-value with \code{B} replicates (default 2000). Set a random seed
#' before the call for reproducibility. Simulation applies only to tables
#' larger than 2 x 2; Fisher's 2 x 2 calculation remains exact.
#' Goodness-of-fit tests, hybrid Fisher approximations, and likelihood-ratio
#' continuity corrections are not implemented.
#'
#' @return An object of class \code{"htest"} containing:
#' \describe{
#'   \item{statistic}{The named test statistic, \code{X-squared} or
#'     \code{G-squared}, for the two asymptotic methods. Fisher's method
#'     does not return this field.}
#'   \item{parameter}{Degrees of freedom, named \code{df}, for the two
#'     asymptotic methods. Fisher's method does not return this field.}
#'   \item{p.value}{The asymptotic, exact, or Monte Carlo p-value.}
#'   \item{method, data.name}{Descriptions of the test and input.}
#'   \item{observed, expected}{Observed and expected count matrices,
#'     retaining the input's dimension names.}
#'   \item{residuals, stdres}{Pearson and standardized Pearson residuals,
#'     respectively, for the two asymptotic methods. These are not deviance
#'     residuals and are not returned by Fisher's method.}
#'   \item{estimate, null.value, alternative, conf.int}{For 2 x 2 Fisher
#'     tests, the odds-ratio inference fields from \code{fisher.test()}.
#'     The confidence interval is omitted when \code{conf.int = FALSE}.}
#' }
#'
#' @references
#' SAS Institute Inc. (2016). \emph{SAS/STAT 14.2 User's Guide}.
#' The FREQ Procedure, pp. 2759--2763: Chi-Square Tests and Statistics.
#' Cary, NC: SAS Institute Inc.
#' \url{https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf}
#'
#' Fisher, R. A. (1935). The logic of inductive inference.
#' \emph{Journal of the Royal Statistical Society}, 98(1), 39--82.
#' \doi{10.2307/2342435}
#'
#' Fisher, R. A. (1962). Confidence limits for a cross-product ratio.
#' \emph{Australian Journal of Statistics}, 4(1), 41.
#' \doi{10.1111/j.1467-842X.1962.tb00285.x}
#'
#' Mehta, C. R., and Patel, N. R. (1983). A network algorithm for performing
#' Fisher's exact test in r x c contingency tables.
#' \emph{Journal of the American Statistical Association}, 78(382), 427--434.
#' \doi{10.1080/01621459.1983.10477989}
#'
#' Clarkson, D. B., Fan, Y., and Joe, H. (1993). A remark on Algorithm 643:
#' FEXACT: An algorithm for performing Fisher's exact test in r x c
#' contingency tables. \emph{ACM Transactions on Mathematical Software},
#' 19(4), 484--488. \doi{10.1145/168173.168412}
#'
#' Patefield, W. M. (1981). Algorithm AS 159: An efficient method of
#' generating R x C tables with given row and column totals.
#' \emph{Applied Statistics}, 30(1), 91--97. \doi{10.2307/2346669}
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
#' small <- matrix(c(7, 77, 2, 90), nrow = 2, byrow = TRUE)
#' indep.test(small, method = "fisher")
#' indep.test(small, method = "fisher", alternative = "greater",
#'            conf.int = FALSE)
#'
#' @export
indep.test <- function(x, method = c("pearson", "lr", "fisher"),
                      correct = FALSE, ..., asymcheck = FALSE) {
  .asymcheck_flag(asymcheck, "asymcheck")
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

  fisher.args <- list(...)
  if (length(fisher.args)) {
    allowed <- c("alternative", "conf.int", "conf.level", "workspace",
                 "simulate.p.value", "B")
    option.names <- names(fisher.args)
    if (method != "fisher" || is.null(option.names) ||
        anyNA(option.names) || any(!option.names %in% allowed) ||
        anyDuplicated(option.names)) {
      stop("'...' accepts only named Fisher options: ",
           paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
  }
  data.name <- paste(deparse(substitute(x)), collapse = "\n")
  if (method == "fisher") {
    is.2x2 <- all(dim(tab) == 2L)
    if (!is.2x2 && !is.null(fisher.args$alternative) &&
        !identical(fisher.args$alternative, "two.sided")) {
      stop("One-sided Fisher alternatives require a 2 x 2 table.",
           call. = FALSE)
    }
    result <- do.call(stats::fisher.test, c(list(x = tab), fisher.args))
    result$data.name <- data.name
    result$observed <- tab
    result$expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
    dimnames(result$expected) <- dimnames(tab)
    if (asymcheck) result$diagnostics <- .asymcheck_indep(tab, method, correct)
    return(result)
  }

  result <- stats::chisq.test(tab, correct = correct)
  result$data.name <- data.name
  if (method == "lr") {
    positive <- tab > 0
    statistic <- 2 * sum(tab[positive] *
                          log(tab[positive] / result$expected[positive]))
    result$statistic <- c("G-squared" = max(0, statistic))
    result$p.value <- stats::pchisq(result$statistic, result$parameter,
                                  lower.tail = FALSE)
    result$method <- "Likelihood-ratio G-squared test of independence"
  }
  if (asymcheck) result$diagnostics <- .asymcheck_indep(tab, method, correct)
  result
}
