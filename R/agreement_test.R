#' Tests of Agreement Between Two Raters
#'
#' Computes Cohen's kappa or weighted kappa, a Wald confidence interval, and
#' an asymptotic test of zero agreement beyond chance. This entry point
#' preserves the inference and defaults of \code{cohen.kappa.test()}.
#'
#' @param x A square matrix or table of non-negative integer paired-rating
#'   counts, with at least two categories and a positive total. Rows and
#'   columns must represent the same categories in the same order. Counts
#'   and their total must not exceed \code{2^53 - 1}.
#' @param weights Agreement weights: \code{"unweighted"} (default),
#'   \code{"linear"}, \code{"quadratic"}, or a numeric matrix. Named weights
#'   use category positions, not numeric labels. A custom matrix must match
#'   \code{x}, be symmetric, have ones on its diagonal, and contain finite
#'   values in [0, 1]. Empty categories are retained when assigning weights.
#' @param conf.level Confidence level strictly between zero and one;
#'   default \code{0.95}.
#' @param alternative Alternative hypothesis: \code{"two.sided"} (default),
#'   \code{"greater"}, or \code{"less"}. Also determines interval direction.
#' @param test.variance Variance for the test: \code{"fce"} (default) uses the
#'   Fleiss--Cohen--Everitt null variance; \code{"delta"} uses the observed
#'   multinomial estimation variance. The CI always uses estimation variance.
#'
#' @details
#' With observed weighted agreement \eqn{P_o} and chance-expected weighted
#' agreement \eqn{P_e} from independent row and column margins, kappa is
#' \eqn{(P_o-P_e)/(1-P_e)}. Unlike association, agreement requires matching
#' rating categories; high association need not imply high agreement.
#' Linear weights are \eqn{1-|i-j|/(K-1)} and quadratic Fleiss--Cohen weights
#' are \eqn{1-\{(i-j)/(K-1)\}^2}. Unweighted kappa uses identity weights.
#'
#' The CI uses the FCE estimation variance, algebraically equivalent to the
#' multinomial delta-method variance. The default test uses the distinct
#' variance under zero kappa; \code{test.variance = "delta"} substitutes the
#' observed estimation variance. Consequently, \code{null.std.error} records
#' the selected test standard error, even when it is an estimation standard
#' error rather than a null standard error.
#'
#' Two-sided Wald intervals are untruncated. One-sided finite limits are
#' restricted to [-1, 1], with the other endpoint at 1 or -1. Intervals need
#' not invert the default test. Undefined kappa, invalid variances, and a
#' nonpositive selected test variance are errors. Exact and multi-rater
#' inference are not implemented. See \code{cohen.kappa.test()} for formulas.
#'
#' @return An \code{htest} with \code{estimate}, \code{null.value = 0}, signed
#'   Z \code{statistic}, \code{p.value}, \code{conf.int}, \code{alternative},
#'   \code{method}, and \code{data.name}. Additional fields are
#'   \code{std.error} (CI standard error), \code{null.std.error} (selected test
#'   standard error), \code{test.variance}, \code{weights} (matrix),
#'   \code{weight.method}, and \code{n}. The estimate is named \code{kappa}
#'   for unweighted inference and \code{weighted kappa} otherwise.
#'
#' @references
#' Fleiss, J. L., Cohen, J., and Everitt, B. S. (1969). Large sample standard
#' errors of kappa and weighted kappa. \emph{Psychological Bulletin}, 72(5),
#' 323--327. \doi{10.1037/h0028106}
#'
#' SAS Institute Inc. (2015). Tests and measures of agreement.
#' \emph{SAS/STAT 14.1 User's Guide}, The FREQ Procedure.
#' \url{https://support.sas.com/documentation/onlinedoc/stat/141/freq.pdf}
#'
#' @seealso \code{\link{cohen.kappa.test}}, \code{\link{assoc.test}}
#'
#' @examples
#' ratings <- matrix(c(9, 1, 2, 8), 2, byrow = TRUE)
#' agreement.test(ratings)
#' agreement.test(ratings, weights = "quadratic")
#' agreement.test(ratings, weights = "linear", test.variance = "delta")
#' agreement.test(ratings, alternative = "greater")
#'
#' @export
agreement.test <- function(x, weights = "unweighted", conf.level = 0.95,
                          alternative = "two.sided", test.variance = "fce") {
  tab <- as.matrix(x)
  if (is.numeric(tab) && all(is.finite(tab)) &&
      (any(tab > 2^53 - 1) || sum(tab) > 2^53 - 1)) {
    stop("Counts and their total must not exceed 2^53 - 1.")
  }
  result <- cohen.kappa.test(x, weights = weights, conf.level = conf.level,
                            alternative = alternative,
                            test.variance = test.variance)
  result$data.name <- paste(deparse(substitute(x)), collapse = "\n")
  result
}
