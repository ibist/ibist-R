#' Nominal Association Measures for a Contingency Table
#'
#' Computes Phi, Pearson's contingency coefficient, and Cramer's V without
#' using any ordering of the row or column categories.
#'
#' @param x A two-dimensional matrix or table of non-negative integer counts,
#'   with at least two rows and two columns and positive marginal totals.
#' @param measure One or more of \code{"phi"}, \code{"contingency"}, or
#'   \code{"cramer-v"}. By default, all three are returned.
#'
#' @details
#' Let \eqn{X^2} be the uncorrected Pearson statistic for independence,
#' \eqn{n} the total count, and \eqn{d=\min(R,C)}. The measures are
#' \deqn{\phi=\sqrt{X^2/n},\qquad
#' C=\sqrt{X^2/(X^2+n)},\qquad
#' V=\sqrt{X^2/\{n(d-1)\}}.}
#' Phi is non-negative, including for 2 x 2 tables; it is not the signed
#' binary correlation. Its maximum is \eqn{\sqrt{d-1}}, so it can exceed
#' one for larger tables. The contingency coefficient's maximum is
#' \eqn{\sqrt{(d-1)/d}}. Cramer's V is between zero and one and equals
#' Phi when \eqn{d=2}. These are the unadjusted definitions used by
#' \code{DescTools::Phi()}, \code{DescTools::ContCoef()}, and
#' \code{DescTools::CramerV()} with their default arguments.
#' The formulas are also given in SAS Institute Inc. (2016, p. 2764),
#' except that SAS retains the sign of Phi and V for 2 x 2 tables;
#' this function reports their absolute values to match DescTools.
#'
#' The Pearson statistic and expected counts are obtained from
#' \code{\link{indep.test}} with no continuity correction. Empty margins
#' are rejected rather than dropped. Sparse tables do not produce the
#' chi-squared approximation warning here because this function returns
#' descriptive point estimates, not asymptotic p-values. Use
#' \code{indep.test()} to test independence and assess that warning.
#' Confidence intervals, bias corrections, and adjustments for table size
#' beyond the defining normalization of Cramer's V are not implemented.
#'
#' @return An object of class \code{"ibist_nominal_assoc"} containing:
#' \describe{
#'   \item{estimate}{Named vector of the selected association measures.}
#'   \item{measure}{Selected measure names in the requested order, with
#'     duplicates removed.}
#'   \item{statistic}{Uncorrected Pearson statistic, named \code{X-squared}.}
#'   \item{n}{Total count.}
#'   \item{observed, expected}{Observed and expected count matrices.}
#'   \item{data.name}{Description of the input.}
#' }
#'
#' @references
#' SAS Institute Inc. (2016). \emph{SAS/STAT 14.2 User's Guide}.
#' The FREQ Procedure, p. 2764: Phi Coefficient, Contingency Coefficient,
#' and Cramer's V. Cary, NC: SAS Institute Inc.
#' \url{https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf}
#'
#' @seealso \code{\link{ordinal.assoc}}, \code{\link{indep.test}}
#'
#' @examples
#' tumor <- matrix(c(109, 68, 14, 72, 79, 37, 9, 42, 139),
#'                 nrow = 3, byrow = TRUE)
#' nominal.assoc(tumor)
#' nominal.assoc(tumor, measure = "cramer-v")
#' nominal.assoc(tumor)$estimate
#'
#' @export
nominal.assoc <- function(x, measure = c("phi", "contingency", "cramer-v")) {
  all.measures <- c("phi", "contingency", "cramer-v")
  measures <- if (missing(measure)) {
    all.measures
  } else {
    unique(match.arg(measure, choices = all.measures, several.ok = TRUE))
  }

  # Only descriptive estimates are returned; no asymptotic inference is made.
  pearson <- suppressWarnings(indep.test(x, method = "pearson",
                                         correct = FALSE))
  n <- sum(pearson$observed)
  chisq <- unname(pearson$statistic)
  d <- min(dim(pearson$observed))
  estimates <- c(phi = sqrt(chisq / n),
                 contingency = sqrt(chisq / (chisq + n)),
                 `cramer-v` = sqrt(chisq / (n * (d - 1))))

  structure(
    list(estimate = estimates[measures], measure = measures,
         statistic = pearson$statistic, n = n,
         observed = pearson$observed, expected = pearson$expected,
         data.name = paste(deparse(substitute(x)), collapse = "\n")),
    class = "ibist_nominal_assoc"
  )
}

#' @rdname nominal.assoc
#' @param digits Number of significant digits used for printing.
#' @param ... Further arguments passed to \code{print.default()}.
#' @export
print.ibist_nominal_assoc <- function(x, digits = 7L, ...) {
  cat("\nNominal association measures\n\ndata:", x$data.name, "\n\n")
  print(cbind(Estimate = x$estimate), digits = digits, ...)
  invisible(x)
}
