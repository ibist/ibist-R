#' CMH Linear-by-Linear Association Test
#'
#' Tests for a zero score correlation in an RxC contingency table using the
#' Cochran-Mantel-Haenszel linear-by-linear statistic.
#'
#' @param x A two-dimensional matrix or table of non-negative integer counts.
#'   Rows and columns must be ordered from the lowest to highest category.
#' @param measure The score correlation to test: \code{"pearson"} or
#'   \code{"spearman"}. Pearson correlation uses the supplied category scores
#'   or, by default, numeric category labels and otherwise category order
#'   numbers. Spearman correlation uses marginal midrank scores.
#' @param alternative Alternative hypothesis: \code{"two.sided"},
#'   \code{"greater"}, or \code{"less"}. A positive association means that
#'   larger row scores tend to occur with larger column scores.
#' @param row.scores,col.scores Optional numeric scores for row and column
#'   categories when \code{measure = "pearson"}. These must be specified
#'   together. They are not used for Spearman correlation.
#'
#' @details
#' Let \eqn{r} be the Pearson correlation between the assigned row and column
#' scores (or their marginal midranks for Spearman correlation). The signed
#' test statistic is \eqn{Z = r\sqrt{n-1}}, and the CMH statistic is
#' \eqn{X^2 = Z^2}, which has an asymptotic chi-squared distribution with one
#' degree of freedom under the null hypothesis. One-sided tests use the signed
#' normal statistic. This function tests zero association only; it does not
#' return a confidence interval. Use \code{ordinal.assoc()} for the
#' multinomial-based Wald confidence interval for Pearson or Spearman
#' correlation. Those intervals are not obtained by inverting this test. The
#' Pearson and Spearman tests in \code{ordinal.assoc()} use their own
#' multinomial null variances and are distinct from this CMH test; their
#' p-values can differ.
#'
#' @return An object of class \code{"htest"}. The \code{statistic} element
#'   contains the signed \eqn{Z} statistic, and the \code{chisq} element
#'   contains its square, the CMH statistic. The \code{estimate} element
#'   contains the score correlation.
#'
#' @references
#' Mantel, N. (1963). Chi-square tests with one degree of freedom: extensions
#' of the Mantel-Haenszel procedure. \emph{Journal of the American Statistical
#' Association}, 58(303), 690--700.
#' \doi{10.1080/01621459.1963.10500879}
#'
#' @examples
#' tab <- matrix(c(109, 68, 14, 72, 79, 37, 9, 42, 139), nrow = 3,
#'               byrow = TRUE)
#' cmh.ll.test(tab, measure = "pearson")
#' cmh.ll.test(tab, measure = "spearman", alternative = "greater")
#' ordinal.assoc(tab, measure = "spearman")
#'
#' @export
cmh.ll.test <- function(
  x,
  measure = c("pearson", "spearman"),
  alternative = c("two.sided", "greater", "less"),
  row.scores = NULL,
  col.scores = NULL
) {
  measure <- match.arg(measure)
  alternative <- match.arg(alternative)

  tab <- as.matrix(x)
  if (length(dim(tab)) != 2L || any(dim(tab) < 2L) || !is.numeric(tab) ||
      anyNA(tab) || any(!is.finite(tab)) || any(tab < 0) ||
      any(tab != floor(tab))) {
    stop("'x' must be a two-dimensional table of non-negative integer counts.")
  }

  if (xor(is.null(row.scores), is.null(col.scores))) {
    stop("'row.scores' and 'col.scores' must be supplied together.")
  }
  if (measure == "spearman" &&
      (!is.null(row.scores) || !is.null(col.scores))) {
    stop("Category scores cannot be supplied for Spearman correlation.")
  }

  keep.rows <- rowSums(tab) > 0
  keep.cols <- colSums(tab) > 0
  if (!is.null(row.scores)) {
    row.scores <- validate_ordinal_scores(
      row.scores,
      nrow(tab),
      "row.scores"
    )[keep.rows]
    col.scores <- validate_ordinal_scores(
      col.scores,
      ncol(tab),
      "col.scores"
    )[keep.cols]
  }
  tab <- tab[keep.rows, keep.cols, drop = FALSE]
  if (any(dim(tab) < 2L)) {
    stop("'x' must have at least two non-empty rows and columns.")
  }

  if (measure == "pearson") {
    if (is.null(row.scores)) {
      row.scores <- default_ordinal_scores(rownames(tab), nrow(tab))
      col.scores <- default_ordinal_scores(colnames(tab), ncol(tab))
    }
    estimate <- ordinal_assoc_score_correlation(
      tab,
      row.scores,
      col.scores
    )$estimate
  } else {
    estimate <- ordinal_assoc_spearman(tab)$estimate
  }

  n <- sum(tab)
  z <- sqrt(n - 1) * estimate
  statistic <- z^2
  p.value <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(-abs(z)),
    greater = stats::pnorm(z, lower.tail = FALSE),
    less = stats::pnorm(z)
  )
  structure(
    list(
      statistic = c(Z = z),
      chisq = c(`CMH X-squared` = statistic),
      p.value = p.value,
      estimate = c(correlation = estimate),
      null.value = c(correlation = 0),
      alternative = alternative,
      method = paste("CMH linear-by-linear association test for",
                     measure, "correlation"),
      data.name = deparse(substitute(x))
    ),
    class = "htest"
  )
}
