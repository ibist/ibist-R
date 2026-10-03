% Generated from R/cmh_ll_test.R: do not edit by hand
\name{cmh.ll.test}
\alias{cmh.ll.test}
\title{CMH Linear-by-Linear Association Test}
\usage{
cmh.ll.test(
  x,
  measure = c("pearson", "spearman"),
  alternative = c("two.sided", "greater", "less"),
  row.scores = NULL,
  col.scores = NULL
)
}
\arguments{
\item{x}{A two-dimensional matrix or table of non-negative integer counts.
Rows and columns must be ordered from the lowest to highest category.}

\item{measure}{The score correlation to test: \code{"pearson"} or
\code{"spearman"}. Pearson correlation uses the supplied category scores
or, by default, numeric category labels and otherwise category order
numbers. Spearman correlation uses marginal midrank scores.}

\item{alternative}{Alternative hypothesis: \code{"two.sided"},
\code{"greater"}, or \code{"less"}. A positive association means that
larger row scores tend to occur with larger column scores.}

\item{row.scores, col.scores}{Optional numeric scores for row and column
categories when \code{measure = "pearson"}. These must be specified
together. They are not used for Spearman correlation.}
}
\value{
An object of class \code{"htest"}. The \code{statistic} element contains
the signed \eqn{Z} statistic, and the \code{chisq} element contains its
square, the CMH statistic. The \code{estimate} element contains the score
correlation.
}
\description{
Tests for a zero score correlation in an RxC contingency table using the
Cochran-Mantel-Haenszel linear-by-linear statistic.
}
\details{
Let \eqn{r} be the Pearson correlation between the assigned row and column
scores (or their marginal midranks for Spearman correlation). The signed
test statistic is \eqn{Z = r\sqrt{n-1}}, and the CMH statistic is
\eqn{X^2 = Z^2}, which has an asymptotic chi-squared distribution with one
degree of freedom under the null hypothesis. One-sided tests use the signed
normal statistic. This function tests zero association only; it does not
return a confidence interval. Use \code{ordinal.assoc()} for the
multinomial-based Wald confidence interval for Pearson or Spearman
correlation. Those intervals are not obtained by inverting this test. The
Pearson and Spearman tests in \code{ordinal.assoc()} use their own
multinomial null variances and are distinct from this CMH test; their
p-values can differ.
}
\examples{
tab <- matrix(c(109, 68, 14, 72, 79, 37, 9, 42, 139), nrow = 3,
              byrow = TRUE)
cmh.ll.test(tab, measure = "pearson")
cmh.ll.test(tab, measure = "spearman", alternative = "greater")
ordinal.assoc(tab, measure = "spearman")
}
\references{
Mantel, N. (1963). Chi-square tests with one degree of freedom: extensions
of the Mantel-Haenszel procedure. \emph{Journal of the American Statistical
Association}, 58(303), 690--700. \doi{10.1080/01621459.1963.10500879}
}
