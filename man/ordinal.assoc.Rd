% Generated from R/ordinal_assoc.R: do not edit by hand
\name{ordinal.assoc}
\alias{ordinal.assoc}
\alias{print.ibist_ordinal_assoc}
\title{Ordinal Association Measures with Asymptotic Confidence Intervals and Tests}
\usage{
ordinal.assoc(
  x,
  measure = c(
    "gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c",
    "pearson", "spearman"
  ),
  conf.level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  row.scores = NULL,
  col.scores = NULL
)
}
\arguments{
\item{x}{A two-dimensional matrix or table of non-negative integer counts.
Rows and columns must be ordered from the lowest to highest category.}

\item{measure}{One or more measures: \code{"gamma"}, \code{"tau-b"},
\code{"tau-c"}, \code{"somers-c|r"}, \code{"somers-r|c"},
\code{"pearson"}, and \code{"spearman"}. By default, all seven are
computed.}

\item{conf.level}{Confidence level for the two-sided asymptotic confidence
intervals.}

\item{alternative}{Alternative hypothesis for the tests:
\code{"two.sided"}, \code{"greater"}, or \code{"less"}. Positive
association means that higher row categories tend to occur with higher
column categories.}

\item{row.scores, col.scores}{Optional numeric scores for row and column
categories when computing Pearson correlation. By default, numeric
category labels are used; otherwise, category order numbers are used.
Spearman correlation always uses marginal midrank scores.}
}
\value{
An object of class \code{"ibist_ordinal_assoc"} containing:
\describe{
  \item{estimate}{Named vector of selected association estimates.}
  \item{std.error}{Named vector of asymptotic standard errors.}
  \item{conf.int}{Matrix of lower and upper confidence limits.}
  \item{statistic}{Named asymptotic test statistics, one per selected
  measure.}
  \item{p.value}{Named p-values, one per selected measure.}
  \item{S, var.S}{The observed concordance contrast and its null variance.}
  \item{concordant, discordant}{Concordant and discordant pair counts.}
}
}
\description{
Computes ordinal association measures and their asymptotic confidence
intervals from an RxC table. It also tests the null hypothesis that each
association measure is zero.
}
\details{
The procedure computes the concordant and discordant pair counts, n_C and
n_D, and sets S = n_C - n_D. Its null variance is
\deqn{\mathrm{Var}_0(S) = \sum_{ij} n_{ij}(A_{ij}-D_{ij})^2 - 4S^2/n,}
where A_ij and D_ij are the concordant and discordant contributions for cell
(i,j). The test statistic is Z = S/sqrt(Var_0(S)).

The concordance-based measures share a test statistic, but Pearson and
Spearman correlations use their own multinomial null variances. Standard
errors and confidence intervals are computed separately for each measure.
Confidence intervals are two-sided regardless of the test alternative.
Rows and columns with zero marginal totals are dropped before computation.

The \code{"somers-c|r"} measure is Somers' Delta(C|R), treating the row
variable as explanatory; \code{"somers-r|c"} is Delta(R|C), treating the
column variable as explanatory.

Pearson correlation uses the supplied category scores or the default
category scores described above (numeric labels such as 1, 2, and 10 are
used as-is). Spearman correlation uses marginal
midranks, so observations tied within a category receive the same rank.
The correlation standard errors and null variances use multinomial
contingency-table formulas. The confidence intervals are Wald intervals.
}
\examples{
tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), nrow = 3, byrow = TRUE)
ordinal.assoc(tab)
ordinal.assoc(
  tab,
  measure = c("tau-b", "somers-c|r", "pearson", "spearman")
)
}
\references{
Brown, M. B., and Benedetti, J. K. (1977). Sampling behavior of tests for
correlation in two-way contingency tables. \emph{Journal of the American
Statistical Association}, 72(358), 309--315.

Goodman, L. A., and Kruskal, W. H. (1954). Measures of association for
cross classifications. \emph{Journal of the American Statistical
Association}, 49(268), 732--764.

Kendall, M. G. (1970). \emph{Rank Correlation Methods} (4th ed.). Griffin.

Somers, R. H. (1962). A new asymmetric measure of association for ordinal
variables. \emph{American Sociological Review}, 27(6), 799--811.

Stuart, A. (1953). The estimation and comparison of strengths of association
in contingency tables. \emph{Biometrika}, 40(1--2), 105--110.
}
\keyword{htest}
