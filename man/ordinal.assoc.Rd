% Generated from R/ordinal_assoc.R: do not edit by hand
\name{ordinal.assoc}
\alias{ordinal.assoc}
\alias{print.ibist_ordinal_assoc}
\title{Ordinal Association Measures with Asymptotic Confidence Intervals and Test}
\usage{
ordinal.assoc(
  x,
  measure = c("gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c"),
  conf.level = 0.95,
  alternative = c("two.sided", "greater", "less")
)
}
\arguments{
\item{x}{A two-dimensional matrix or table of non-negative integer counts.
Rows and columns must be ordered from the lowest to highest category.}

\item{measure}{One or more measures: \code{"gamma"}, \code{"tau-b"},
\code{"tau-c"}, \code{"somers-c|r"}, and \code{"somers-r|c"}. By
default, all five are computed.}

\item{conf.level}{Confidence level for the two-sided asymptotic confidence
intervals.}

\item{alternative}{Alternative hypothesis for the shared test:
\code{"two.sided"}, \code{"greater"}, or \code{"less"}. Positive
association means that higher row categories tend to occur with higher
column categories.}
}
\value{
An object of class \code{"ibist_ordinal_assoc"} containing:
\describe{
  \item{estimate}{Named vector of selected association estimates.}
  \item{std.error}{Named vector of asymptotic standard errors.}
  \item{conf.int}{Matrix of lower and upper confidence limits.}
  \item{statistic}{The shared asymptotic Z statistic.}
  \item{p.value}{P-value for the selected alternative.}
  \item{S, var.S}{The observed concordance contrast and its null variance.}
  \item{concordant, discordant}{Concordant and discordant pair counts.}
}
}
\description{
Computes concordance-based ordinal association measures and their asymptotic
confidence intervals from an RxC table. It also tests the common null
hypothesis that the association measure is zero.
}
\details{
The procedure computes the concordant and discordant pair counts, n_C and
n_D, and sets S = n_C - n_D. Its null variance is
\deqn{\mathrm{Var}_0(S) = \sum_{ij} n_{ij}(A_{ij}-D_{ij})^2 - 4S^2/n,}
where A_ij and D_ij are the concordant and discordant contributions for cell
(i,j). The test statistic is Z = S/sqrt(Var_0(S)).

The selected measures have different sampling variances, so their standard
errors and confidence intervals differ. Their tests of a zero association
nevertheless reduce to the same Z statistic and p-value. Confidence intervals
are two-sided regardless of the test alternative. Rows and columns with zero
marginal totals are dropped before computation.

The \code{"somers-c|r"} measure is Somers' Delta(C|R), treating the row
variable as explanatory; \code{"somers-r|c"} is Delta(R|C), treating the
column variable as explanatory.
}
\examples{
tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), nrow = 3, byrow = TRUE)
ordinal.assoc(tab)
ordinal.assoc(tab, measure = c("tau-b", "somers-c|r"))
}
\references{
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
