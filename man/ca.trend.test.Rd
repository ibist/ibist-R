% Generated from R/ca_trend_test.R: do not edit by hand
\name{ca.trend.test}
\alias{ca.trend.test}
\title{Cochran-Armitage Trend Test}
\usage{
ca.trend.test(
  x,
  scores = NULL,
  alternative = c("two.sided", "greater", "less"),
  exact = FALSE
)
}
\arguments{
\item{x}{An \eqn{R \times 2} matrix or table of non-negative integer
counts. Rows are ordered groups and columns are outcome categories. The
first column is treated as the event of interest.}

\item{scores}{Optional numeric scores for the ordered groups. Defaults to
\code{1:R}; scores must be finite and strictly increasing.}

\item{alternative}{Alternative hypothesis: \code{"two.sided"},
\code{"greater"} (increasing trend), or \code{"less"} (decreasing
trend).}

\item{exact}{Logical; if \code{TRUE}, compute the exact conditional
p-value given both margins. Otherwise, use the asymptotic normal
approximation.}
}
\value{
An object of class \code{"htest"}. The signed Cochran-Armitage statistic
is returned along with its square, the score statistic \eqn{T}, the
p-value, and the alternative hypothesis.
}
\description{
Tests for a linear trend in a binary outcome across ordered groups using
the Cochran-Armitage statistic. The test can use a normal approximation
or an exact conditional distribution.
}
\details{
A positive signed statistic indicates that the event proportion tends to
increase with the group scores. The exact test conditions on the row totals
and the total number of events and uses the conditional distribution of
\eqn{T = \sum_i x_i n_{i1}}. For a two-sided exact test, the p-value is
twice the smaller one-sided tail probability, truncated at one. The exact
calculation is performed by dynamic programming and can become slow when
the table has many rows or large margins.
}
\examples{
tab <- matrix(c(106, 109, 128, 130, 99, 95, 76, 74), nrow = 4)
ca.trend.test(tab)
ca.trend.test(tab, scores = c(1, 2, 4, 8), alternative = "greater")
small.tab <- matrix(c(2, 0, 1, 1, 0, 2), nrow = 3, byrow = TRUE)
ca.trend.test(small.tab, exact = TRUE)
}
\references{
Cochran, W. G. (1954). Some methods for strengthening the common
chi-squared tests. \emph{Biometrics}, 10(4), 417--451.
\doi{10.2307/3001616}

Armitage, P. (1955). Tests for linear trends in proportions and
frequencies. \emph{Biometrics}, 11(3), 375--386.
\doi{10.2307/3001775}

Portier, C., and Hoel, D. G. (1984). Type 1 error of trend tests in
proportions and the design of cancer screens.
\emph{Communications in Statistics - Theory and Methods}, 13(1), 1--14.
}
