% Generated from R/kappa_test.R: do not edit by hand
\name{kappa.test}
\alias{kappa.test}
\title{Asymptotic Inference for Cohen's Kappa}
\usage{
kappa.test(
  x,
  weights = "unweighted",
  conf.level = 0.95,
  alternative = c("two.sided", "greater", "less")
)
}
\arguments{
\item{x}{A square matrix or table of non-negative integer counts. Rows and
columns represent the same ordered rating categories.}

\item{weights}{The agreement weights: \code{"unweighted"}, \code{"linear"},
\code{"quadratic"}, or a numeric matrix of user-specified weights. The
named choices use category positions in \code{x}; the quadratic weights
are Fleiss--Cohen weights. A supplied matrix must have the same dimensions
as \code{x}, be symmetric, have ones on the diagonal, and contain values
between zero and one.}

\item{conf.level}{Confidence level for the two-sided Wald confidence
interval.}

\item{alternative}{Alternative hypothesis: \code{"two.sided"},
\code{"greater"}, or \code{"less"}.}
}
\value{
An object of class \code{"htest"} containing the estimate, confidence
interval, Z statistic, and p-value. It also includes the standard error
used for the confidence interval, the standard error under the null
hypothesis, the weight matrix, and the sample size.
}
\description{
Computes Cohen's kappa or weighted kappa from a square table of paired
ratings, with a Wald confidence interval and an asymptotic test of
\eqn{H_0: \kappa = 0}.
}
\details{
The confidence interval uses the large-sample variance of the kappa
estimate. The test uses the variance under the null hypothesis
\eqn{\kappa=0}; these variances are not interchangeable. The formulas
follow Fleiss, Cohen, and Everitt (1969) and the corresponding SAS
\code{PROC FREQ} formulas. The interval is an untruncated Wald interval.

This function provides asymptotic inference only. SAS \code{PROC FREQ}
also provides exact tests for simple and weighted kappa; exact inference is
not implemented here.
}
\examples{
peff <- matrix(
  c(581, 134, 12, 6, 0, 2,
    91, 73, 18, 31, 3, 1,
    4, 18, 15, 19, 3, 3,
    1, 8, 7, 24, 11, 6,
    0, 0, 2, 4, 3, 5,
    0, 0, 0, 9, 11, 13),
  nrow = 6,
  byrow = TRUE
)
kappa.test(peff)
kappa.test(peff, weights = "quadratic")
}
\references{
Fleiss, J. L., Cohen, J., and Everitt, B. S. (1969). Large sample standard
errors of kappa and weighted kappa. \emph{Psychological Bulletin}, 72(5),
323--327. \doi{10.1037/h0028106}

SAS Institute Inc. (2015). Tests and measures of agreement. In
\emph{Base SAS 9.4 Procedures Guide: Statistical Procedures, Third Edition}.
\url{https://support.sas.com/documentation/cdl/en/procstat/67528/HTML/default/procstat_freq_details76.htm}
}
