% Generated from R/data-asah.R: do not edit by hand
\docType{data}
\name{asah}
\alias{asah}
\title{Observer Ratings for Aneurysmal Subarachnoid Hemorrhage Scales}
\format{
A data frame with 309 rows (103 paired ratings for each of 3 scales)
and 4 variables:
\describe{
  \item{scale}{Grading scale: WFNS, HH, or PAASH.}
  \item{assessment}{Sequential row number within scale. These numbers are
  generated for the reconstructed data and are not study participant IDs.}
  \item{observer1}{Ordered factor giving the first observer's grade (1--5).}
  \item{observer2}{Ordered factor giving the second observer's grade (1--5).}
}
}
\description{
Reconstructed paired observer ratings for three clinical grading scales
used after aneurysmal subarachnoid hemorrhage: WFNS, Hunt and Hess (HH),
and PAASH.
}
\details{
The paired ratings are expanded from the published cross-tabulations. They
reproduce the observed cell counts for each scale and can be used to
illustrate agreement summaries such as kappa. The original participant-level
records are not available here; generated rows cannot be linked across scales.
}
\source{
Degen, L. A. R., Dorhout Mees, S. M., Algra, A., and Rinkel, G. J. E. (2011).
Interobserver variability of grading scales for aneurysmal subarachnoid
hemorrhage. \emph{Stroke}, 42(6), 1546--1549.
\doi{10.1161/STROKEAHA.110.601211}
}
\examples{
data(asah)
wfns <- subset(asah, scale == "WFNS")
table(wfns$observer1, wfns$observer2)
# irr::kappa2(wfns[c("observer1", "observer2")], weight = "squared")
}
\keyword{datasets}
