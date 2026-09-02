#' Blood Pressure Cohort Data (midbp)
#'
#' Aggregated data of cases and person-years by diastolic blood pressure
#' category and gender.
#'
#' @format A data frame with 16 rows and 4 variables:
#' \describe{
#'   \item{dbp}{Diastolic blood pressure category (1--8).}
#'   \item{gender}{Gender (Men, Women).}
#'   \item{cases}{Number of observed cases.}
#'   \item{pyears}{Person-years of follow-up.}
#' }
#'
#' @details
#' The dataset is aggregated by DBP category and gender, suitable for
#' rate modeling (e.g., Poisson regression with offset).
#'
#' The DBP categories are coded as:
#' 1 = <70 mmHg, 2 = 70--74 mmHg, 3 = 75--79 mmHg,
#' 4 = 80--84 mmHg, 5 = 85--89 mmHg, 6 = 90--94 mmHg,
#' 7 = 95--99 mmHg, and 8 = >=100 mmHg.
#'
#' @examples
#' data(midbp)
#' head(midbp)
#'
"midbp"
