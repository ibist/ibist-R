#include <Rcpp.h>
#include "wang_core.h"

using Rcpp::NumericVector;


// [[Rcpp::export]]
NumericVector wang_paired_ci_cpp(int n10,
                                 int t,
                                 int n01,
                                 double conf_level = 0.95,
                                 std::string ci_type = "Two.sided",
                                 double precision = 0.00001,
                                 int grid_one = 30,
                                 int grid_two = 20) {
  if (n10 < 0 || t < 0 || n01 < 0) Rcpp::stop("Counts must be non-negative.");
  if (!(conf_level > 0.0 && conf_level < 1.0)) {
    Rcpp::stop("conf_level must be between 0 and 1.");
  }
  if (!(precision > 0.0)) Rcpp::stop("precision must be positive.");
  if (grid_one < 2 || grid_two < 2) {
    Rcpp::stop("Grid sizes must be at least 2.");
  }

  const int n = n10 + t + n01;
  const double estimate = n == 0 ? 0.0 :
    (static_cast<double>(n10) - static_cast<double>(n01)) /
      static_cast<double>(n);

  double lower = -1.0;
  double upper = 1.0;

  if (ci_type == "Lower") {
    lower = ibist_wang::wang_lower_one_sided(
      n10, t, n01, conf_level, precision, grid_one, grid_two);
  } else if (ci_type == "Upper") {
    upper = -ibist_wang::wang_lower_one_sided(
      n01, t, n10, conf_level, precision, grid_one, grid_two);
  } else if (ci_type == "Two.sided") {
    const double one_sided_level = 1.0 - (1.0 - conf_level) / 2.0;
    lower = ibist_wang::wang_lower_one_sided(
      n10, t, n01, one_sided_level, precision,
      grid_one, grid_two);
    upper = -ibist_wang::wang_lower_one_sided(
      n01, t, n10, one_sided_level, precision,
      grid_one, grid_two);
  } else {
    Rcpp::stop("ci_type must be 'Lower', 'Upper', or 'Two.sided'.");
  }

  NumericVector out = NumericVector::create(
    Rcpp::Named("estimate") = estimate,
    Rcpp::Named("lower") = lower,
    Rcpp::Named("upper") = upper
  );
  return out;
}

// [[Rcpp::export]]
NumericVector wang_paired_lower_all_cpp(int n,
                                        double conf_level = 0.95,
                                        double precision = 0.00001,
                                        int grid_one = 30,
                                        int grid_two = 20) {
  if (n < 0) Rcpp::stop("n must be non-negative.");
  if (!(conf_level > 0.0 && conf_level < 1.0)) {
    Rcpp::stop("conf_level must be between 0 and 1.");
  }
  if (!(precision > 0.0)) Rcpp::stop("precision must be positive.");
  if (grid_one < 2 || grid_two < 2) {
    Rcpp::stop("Grid sizes must be at least 2.");
  }

  NumericVector result = Rcpp::wrap(ibist_wang::wang_lower_all_one_sided(
    n, conf_level, precision, grid_one, grid_two));
  for (R_xlen_t i = 0; i < result.size(); ++i) {
    if (std::isnan(result[i])) result[i] = NA_REAL;
  }
  return result;
}
