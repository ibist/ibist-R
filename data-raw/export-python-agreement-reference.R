# Run from ibist-R; Python also checks independent variance calculations.
source("R/cohen_kappa_test.R")
source("R/agreement_test.R")
examples <- list(
  matrix(c(10, 4, 1, 0, 5, 10, 12, 2, 2, 4, 12, 5, 0, 2, 6, 13),
         4, byrow = TRUE),
  matrix(c(4, 0, 0, 0, 4, 0, 0, 4, 0), 3, byrow = TRUE),
  matrix(c(9, 1, 2, 8), 2, byrow = TRUE),
  matrix(c(1, 8, 9, 2), 2, byrow = TRUE),
  matrix(c(581, 134, 12, 6, 0, 2, 91, 73, 18, 31, 3, 1,
           4, 18, 15, 19, 3, 3, 1, 8, 7, 24, 11, 6,
           0, 0, 2, 4, 3, 5, 0, 0, 0, 9, 11, 13), 6, byrow = TRUE),
  matrix(4, 3, 3), diag(c(4, 4)), diag(c(1, 2, 3)),
  matrix(c(4, 0, 1, 0, 0, 0, 2, 0, 5), 3, byrow = TRUE),
  matrix(c(0, 4, 0, 2, 0, 0, 0, 0, 0), 3, byrow = TRUE),
  matrix(c(1, 0, 0, 0), 2), matrix(0, 2, 2)
)
set.seed(20261108)
for (i in 1:25) {
  size <- sample(2:6, 1)
  examples[[length(examples) + 1L]] <- matrix(sample(0:30, size^2, TRUE), size)
}
# Cover dimensions 2--6 and key boundaries without a full Cartesian grid.
size.five <- which(vapply(examples, nrow, integer(1)) == 5L)[1L]
examples <- examples[c(1, 2, 4, 5, size.five, 7, 8, 9, 10, 11, 12)]
cases <- list()
for (j in seq_along(examples)) {
  tab <- examples[[j]]
  size <- nrow(tab)
  custom <- diag(size)
  custom[upper.tri(custom)] <- rep(c(0.2, 0.6, 0.1),
                                  length.out = choose(size, 2))
  custom[lower.tri(custom)] <- t(custom)[lower.tri(custom)]
  for (weights in list("unweighted", "linear", "quadratic", custom)) {
    for (variance in c("fce", "delta")) {
      for (alternative in c("two.sided", "greater", "less")) {
        name <- if (is.character(weights)) weights else "custom"
        keep <- if (j <= 5L) {
          alternative == "two.sided" ||
            (j == 1L && name == "linear") ||
            (j == 4L && name == "quadratic") ||
            (j == 2L && name == "custom")
        } else {
          alternative == "two.sided" && switch(as.character(j),
            "6" = name == "unweighted",
            "7" = name == "quadratic" && variance == "delta",
            "8" = name %in% c("linear", "custom") && variance == "fce",
            "9" = name == "unweighted" && variance == "fce",
            "10" = name == "unweighted" && variance == "fce",
            "11" = name == "unweighted" && variance == "fce")
        }
        if (!keep) next
        input <- list(
          x = lapply(seq_len(size), function(i) unname(tab[i, ])),
          weights = if (is.matrix(weights)) {
            lapply(seq_len(size), function(i) unname(weights[i, ]))
          } else weights,
          conf_level = 0.9, alternative = alternative,
          test_variance = variance)
        result <- tryCatch(
          agreement.test(tab, weights = weights, conf.level = 0.9,
                         alternative = alternative, test.variance = variance),
          error = function(error) error)
        if (inherits(result, "error")) {
          cases[[length(cases) + 1L]] <- c(input, list(error = TRUE))
        } else {
          cases[[length(cases) + 1L]] <- c(input, list(
            estimate = unname(result$estimate),
            estimate_name = names(result$estimate),
            statistic = unname(result$statistic), p_value = result$p.value,
            conf_int = as.list(as.numeric(result$conf.int)),
            std_error = result$std.error,
            null_std_error = result$null.std.error,
            weight_method = result$weight.method, method = result$method,
            weight_matrix = lapply(seq_len(size), function(i) {
              as.list(unname(result$weights[i, ]))
            }), n = result$n, null_value = unname(result$null.value)))
        }
      }
    }
  }
}
files <- c("R/cohen_kappa_test.R", "R/agreement_test.R")
jsonlite::write_json(
  list(source = "ibist-R working sources; independent Python checks also used",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(files)), cases = cases),
  "../ibist-py/tests/reference/agreement_r.json", pretty = TRUE,
  auto_unbox = TRUE, digits = NA, na = "string")
message("Exported ", length(cases), " agreement cases.")
