# R common-OR tests supplement independent Python tail/enumeration checks.
# Run from ibist-R with Rscript --vanilla.
for (file in c("or_2x2_ci", "or_common")) source(paste0("R/", file, ".R"))
examples <- list(
  list(c(3, 1, 1, 8)),
  list(c(7, 27, 1, 33)),
  list(c(2, 3, 1, 4), c(1, 2, 3, 2)),
  list(c(8, 2, 3, 1), c(9, 1, 4, 2)),
  list(c(10, 1, 12, 1), c(9, 0, 11, 1), c(8, 0, 7, 3)),
  list(c(0, 3, 2, 0), c(0, 4, 1, 0)),
  list(c(3, 0, 0, 2)), list(c(0, 0, 0, 0)),
  list(c(5, 0, 0, 0)), list(c(0, 0, 5, 2), c(7, 3, 2, 8))
)
set.seed(20261107)
for (i in 1:40) {
  examples[[length(examples) + 1L]] <- replicate(
    sample(1:5, 1), sample(0:20, 4, TRUE), simplify = FALSE)
}
cases <- list()
add <- function(example, theta) {
  tabs <- lapply(example, matrix, nrow = 2, byrow = TRUE)
  rows <- function(tab) lapply(1:2, function(i) unname(tab[i, ]))
  for (alternative in c("two.sided", "greater", "less")) {
    result <- or.common.test(tabs, theta = theta, alternative = alternative)
    cases[[length(cases) + 1L]] <<- list(
      x = lapply(tabs, rows), theta = theta, alternative = alternative,
      statistic = unname(result$statistic), p_value = result$p.value,
      estimate = unname(result$estimate),
      null_value = unname(result$null.value),
      parameter = unname(result$parameter), n_strata = length(tabs),
      method = result$method
    )
  }
}
for (example in examples) for (theta in c(0.25, 1, 4)) add(example, theta)
add(rep(list(c(1, 0, 0, 1)), 2), 1e-150)
add(rep(list(c(0, 1, 1, 0)), 2), 1e150)
files <- c("R/or_2x2_ci.R", "R/or_common.R")
jsonlite::write_json(
  list(source = "ibist-R working sources; independent Python checks also used",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(files)), cases = cases),
  "../ibist-py/tests/reference/common_test_r.json", pretty = TRUE,
  auto_unbox = TRUE, digits = NA, na = "string"
)
message("Exported ", length(cases), " common-OR test cases.")
