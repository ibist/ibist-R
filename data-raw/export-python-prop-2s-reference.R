# Explicitly regenerate cross-language regression cases, then review changes.
# Python also checks independent formulas and constrained likelihood fits.
# Run from the R package root:
# Rscript data-raw/export-python-prop-2s-reference.R

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite to export development reference cases.")
}
source("R/ci.R")
source("R/prop_1s_ci.R")
source("R/prop_2s_ci.R")
methods <- c("newcombe", "wald", "miettinen-nurminen", "mee",
             "agresti-caffo", "brown-li")
cases <- list()
add <- function(x, n, level = 0.95, alternative = "two.sided",
                correct = FALSE, selected = methods) {
  result <- prop.2s.ci(x, n, conf.level = level, method = selected,
                        alternative = alternative, correct = correct)
  cases[[length(cases) + 1L]] <<- list(
    x = x, n = n,
    options = list(conf_level = level, method = selected,
                   alternative = alternative, correct = correct),
    estimate = unname(result$estimate), methods = as.list(result$method),
    sample_proportions = unname(result$sample.proportions),
    conf_int = lapply(seq_len(nrow(result$conf.int)), function(i) {
      unname(result$conf.int[i, ])
    })
  )
}
examples <- list(
  list(x = c(1366, 383), n = c(15571, 5366)),
  list(x = c(38, 267), n = c(100, 984)),
  list(x = c(56, 48), n = c(70, 80)),
  list(x = c(9, 3), n = c(10, 10)),
  list(x = c(10, 0), n = c(10, 20)),
  list(x = c(12, 4), n = c(40, 30)),
  list(x = c(0, 0), n = c(3, 5)),
  list(x = c(0, 5), n = c(3, 5)),
  list(x = c(3, 0), n = c(3, 5)),
  list(x = c(3, 5), n = c(3, 5)),
  list(x = c(2, 4), n = c(3, 5)),
  list(x = c(0, 1), n = c(1, 1))
)
for (case in examples) {
  for (level in c(0.8, 0.95, 0.99)) {
    for (alternative in c("two.sided", "less", "greater")) {
      add(case$x, case$n, level, alternative)
      add(case$x, case$n, level, alternative, correct = TRUE,
          selected = c("newcombe", "wald"))
    }
  }
}
set.seed(20261106)
for (i in 1:30) {
  n <- sample(5:100, 2, TRUE)
  x <- c(sample(0:n[1], 1), sample(0:n[2], 1))
  add(x, n)
}
for (level in c(0.2, 0.5)) {
  for (alternative in c("less", "greater")) {
    add(c(4, 1), c(6, 5), level, alternative)
    add(c(4, 1), c(6, 5), level, alternative, correct = TRUE)
  }
}
add(c(9, 3), c(10, 10), selected = c("mee", "newcombe", "mee"))
output <- "../ibist-py/tests/reference/prop_2s_r.json"
jsonlite::write_json(
  list(source = "R source in ibist-R; independently checked in Python tests",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(c("R/ci.R", "R/prop_1s_ci.R",
                                          "R/prop_2s_ci.R"))), cases = cases),
  output, pretty = TRUE, auto_unbox = TRUE, digits = NA
)
message("Exported ", length(cases), " R regression cases to ", output)
