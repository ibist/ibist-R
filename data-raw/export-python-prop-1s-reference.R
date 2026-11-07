# Explicitly regenerate cross-language regression cases, then review changes.
# Independent formula and coverage tests are maintained in Python.
# Run from the R package root:
# Rscript data-raw/export-python-prop-1s-reference.R

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite to export development reference cases.")
}
source("R/ci.R")
source("R/prop_1s_ci.R")
methods <- c("wilson", "wald", "agresti-coull", "clopper-pearson", "blaker")
cases <- list()
add <- function(x, n, level = 0.95, alternative = "two.sided",
                selected = methods, correct = FALSE, midp = FALSE) {
  result <- prop.1s.ci(x, n, conf.level = level, method = selected,
                        alternative = alternative, correct = correct,
                        midp = midp)
  inputs <- result$inputs
  inputs$input <- inputs$input - 1L
  cases[[length(cases) + 1L]] <<- list(
    x = x, n = n,
    options = list(conf_level = level, method = selected, correct = correct,
                   midp = midp, alternative = alternative),
    estimate = unname(result$estimate),
    conf_int = lapply(seq_len(nrow(result$conf.int)), function(i) {
      unname(result$conf.int[i, ])
    }),
    methods = as.list(sub("\\[[0-9]+\\]$", "", result$method)),
    inputs = inputs
  )
}
for (n in c(1, 2, 10, 40)) {
  for (x in unique(c(0, 1, floor(n / 2), n))) {
    for (level in c(0.8, 0.95, 0.99)) {
      for (alternative in c("two.sided", "greater", "less")) {
        add(x, n, level, alternative)
        add(x, n, level, alternative, selected = c("wilson", "wald"),
            correct = TRUE)
        add(x, n, level, alternative, selected = "clopper-pearson",
            midp = TRUE)
      }
    }
  }
}
add(341, 762, level = 0.9)
add(c(0, 5, 40), 40, correct = TRUE, midp = TRUE)
add(c(0, 3, 12), c(1, 4, 12))
add(c(0, 3, 12), c(1, 4, 12), alternative = "greater")
add(3, 10, selected = c("wilson", "wilson", "blaker"))
for (alternative in c("greater", "less")) {
  add(c(0, 1, 5, 10), 10, level = 0.2, alternative = alternative,
      correct = TRUE, midp = TRUE)
}
output <- "../ibist-py/tests/reference/prop_1s_r.json"
jsonlite::write_json(
  list(source = "R source in ibist-R; independently checked in Python tests",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(c("R/ci.R", "R/prop_1s_ci.R"))),
       cases = cases), output, pretty = TRUE, auto_unbox = TRUE, digits = NA
)
message("Exported ", length(cases), " R regression cases to ", output)
