# Explicitly regenerate R comparison cases, then review changes before use.
# These are cross-language regression values, not independent validation.
# Run from the R package root: Rscript data-raw/export-python-reference.R [output]

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite to export development reference cases.")
}
source("R/ci.R")
source("R/or_2x2_ci.R")
cases <- list(
  fagerland = c(7, 27, 1, 33),
  balanced = c(10, 10, 10, 10),
  fluorosis = c(38, 267, 62, 717),
  sparse = c(1, 4, 3, 2),
  zero_a = c(0, 5, 4, 6),
  zero_b = c(4, 0, 2, 5),
  zero_c = c(4, 2, 0, 5),
  zero_d = c(4, 2, 5, 0),
  diagonal = c(5, 0, 0, 6),
  antidiagonal = c(0, 5, 6, 0)
)
references <- list()
for (name in names(cases)) {
  tab <- matrix(cases[[name]], nrow = 2, byrow = TRUE)
  for (level in c(0.90, 0.95, 0.99)) {
    for (midp in c(FALSE, TRUE)) {
      methods <- c(if (all(tab > 0)) "wald", "gart", "cornfield", "bp")
      result <- or.2x2.ci(tab, method = methods, conf.level = level,
                          midp = midp)
      # JSON strings preserve unbounded limits without nonstandard Infinity.
      number <- function(x) {
        if (is.infinite(x)) "inf" else if (is.na(x)) "nan" else x
      }
      intervals <- lapply(seq_along(result$method), function(i) {
        list(method = result$method[i], lower = number(result$conf.int[i, 1]),
             upper = number(result$conf.int[i, 2]))
      })
      references[[length(references) + 1L]] <- list(
        name = name, table = unname(split(tab, row(tab))), conf_level = level,
        midp = midp, estimate = number(unname(result$estimate)),
        intervals = intervals
      )
    }
  }
}
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else {
  "../ibist-py/tests/reference/or_2x2_r.json"
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  list(source = "R source in ibist-R; independently checked in Python tests",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(c("R/or_2x2_ci.R", "R/ci.R"))),
       cases = references), output, pretty = TRUE, auto_unbox = TRUE,
  digits = NA
)
message("Exported ", length(references), " R comparison cases to ", output)
