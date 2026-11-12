# Export shared prepared data and odds-ratio references for IBist.jl.
# Run from ibist-R: Rscript data-raw/export-julia.R [../IBist.jl]
# This reuses the portable exporters; it never repeats data-cleaning rules.

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "../IBist.jl"
rscript <- file.path(R.home("bin"), "Rscript")
run_exporter <- function(script, target) {
  status <- system2(rscript, c("--vanilla", shQuote(script), shQuote(target)))
  if (status != 0L) stop("Exporter failed: ", script)
}
run_exporter("data-raw/export-python.R", file.path(output, "data"))
run_exporter("data-raw/export-python-reference.R",
             file.path(output, "test", "reference", "or_2x2_r.json"))

# One-sided and low-confidence comparisons supplement the shared fixtures.
source("R/ci.R")
source("R/or_2x2_ci.R")
tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)
cases <- list()
number <- function(x) {
  if (is.infinite(x)) "inf" else if (is.na(x)) "nan" else unname(x)
}
for (side in c("lower", "upper")) {
  for (level in c(0.2, 0.95)) {
    for (midp in c(FALSE, TRUE)) {
      methods <- c("wald", "gart", "cornfield", "bp", "cornfield")
      result <- or.2x2.ci(tab, method = methods, conf.level = level,
                          midp = midp, side = side)
      cases[[length(cases) + 1L]] <- list(
        table = lapply(1:2, function(i) unname(tab[i, ])),
        side = side, conf_level = level, midp = midp,
        estimate = number(result$estimate), methods = I(methods),
        labels = I(result$method),
        conf_int = lapply(seq_len(nrow(result$conf.int)), function(i) {
          unname(lapply(result$conf.int[i, ], number))
        })
      )
    }
  }
}
jsonlite::write_json(
  list(source = "R one-sided OR limits; supplemental regression checks",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(c("R/or_2x2_ci.R", "R/ci.R"))),
       cases = cases),
  file.path(output, "test", "reference", "or_sides_r.json"), pretty = TRUE,
  auto_unbox = TRUE, digits = NA
)
message("Exported ", length(cases), " one-sided OR reference calls.")
