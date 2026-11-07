# Shared CI interface cases, generated from R first. Run from ibist-R.
files <- c("ci", "prop_1s_ci", "prop_2s_ci", "rr_2x2_ci",
           "or_2x2_ci", "rate_1s_ci", "rate_2s_ci", "or_paired_ci",
           "prop_paired_ci", "or_common")
for (file in files) source(paste0("R/", file, ".R"))
tab <- matrix(c(7, 27, 1, 33), 2, byrow = TRUE)
specs <- list(
  list(name = "prop_1s_ci", args = list(x = 5, n = 40),
       methods = c("wilson", "wald", "agresti-coull", "clopper-pearson",
                   "blaker")),
  list(name = "prop_2s_ci", args = list(x = c(12, 4), n = c(40, 40)),
       methods = c("newcombe", "wald", "miettinen-nurminen", "mee",
                   "agresti-caffo", "brown-li")),
  list(name = "rr_2x2_ci", args = list(x = tab),
       methods = c("koopman", "katz")),
  list(name = "or_2x2_ci", args = list(x = tab),
       methods = c("wald", "gart", "cornfield", "bp")),
  list(name = "rate_1s_ci", args = list(x = 5, exposure = 10),
       methods = c("exact", "score", "wh", "wald", "log")),
  list(name = "rate_2s_ci", args = list(x = c(9, 12), exposure = c(2, 7)),
       methods = c("log", "score", "exact"))
)
cases <- list()
for (spec in specs) {
  fun <- get(gsub("_", ".", spec$name))
  for (level in c(0.2, 0.8, 0.95, 0.99)) {
    for (side in c("two.sided", "lower", "upper")) {
      for (method in spec$methods) {
        options <- c(spec$args, list(conf.level = level, side = side,
                                     method = method))
        result <- do.call(fun, options)
        names(options) <- gsub("[.]", "_", names(options))
        if (is.matrix(options$x)) {
          options$x <- lapply(1:2, function(i) unname(options$x[i, ]))
        }
        estimate <- if (spec$name == "rate_2s_ci") {
          result$estimate["rate ratio"]
        } else result$estimate
        cases[[length(cases) + 1L]] <- list(
          function_name = spec$name, options = options,
          side = result$side, estimate = unname(estimate),
          conf_level = result$conf.level, methods = as.list(result$method),
          conf_int = lapply(seq_len(nrow(result$conf.int)), function(i) {
            unname(result$conf.int[i, ])
          })
        )
      }
    }
  }
}
# Compare every public CI entry point's default metadata as well as limits.
default.specs <- c(specs, list(
  list(name = "or_paired_ci", args = list(x = tab)),
  list(name = "prop_paired_ci", args = list(b = 8, c = 25, n = 180)),
  list(name = "or_common_ci", args = list(x = list(tab, tab)))
))
defaults <- lapply(default.specs, function(spec) {
  result <- do.call(get(gsub("_", ".", spec$name)), spec$args)
  options <- spec$args
  rows <- function(tab) lapply(1:2, function(i) unname(tab[i, ]))
  if (is.matrix(options$x)) options$x <- rows(options$x)
  if (spec$name == "or_common_ci") options$x <- lapply(options$x, rows)
  list(function_name = spec$name, options = options,
       method = unname(result$method), side = result$side,
       conf_level = result$conf.level, parameter = result$parameter,
       null_value = if (is.null(result$null.value)) NULL else {
         unname(result$null.value)
       })
})
jsonlite::write_json(
  list(source = "R CI interface; supplemental to independent formula tests",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(paste0("R/", files, ".R"))),
       defaults = defaults, cases = cases),
  "../ibist-py/tests/reference/ci_interface_r.json", pretty = TRUE,
  auto_unbox = TRUE, digits = NA, na = "string", null = "null"
)
message("Exported ", length(cases), " CI interface cases.")
