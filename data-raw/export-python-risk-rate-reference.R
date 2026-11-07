# Regenerate risk-ratio and rate regression references from the R sources.
# Run from the R package root; optional argument selects the output directory.
# Rscript data-raw/export-python-risk-rate-reference.R
# Python tests also validate independent formulas and likelihood fits.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite to export development reference cases.")
}
for (file in c("ci", "prop_1s_ci", "prop_2s_ci", "rr_2x2_ci",
               "rate_1s_ci", "rate_2s_ci", "rate_test")) {
  source(paste0("R/", file, ".R"))
}
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "../ibist-py/tests/reference"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
sources <- c("R/ci.R", "R/prop_1s_ci.R", "R/prop_2s_ci.R",
             "R/rr_2x2_ci.R", "R/rate_1s_ci.R", "R/rate_2s_ci.R",
             "R/rate_test.R")
write_cases <- function(cases, file) {
  jsonlite::write_json(
    list(source = "ibist-R source; supplemental to independent Python checks",
         r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
         source_md5 = as.list(tools::md5sum(sources)), cases = cases),
    file.path(output, file), pretty = TRUE, auto_unbox = TRUE, digits = NA,
    na = "string"
  )
  message("Exported ", length(cases), " cases to ", file)
}
limits <- function(result) {
  lapply(seq_len(nrow(result$conf.int)), function(i) {
    unname(result$conf.int[i, ])
  })
}

risk <- list()
examples <- list(c(38, 62, 267, 717), c(11, 4, 2, 6), c(4, 11, 6, 2),
                 c(0, 3, 0, 5), c(0, 3, 5, 0), c(3, 0, 0, 5),
                 c(3, 0, 5, 0), c(2, 1, 4, 1))
set.seed(20261107)
for (i in 1:20) examples[[length(examples) + 1L]] <- sample(1:100, 4, TRUE)
for (values in examples) {
  tab <- matrix(values, 2, byrow = TRUE)
  methods <- if (all(tab[, 1] > 0)) c("koopman", "katz") else "koopman"
  for (level in c(0.2, 0.8, 0.95, 0.99)) {
    for (alternative in c("two.sided", "less", "greater")) {
      result <- rr.2x2.ci(tab, conf.level = level, method = methods,
                          alternative = alternative)
      risk[[length(risk) + 1L]] <- list(
        x = lapply(1:2, function(i) unname(tab[i, ])),
        options = list(conf_level = level, method = as.list(methods),
                       alternative = alternative),
        estimate = unname(result$estimate), conf_int = limits(result),
        sample_proportions = unname(result$sample.proportions)
      )
    }
  }
}
write_cases(risk, "rr_2x2_r.json")

rates <- list()
for (x in c(0, 1, 5, 411)) {
  for (exposure in c(1, 10, 25800)) {
    for (level in c(0.8, 0.95, 0.99)) {
      for (correct in c(FALSE, TRUE)) {
        methods <- c("exact", "score", "wh", "wald", "log")
        result <- rate.1s.ci(x, exposure, conf.level = level,
                             method = methods, correct = correct)
        rates[[length(rates) + 1L]] <- list(
          function_name = "rate_1s_ci", x = x, exposure = exposure,
          options = list(conf_level = level, method = as.list(methods),
                         correct = correct),
          estimate = unname(result$estimate), conf_int = limits(result),
          negative_wh_lower = x == 0 && correct
        )
      }
    }
  }
}
examples <- list(list(x = c(151, 55), exposure = c(57518.1, 74573.5)),
                 list(x = c(9, 12), exposure = c(1817.6, 7496.3)),
                 list(x = c(0, 2), exposure = c(1, 3)),
                 list(x = c(2, 0), exposure = c(3, 1)),
                 list(x = c(12, 5), exposure = c(100, 80)))
for (case in examples) {
  methods <- if (all(case$x > 0)) c("log", "score", "exact") else {
    c("score", "exact")
  }
  for (level in c(0.8, 0.95, 0.99)) {
    result <- rate.2s.ci(case$x, case$exposure, conf.level = level,
                         method = methods)
    rates[[length(rates) + 1L]] <- list(
      function_name = "rate_2s_ci", x = case$x, exposure = case$exposure,
      options = list(conf_level = level, method = as.list(methods)),
      estimate = unname(result$estimate["rate ratio"]),
      conf_int = limits(result)
    )
  }
}
tests <- c(examples[c(1, 2, 5)],
           list(list(x = 411, exposure = 25800, r = 0.0119),
                list(x = 0, exposure = 10, r = 0.1),
                list(x = 1, exposure = 1, r = 1.1)))
for (case in tests) {
  for (correct in c(FALSE, TRUE)) {
    for (alternative in c("two.sided", "less", "greater")) {
      r <- if (is.null(case$r)) 1 else case$r
      result <- suppressWarnings(rate.test(case$x, case$exposure, r = r,
                                            correct = correct,
                                            alternative = alternative))
      estimate <- if (length(case$x) == 1L) result$estimate else {
        result$estimate["rate ratio"]
      }
      rates[[length(rates) + 1L]] <- list(
        function_name = "rate_test", x = case$x, exposure = case$exposure,
        options = list(r = r, correct = correct, alternative = alternative),
        estimate = unname(estimate),
        statistic = unname(result$statistic), p_value = result$p.value,
        conf_int = unname(result$conf.int)
      )
    }
  }
}
write_cases(rates, "rates_r.json")
