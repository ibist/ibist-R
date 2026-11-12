# Representative interface cases; independent formulas are checked in Python.
for (file in c("or_common", "ca_trend_test", "trend_test", "power_p1s",
               "power_p2s")) source(paste0("R/", file, ".R"))
rows <- function(tab) lapply(seq_len(nrow(tab)), function(i) unname(tab[i, ]))
cases <- list()
add <- function(kind, args) {
  fun <- switch(kind, trend = trend.test, one = power.p1s.test,
                two = power.p2s.test)
  result <- do.call(fun, args)
  if (kind == "trend") {
    args$x <- rows(args$x)
    result <- list(statistic = unname(result$statistic),
                   chisq = unname(result$chisq),
                   score_statistic = unname(result$score.statistic),
                   p_value = result$p.value, method = result$method,
                   alternative = result$alternative)
  } else {
    result <- unclass(result)
    names(result) <- gsub("\\.", "_", names(result))
  }
  names(args) <- gsub("\\.", "_", names(args))
  cases[[length(cases) + 1L]] <<- list(kind = kind, args = args,
                                     result = result)
}
add("trend", list(x = matrix(c(106, 109, 128, 130, 99, 95, 76, 74), 4)))
add("trend", list(x = rbind(c(8, 12), c(12, 8), c(16, 4)),
                  scores = c(1, 2, 4), alternative = "greater"))
for (alt in c("two.sided", "greater", "less")) {
  add("trend", list(x = rbind(c(2, 0), c(1, 1), c(0, 2)),
                    scores = c(0, 1, 3), exact = TRUE, alternative = alt))
}
add("trend", list(x = rbind(c(8, 12), c(0, 0), c(16, 4))))
add("one", list(n = 50, p0 = 0.1, p1 = 0.25))
add("one", list(n = 40, p0 = 0.4, p1 = 0.2,
                 alternative = "less", correct = TRUE))
add("one", list(n = 40, p0 = 0.3, p1 = 0.5,
                 strict = TRUE, correct = TRUE))
add("one", list(p0 = 0.1, p1 = 0.25, power = 0.8))
add("one", list(n = 50, p1 = 0.25, power = 0.8))
add("one", list(n = 50, p0 = 0.25, power = 0.8, alternative = "less"))
add("one", list(n = 50, p0 = 0.1, p1 = 0.25, sig.level = NULL, power = 0.8))
for (method in c("quantile", "cp", "midp", "midp-rand")) {
  add("one", list(n = 35, p0 = 0.1, p1 = 0.25, exact = TRUE,
                   exact.method = method))
}
add("one", list(n = 35, p0 = 0.9, p1 = 0.75, exact = TRUE,
                 exact.method = "midp-rand", alternative = "less"))
add("one", list(p0 = 0.01, p1 = 0.99, power = 0.8, exact = TRUE,
                 alternative = "greater", size.rule = "minimal", max_n = 100))
add("one", list(p0 = 0.1, p1 = 0.25, power = 0.8, exact = TRUE, max_n = 500))
add("one", list(n = 50, p0 = 0.1, power = 0.8, exact = TRUE,
                 exact.method = "cp"))
add("one", list(n = 50, p1 = 0.25, power = 0.8, exact = TRUE,
                 exact.method = "cp"))
add("one", list(n = 50, p0 = 0.1, p1 = 0.25, power = 0.8,
                 sig.level = NULL, exact = TRUE, exact.method = "cp"))
add("two", list(n = 100, p1 = 0.3, p2 = 0.5, group.rate = 2))
add("two", list(n = 100, p1 = 0.3, p2 = 0.5, group.rate = 1.5,
                 correct = TRUE, strict = TRUE, alternative = "one.sided"))
add("two", list(n = 100, p1 = 0.3, p2 = 0.5, group.rate = 1.5,
                 correct = TRUE, strict = TRUE))
add("two", list(p1 = 0.3, p2 = 0.5, power = 0.8, group.rate = 1.5))
add("two", list(n = 100, p2 = 0.5, power = 0.8))
add("two", list(n = 100, p1 = 0.3, power = 0.8))
add("two", list(n = 100, p1 = 0.3, p2 = 0.5, group.rate = 2,
                 sig.level = NULL, power = 0.8))
files <- paste0("R/", c("ca_trend_test", "trend_test", "power_p1s",
                        "power_p2s"), ".R")
output <- "../ibist-py/tests/reference/trend_power_r.json"
jsonlite::write_json(
  list(source = "ibist-R working sources; independent Python checks also used",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(files)), cases = cases),
  output, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null")
# Wrap long property values without changing JSON string contents.
lines <- readLines(output)
for (i in which(nchar(lines) > 80)) {
  line <- lines[i]
  colon <- regexpr(": ", line, fixed = TRUE)[1L]
  value <- substring(line, colon + 2L)
  indent <- min(nchar(line) - nchar(trimws(line)) + 2L, 80L - nchar(value))
  lines[i] <- paste0(substr(line, 1L, colon), "\n",
                      strrep(" ", max(0L, indent)), value)
}
writeLines(lines, output)
message("Exported ", length(cases), " trend/power cases.")
