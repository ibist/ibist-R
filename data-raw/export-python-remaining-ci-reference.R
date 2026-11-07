# Native Python CI regression references; run from the R package root.
# Uses the working R sources and records their hashes, not only Git HEAD.
pkgload::load_all(".", quiet = TRUE)
sources <- c("R/ci.R", "R/prop_paired_ci.R", "R/or_paired_ci.R",
             "R/or_common.R", "R/prop_1s_ci.R", "R/prop_2s_ci.R",
             "R/prop_2s_exact.R", "src/wang_paired_ci.cpp", "src/wang_core.h")
cases <- list()
add <- function(name, options) {
  r.options <- options
  names(r.options) <- gsub("_", ".", names(r.options), fixed = TRUE)
  if (!is.null(r.options$method)) r.options$method <- unlist(r.options$method)
  result <- suppressWarnings(do.call(get(gsub("_", ".", name)), r.options))
  portable <- options
  rows <- function(tab) lapply(1:2, function(i) unname(tab[i, ]))
  if (is.matrix(portable$x)) portable$x <- rows(portable$x)
  if (name == "or_common_ci") portable$x <- lapply(portable$x, rows)
  cases[[length(cases) + 1L]] <<- list(
    function_name = name, options = portable,
    estimate = unname(result$estimate),
    side = result$side, methods = as.list(result$method),
    conf_int = lapply(seq_len(nrow(result$conf.int)), function(i) {
      unname(result$conf.int[i, ])
    })
  )
}
for (counts in list(c(8, 25, 180), c(160, 3, 320), c(7, 3, 20),
                    c(0, 0, 8), c(0, 8, 8), c(8, 0, 8), c(2, 1, 3))) {
  for (level in c(0.2, 0.8, 0.95, 0.99)) {
    for (side in c("two.sided", "lower", "upper")) {
      add("prop_paired_ci", list(b = counts[1], c = counts[2], n = counts[3],
          conf_level = level, side = side,
          method = as.list(c("score", "wald", "waldcc", "agresti-min"))))
    }
  }
}
# Every collapsed table at n=3, plus the n=4 regression examples.
wang <- expand.grid(b = 0:3, c = 0:3)
wang <- wang[rowSums(wang) <= 3, ]
for (i in seq_len(nrow(wang))) {
  for (side in c("two.sided", "lower", "upper")) {
    add("prop_paired_ci", list(b = wang$b[i], c = wang$c[i], n = 3,
        method = "wang", precision = 0.0001, side = side))
  }
}
for (counts in list(c(3, 0, 4), c(2, 2, 4), c(1, 2, 4))) {
  add("prop_paired_ci", list(b = counts[1], c = counts[2], n = counts[3],
      method = "wang", precision = 0.0001))
}
for (counts in list(c(51, 30, 10, 41), c(5, 4, 0, 6), c(5, 0, 4, 6),
                    c(5, 0, 0, 6), c(0, 2, 1, 0), c(1, 1, 1, 1))) {
  methods <- c("wilson", "clopper-pearson", "midp", "blaker")
  if (all(counts[2:3] > 0)) methods <- c("wald", methods)
  for (level in c(0.2, 0.8, 0.95, 0.99)) {
    for (side in c("two.sided", "lower", "upper")) {
      add("or_paired_ci", list(x = matrix(counts, 2, byrow = TRUE),
          conf_level = level, side = side, method = as.list(methods)))
    }
  }
}
examples <- list(list(c(7, 3, 2, 8), c(4, 6, 1, 9)),
                 list(c(2, 3, 1, 4)), list(c(0, 3, 2, 0), c(0, 4, 1, 0)),
                 list(c(3, 0, 0, 2)), list(c(0, 0, 0, 0)),
                 list(c(0, 0, 4, 2), c(3, 4, 1, 2)))
set.seed(20261107)
for (i in 1:8) {
  examples[[length(examples) + 1L]] <- replicate(
    sample(1:4, 1), sample(0:10, 4, TRUE), simplify = FALSE)
}
for (example in examples) {
  for (level in c(0.2, 0.8, 0.95, 0.99)) {
    for (side in c("two.sided", "lower", "upper")) {
      for (algorithm in c("network", "convolution")) {
        add("or_common_ci", list(x = lapply(example, matrix, nrow = 2,
                                           byrow = TRUE),
            conf_level = level, side = side, algorithm = algorithm))
      }
    }
  }
}
for (case in list(list(x = c(4, 1), n = c(6, 5)),
                 list(x = c(0, 0), n = c(3, 5)),
                 list(x = c(3, 5), n = c(3, 5)),
                 list(x = c(3, 0), n = c(3, 5)),
                 list(x = c(0, 5), n = c(3, 5)),
                 list(x = c(9, 3), n = c(10, 10)))) {
  for (side in c("two.sided", "lower", "upper")) {
    add("prop_2s_ci", c(case, list(method = "chan-zhang", side = side)))
  }
}
add("prop_2s_ci", list(x = c(4, 1), n = c(6, 5), method = "chan-zhang",
                        conf_level = 0.2, side = "lower"))
add("prop_2s_ci", list(x = c(5, 10), n = c(9, 19), method = "chan-zhang",
                        conf_level = 0.12))
jsonlite::write_json(
  list(source = "ibist-R working sources; independent checks also required",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(sources)), cases = cases),
  "../ibist-py/tests/reference/remaining_ci_r.json", pretty = TRUE,
  auto_unbox = TRUE, digits = NA, na = "string"
)
message("Exported ", length(cases), " remaining CI cases.")
