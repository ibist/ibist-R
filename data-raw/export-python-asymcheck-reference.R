# Representative diagnostics only; no inference or large Cartesian grids.
for (file in c("ci", "asymcheck", "indep_test", "rate_1s_ci",
               "mantelhaen_asymcheck")) source(paste0("R/", file, ".R"))
cases <- list()
rows <- function(x) lapply(seq_len(nrow(x)), function(i) unname(x[i, ]))
add <- function(kind, args) {
  check <- if (kind == "mh") .asymcheck_mh(args$x) else {
    fun <- if (kind == "indep") indep.test else rate.1s.ci
    do.call(asymcheck, c(list(object = fun), args))
  }
  checks <- check$checks
  names(checks) <- gsub("\\.", "_", names(checks))
  result <- list(procedure = check$procedure, method = check$method,
    applicable = check$applicable, condition_met = check$condition.met,
    checks = checks, notes = check$notes)
  if (kind == "indep") {
    args$x <- rows(args$x)
    result$expected <- rows(check$details$expected)
  }
  if (kind == "mh") {
    args$x <- lapply(args$x, rows)
    result$variance_sum <- check$details$variance.sum
    result$variance <- unname(check$details$variance)
  }
  names(args) <- gsub("\\.", "_", names(args))
  cases[[length(cases) + 1L]] <<- list(kind = kind, args = args,
                                     result = result)
}
tab <- matrix(c(12, 8, 5, 15), 2, byrow = TRUE)
for (method in c("pearson", "lr", "fisher")) {
  add("indep", list(x = tab, method = method))
}
add("indep", list(x = rbind(c(1, 5, 5, 5, 5), c(1, 5, 5, 5, 5))))
add("indep", list(x = matrix(c(0, 1, 1, 20), 2)))
add("mh", list(x = list(matrix(c(88, 469, 48, 653), 2),
                         matrix(c(27, 480, 34, 1324), 2))))
add("mh", list(x = list(matrix(1, 2, 2))))
add("mh", list(x = list(matrix(0, 2, 2))))
add("rate", list(x = 12, exposure = 100, method = "wald"))
add("rate", list(x = 12, exposure = 100, method = c("exact", "wald", "wald"),
                  conf.level = 0.99, side = "lower", correct = FALSE))
add("rate", list(x = 0, method = "exact"))
add("rate", list(x = 0, method = "log"))
add("rate", list(x = 0, method = "wald", correct = FALSE))
add("rate", list(x = 0, method = "wald", correct = TRUE))
path <- "R/asymcheck.R"
json <- jsonlite::toJSON(list(source = path,
  source_md5 = unname(tools::md5sum(path)), cases = cases),
  auto_unbox = TRUE, null = "null", na = "null", digits = NA)
json <- sub("\n$", "", jsonlite::prettify(json, indent = 1))
writeLines(json,
  "../ibist-py/tests/reference/asymcheck_r.json")
