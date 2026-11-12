# Explicitly regenerate cross-language regression cases, then review changes.
# Independent formula and published-result tests are maintained in Python.
# Run from the R package root:
# Rscript data-raw/export-python-contingency-reference.R

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite to export development reference cases.")
}
source("R/asymcheck.R")
source("R/indep_test.R")
source("R/nominal_assoc.R")
source("R/or_2x2_ci.R")
source("R/or_common.R")
rows <- function(x) lapply(seq_len(nrow(x)), function(i) unname(x[i, ]))
number <- function(x) {
  if (is.infinite(x)) "inf" else if (is.na(x)) "nan" else unname(x)
}
cases <- list()
add <- function(case) cases[[length(cases) + 1L]] <<- case
set.seed(20261106)
tables <- list(
  screening = matrix(c(1366, 14205, 383, 4983), 2, byrow = TRUE),
  tumor = matrix(c(109, 68, 14, 72, 79, 37, 9, 42, 139), 3, byrow = TRUE),
  zero = matrix(c(0, 20, 10, 10), 2, byrow = TRUE),
  capped = matrix(c(11, 10, 10, 10), 2, byrow = TRUE)
)
for (i in 1:40) {
  nr <- sample(2:5, 1); nc <- sample(2:5, 1)
  tables[[paste0("random", i)]] <- matrix(sample(1:100, nr * nc, TRUE), nr)
}
suppressWarnings({
  for (name in names(tables)) {
    tab <- tables[[name]]
    for (method in c("pearson", "lr")) {
      for (correct in c(FALSE, if (method == "pearson" && nrow(tab) == 2 &&
                                   ncol(tab) == 2) TRUE)) {
        result <- indep.test(tab, method = method, correct = correct)
        add(list(kind = "independence", name = name, table = rows(tab),
                 method = method, options = list(correct = correct),
                 statistic = unname(result$statistic), p_value = result$p.value,
                 expected = rows(result$expected),
                 residuals = rows(result$residuals),
                 stdres = rows(result$stdres)))
      }
    }
    result <- nominal.assoc(tab)
    add(list(kind = "nominal", name = name, table = rows(tab),
             estimate = unname(result$estimate)))
  }
  for (tab in list(matrix(c(7, 77, 2, 90), 2, byrow = TRUE),
                    matrix(c(0, 10, 3, 7), 2, byrow = TRUE),
                    matrix(c(12, 8, 5, 15), 2, byrow = TRUE))) {
    for (alternative in c("two.sided", "less", "greater")) {
      result <- indep.test(tab, method = "fisher", alternative = alternative)
      add(list(kind = "independence", table = rows(tab), method = "fisher",
               options = list(alternative = alternative),
               p_value = result$p.value,
               estimate = number(result$estimate),
               conf_int = lapply(result$conf.int, number)))
    }
  }
  for (tab in list(matrix(c(2, 1, 0, 0, 1, 2), 2, byrow = TRUE),
                    diag(2, 3), matrix(c(1, 2, 1, 0, 1, 2, 2, 1, 1), 3))) {
    result <- indep.test(tab, method = "fisher")
    add(list(kind = "independence", table = rows(tab), method = "fisher",
             p_value = result$p.value))
  }
  for (i in 1:8) {
    repeat {
      tab <- matrix(sample(0:2, 9, TRUE), 3)
      if (all(rowSums(tab) > 0) && all(colSums(tab) > 0)) break
    }
    result <- indep.test(tab, method = "fisher")
    add(list(kind = "independence", table = rows(tab), method = "fisher",
             p_value = result$p.value))
  }
  strata <- list(
    book = array(c(88, 469, 48, 653, 27, 480, 34, 1324), c(2, 2, 2)),
    unit = array(c(30, 20, 20, 30, 20, 30, 30, 20), c(2, 2, 2)),
    zero = array(c(0, 20, 10, 30, 5, 25, 15, 35), c(2, 2, 2))
  )
  for (i in 1:30) {
    k <- sample(2:5, 1)
    strata[[paste0("random", i)]] <-
      array(sample(1:100, 4 * k, TRUE), c(2, 2, k))
  }
  for (name in names(strata)) {
    x <- strata[[name]]
    for (method in c("breslow-day", "woolf")) {
      for (correct in c(FALSE, if (method == "breslow-day") TRUE)) {
        result <- or.homogeneity.test(x, method = method, correct = correct)
        add(list(kind = "homogeneity", name = name,
                 tables = lapply(seq_len(dim(x)[3]),
                                  function(i) rows(x[, , i])),
                 method = method, correct = correct,
                 statistic = unname(result$statistic),
                 p_value = result$p.value))
      }
    }
  }
  for (tabs in list(list(matrix(c(2, 3, 1, 4), 2, byrow = TRUE),
                         matrix(c(1, 2, 3, 2), 2, byrow = TRUE)),
                    list(matrix(c(0, 2, 1, 1), 2, byrow = TRUE)))) {
    result <- or.homogeneity.test(tabs)
    add(list(kind = "homogeneity", tables = lapply(tabs, rows),
             method = "zelen",
             statistic = unname(result$statistic), p_value = result$p.value))
  }
})
files <- c("R/indep_test.R", "R/nominal_assoc.R", "R/or_common.R")
output <- "../ibist-py/tests/reference/contingency_r.json"
jsonlite::write_json(
  list(source = "R source in ibist-R; independently checked in Python tests",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(files)), cases = cases), output,
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)
message("Exported ", length(cases), " R regression cases to ", output)
