# Run from ibist-R. Python also validates independent test formulas.
for (file in c("asymcheck", "indep_test", "nominal_assoc", "ordinal_assoc",
               "cmh_ll_test",
               "assoc_test")) source(paste0("R/", file, ".R"))
examples <- list(
  list(x = rbind(c(12, 3, 1), c(4, 9, 2), c(1, 5, 13))),
  list(x = rbind(c(26, 26, 23, 18, 9), c(6, 7, 9, 14, 23))),
  list(x = matrix(c(40, 10, 15, 35), 2, byrow = TRUE)),
  list(x = matrix(c(1, 0, 0, 1), 2, byrow = TRUE)),
  list(x = matrix(c(8, 1, 4, 2, 5, 9), 2, byrow = TRUE),
       row_labels = c("1", "2"), col_labels = c("1", "2", "10")),
  list(x = rbind(c(2, 0, 1), c(0, 0, 0), c(1, 0, 4)))
)
set.seed(20261107)
for (i in 1:15) {
  nr <- sample(2:5, 1)
  nc <- sample(2:5, 1)
  examples[[length(examples) + 1L]] <- list(
    x = matrix(sample(1:30, nr * nc, TRUE), nr))
}
cases <- list()
fields <- c("statistic", "p.value", "estimate", "null.value", "parameter",
            "std.error", "chisq")
scalar <- function(fit, name) {
  if (is.null(fit[[name]])) NULL else unname(fit[[name]][1L])
}
for (example in examples) {
  tab <- example$x
  if (!is.null(example$row_labels)) rownames(tab) <- example$row_labels
  if (!is.null(example$col_labels)) colnames(tab) <- example$col_labels
  input <- example
  input$x <- lapply(seq_len(nrow(tab)), function(i) unname(tab[i, ]))
  for (method in c("multinomial", "pearson", "cmh")) {
    if (method == "pearson" &&
        (any(rowSums(tab) == 0) || any(colSums(tab) == 0))) next
    measures <- switch(method,
      multinomial = c("gamma", "tau-b", "tau-c", "somers-c|r",
                       "somers-r|c", "pearson", "spearman", "gamma"),
      pearson = c("phi", "contingency", "cramer-v", "phi"),
      cmh = c("pearson", "spearman", "pearson"))
    alternatives <- if (method == "pearson") "two.sided" else {
      c("two.sided", "greater", "less")
    }
    for (alternative in alternatives) {
      for (scored in if (method == "pearson") FALSE else c(FALSE, TRUE)) {
        for (correct in if (method == "pearson" && all(dim(tab) == 2)) {
          c(FALSE, TRUE)
        } else FALSE) {
          args <- list(measure = measures, method = method,
                       alternative = alternative, conf.level = 0.9,
                       correct = correct)
          if (scored) {
            args$row.scores <- seq_len(nrow(tab))^2
            args$col.scores <- seq_len(ncol(tab))^2
          }
          fits <- suppressWarnings(do.call(assoc.test, c(list(x = tab), args)))
          results <- unname(lapply(fits, function(fit) {
            values <- lapply(fields, function(name) scalar(fit, name))
            names(values) <- gsub("\\.", "_", fields)
            c(values, list(
              conf_int = if (is.null(fit$conf.int)) NULL else {
                as.list(as.numeric(fit$conf.int))
              }, conf_level = attr(fit$conf.int, "conf.level"),
              measure = fit$measure, method = fit$method,
              test_method = fit$test.method, alternative = fit$alternative))
          }))
          names(args) <- gsub("\\.", "_", names(args))
          args$measure <- as.list(args$measure)
          cases[[length(cases) + 1L]] <- c(input, args, list(results = results))
        }
      }
    }
  }
}
files <- paste0("R/", c("indep_test", "nominal_assoc", "ordinal_assoc",
                        "cmh_ll_test", "assoc_test"), ".R")
jsonlite::write_json(
  list(source = "ibist-R working sources; independent Python checks also used",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = as.list(tools::md5sum(files)), cases = cases),
  "../ibist-py/tests/reference/assoc_test_r.json", pretty = TRUE,
  auto_unbox = TRUE, digits = NA, na = "string", null = "null"
)
message("Exported ", length(cases), " association-test calls.")
