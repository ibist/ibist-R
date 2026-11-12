# Run from ibist-R; Python also checks formulas independently.
source("R/ordinal_assoc.R")
examples <- list(
  list(x = rbind(c(12, 3, 1), c(4, 9, 2), c(1, 5, 13))),
  list(x = rbind(c(26, 26, 23, 18, 9), c(6, 7, 9, 14, 23)),
       row_labels = c("No", "Yes"), col_labels = as.character(0:4)),
  list(x = rbind(c(8, 1, 4), c(2, 5, 9)),
       row_labels = c("1", "2"), col_labels = c("1", "2", "10")),
  list(x = rbind(c(8, 1, 4), c(2, 5, 9)),
       row_scores = c(-1, 4), col_scores = c(0, 2, 10)),
  list(x = matrix(c(2, 1, 3, 4, 5, 0, 1, 2), 2, byrow = TRUE),
       row_scores = c(2, -1), col_scores = c(4, 0, 1, 3)),
  list(x = diag(c(1, 2, 3))), list(x = diag(c(4, 4))),
  list(x = matrix(c(2, 3, 3, 2), 2, byrow = TRUE)),
  list(x = matrix(c(1, 0, 0, 1), 2, byrow = TRUE)),
  list(x = rbind(c(3, 0, 1), c(0, 0, 0), c(1, 0, 4))),
  list(x = rbind(c(3, 0, 1), c(0, 0, 0), c(1, 0, 4)),
       row_scores = c(0, 2, 10), col_scores = c(1, 3, 4)),
  list(x = rbind(c(3, 0, 1), c(0, 0, 0), c(1, 0, 4)),
       row_labels = c("0", "2", "10"), col_labels = c("1", "3", "4")),
  list(x = matrix(10, 3, 4)),
  list(x = rbind(c(1, 1, 0), c(0, 1, 1)),
       measure = c("gamma", "spearman", "gamma"))
)
set.seed(20261107)
for (i in 1:60) {
  nr <- sample(2:6, 1)
  nc <- sample(2:6, 1)
  tab <- matrix(sample(1:30, nr * nc, TRUE), nr)
  examples[[length(examples) + 1L]] <- list(x = tab)
}
cases <- list()
for (example in examples) {
  tab <- example$x
  if (!is.null(example$row_labels)) rownames(tab) <- example$row_labels
  if (!is.null(example$col_labels)) colnames(tab) <- example$col_labels
  for (alternative in c("two.sided", "greater", "less")) {
    level <- if (alternative == "two.sided") 0.95 else 0.9
    args <- list(x = tab, conf.level = level, alternative = alternative)
    if (!is.null(example$row_scores)) args$row.scores <- example$row_scores
    if (!is.null(example$col_scores)) args$col.scores <- example$col_scores
    if (!is.null(example$measure)) args$measure <- example$measure
    result <- do.call(ordinal.assoc, args)
    example$x <- lapply(seq_len(nrow(tab)), function(i) unname(tab[i, ]))
    cases[[length(cases) + 1L]] <- c(
      example, list(
        conf_level = level, alternative = alternative,
        measures = unname(as.list(result$measure)),
        estimate = unname(as.list(result$estimate)),
        std_error = unname(as.list(result$std.error)),
        lower = unname(as.list(result$conf.int[, 1])),
        upper = unname(as.list(result$conf.int[, 2])),
        statistic = unname(as.list(result$statistic)),
        p_value = unname(as.list(result$p.value)), S = result$S,
        var_S = result$var.S, concordant = result$concordant,
        discordant = result$discordant, n = result$n, method = result$method))
  }
}
jsonlite::write_json(
  list(source = "ibist-R; independent Python validation also required",
       r_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_md5 = unname(tools::md5sum("R/ordinal_assoc.R")), cases = cases),
  "../ibist-py/tests/reference/ordinal_r.json", pretty = TRUE,
  auto_unbox = TRUE, digits = NA, na = "string"
)
message("Exported ", length(cases), " ordinal association cases.")
