# Book examples and threshold boundaries; no inference is run.
files <- paste0("R/", c("ci", "asymcheck", "asymcheck_book", "prop_1s_ci",
  "prop_2s_ci", "rr_2x2_ci", "or_2x2_ci", "rate_1s_ci", "rate_2s_ci",
  "rate_test", "indep_test"), ".R")
for (file in files) source(file)
cases <- list()
python.names <- function(x) {
  if (is.list(x)) {
    for (i in seq_along(x)) x[i] <- list(python.names(x[[i]]))
    if (!is.null(names(x))) names(x) <- gsub("\\.", "_", names(x))
  }
  x
}
add <- function(kind, args) {
  fun <- get(switch(kind, p1 = "prop.1s.ci", p2 = "prop.2s.ci",
    rr = "rr.2x2.ci", or = "or.2x2.ci", rate1 = "rate.1s.ci",
    rate2 = "rate.2s.ci", `rate-test` = "rate.test"))
  check <- do.call(asymcheck, c(list(object = fun), args))
  for (name in intersect(c("input", "sample"), names(check$checks))) {
    check$checks[[name]] <- check$checks[[name]] - 1L
  }
  if (is.matrix(args$x)) args$x <-
    lapply(seq_len(nrow(args$x)), function(i) unname(args$x[i, ]))
  cases[[length(cases)+1L]] <<- list(kind = kind, args = python.names(args),
    result = python.names(unclass(check)))
}
add("p1", list(x = c(0, 16), n = 40, correct = TRUE, midp = TRUE,
  method = c("wald", "wilson", "agresti-coull", "clopper-pearson", "wald"),
  side = "lower"))
add("p1", list(x = 10, n = 20))
vf <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
add("p2", list(x = vf, method = c("newcombe", "wald", "miettinen-nurminen",
  "mee", "agresti-caffo", "brown-li", "chan-zhang"), side = "upper"))
add("rr", list(x = vf, method = c("katz", "koopman")))
add("or", list(x = vf, method = c("wald", "gart", "cornfield", "bp"),
  midp = TRUE, side = "upper", tol = 1e-7))
add("or", list(x = matrix(2, 2, 2), method = "gart"))
add("rate2", list(x = c(151, 55), exposure = c(57518.1, 74573.5),
  method = c("log", "score", "exact")))
add("rate2", list(x = c(10, 11)))
add("rate-test", list(x = 411, exposure = 25800, r = 0.0119))
add("rate-test", list(x = c(9, 12), exposure = c(1817.6, 7496.3)))
for (x in c(49, 50)) add("rate1", list(x = x, method = "wald"))
for (x in c(9, 10)) add("rate1", list(x = x, method = "log"))
add("rate1", list(x = 3, method = c("wh", "score", "exact")))
json <- jsonlite::toJSON(list(r_revision = system2("git",
    c("rev-parse", "HEAD"), stdout = TRUE),
  source_md5 = as.list(tools::md5sum(files)), cases = cases),
  auto_unbox = TRUE, null = "null", na = "null", digits = NA)
writeLines(sub("\n$", "", jsonlite::prettify(json, indent = 1)),
  "../ibist-py/tests/reference/book_asymcheck_r.json")
