# R is the shared preparation source; never regenerate published expectations.
for (file in c("asymcheck", "glm_gof")) source(paste0("R/", file, ".R"))
root <- "../ibist-py/tests/reference"
dir.create("tests/testthat/fixtures", recursive = TRUE, showWarnings = FALSE)

# Reproduce SAS LOGIEX19's training sample from the Wooldridge Mroz data.
# Split multiplication keeps the legacy SAS RANUNI recurrence exact in doubles.
env <- new.env()
utils::data("mroz", package = "wooldridge", envir = env)
seed <- 3939
keep <- logical(nrow(env$mroz))
for (i in seq_along(keep)) {
  hi <- seed %/% 65536
  lo <- seed %% 65536
  seed <- ((((hi*397204094) %% 2147483647)*65536) +
            lo*397204094) %% 2147483647
  keep[i] <- seed/2147483647 >= 0.3
}
stopifnot(sum(keep) == 552L)
raw <- env$mroz[keep, ]
training <- data.frame(observed = 1-raw$inlf,
  income = signif(raw$nwifeinc, 7), educ = raw$educ, exper = raw$exper,
  exper2 = raw$expersq, age = raw$age, kidslt6 = raw$kidslt6,
  kidsge6 = raw$kidsge6)
csv <- "tests/testthat/fixtures/glm-gof-sas-mroz.csv"
utils::write.table(training, csv, sep = ",", row.names = FALSE, quote = FALSE)
file.copy(csv, file.path(root, "glm_gof_sas_mroz.csv"), overwrite = TRUE)

cases <- list()
rows <- function(x) lapply(seq_len(nrow(x)), function(i) unname(x[i, ]))
add <- function(label, x, y, n, family = "binomial", link = NULL,
                method = c("pearson", "deviance"), offset = rep(0, length(y)),
                g = 10) {
  if (is.null(link)) link <- if (family == "binomial") "logit" else "log"
  fam <- if (family == "binomial") stats::binomial(link) else
    stats::poisson(link)
  response <- if (family == "binomial") y/n else y
  fit <- stats::glm(response ~ x-1, family = fam, weights = n, offset = offset,
                    control = stats::glm.control(epsilon = 1e-12))
  result <- suppressWarnings(glm.gof(fit, method, g))
  tests <- if (inherits(result, "htest")) list(result) else result
  output <- lapply(unname(tests), function(z) list(
    statistic = unname(z$statistic),
    df = unname(z$parameter), p_value = z$p.value, method = z$test.method,
    groups = rows(as.matrix(z$groups))))
  cases[[length(cases) + 1L]] <<- list(label = label,
    model = list(x = rows(x), y = y, n = n, family = family, link = link,
                 offset = offset), args = list(method = as.list(method), g = g),
    result = output)
}
x <- cbind(1, seq(25, 75, 10))
y <- c(0, 0, 0, 7, 6, 17)
n <- c(14, 20, 19, 18, 12, 17)
add("sas-binomial", x, y, n)
add("logistic-extra-methods", x, y, n,
    method = c("hosmer-lemeshow", "osius-rojek", "stukel"), g = 5)
add("probit", x, y, n, link = "probit")
add("aliased-and-repeated", cbind(x, x[, 2]), y, n,
    method = c("deviance", "pearson", "deviance"))
add("offsets-with-same-design", cbind(1, rep(0:2, each = 2)),
    c(3, 5, 9, 10, 13, 16), rep(20, 6), offset = rep(c(0, 0.5), 3))
add("one-stukel-direction", cbind(1, 1:6), c(1, 3, 4, 2, 5, 7), rep(50, 6),
    method = "stukel")
x <- cbind(1, c(1, 0, 0, 1, 0, 0), c(0, 1, 0, 0, 1, 0),
            c(0, 0, 0, 1, 1, 1))
add("sas-poisson-offset", x, c(42, 37, 1, 101, 73, 14), rep(1, 6),
    family = "poisson", offset = log(c(500, 1200, 100, 400, 500, 300)))
add("poisson-frequencies", cbind(1, 0:5), c(2, 4, 3, 8, 7, 11),
    c(1, 2, 1, 3, 1, 2), family = "poisson")
paths <- c("R/glm_gof.R", "R/asymcheck.R")
json <- jsonlite::toJSON(list(source_md5 = as.list(tools::md5sum(paths)),
  training_md5 = unname(tools::md5sum(csv)), cases = cases),
  auto_unbox = TRUE, digits = NA, null = "null")
json <- sub("\n$", "", jsonlite::prettify(json, indent = 1))
writeLines(json, file.path(root, "glm_gof_r.json"))
