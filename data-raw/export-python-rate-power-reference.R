# Small representative export; independent formulas are checked in tests.
source("R/power_rates.R")
cases <- list()
add <- function(kind, args) {
  fun <- if (kind == "one") power.r1s.test else power.r2s.test
  result <- unclass(do.call(fun, args))
  names(result) <- gsub("\\.", "_", names(result))
  names(args) <- gsub("\\.", "_", names(args))
  cases[[length(cases) + 1L]] <<- list(kind = kind, args = args,
                                     result = result)
}
for (kind in c("one", "two")) {
  args <- if (kind == "one") list(n = 100, r0 = 0.1, r1 = 0.2) else
    list(n = 100, r1 = 0.1, r2 = 0.2)
  add(kind, args)
  for (name in c("n", if (kind == "one") c("r0", "r1") else
                c("r1", "r2"), "sig.level")) {
    inverse <- args
    inverse[name] <- list(NULL)
    inverse$power <- 0.8
    add(kind, inverse)
  }
}
add("one", list(n = 31.5, r0 = 0.2, r1 = 0.1, exposure = 3,
                 alternative = "less", correct = TRUE))
add("one", list(n = 40, r0 = 0.1, r1 = 0.2, exposure = 2,
                 alternative = "greater"))
add("one", list(n = 15, r0 = 1.5, r1 = 2.5, strict = TRUE,
                 correct = TRUE))
add("one", list(n = 10, r0 = 0.1, r1 = 0.1, strict = TRUE))
add("one", list(n = 50, r0 = 0.2, power = 0.8, alternative = "less"))
add("two", list(n = 50, r1 = 0.1, r2 = 0.2, group.rate = 2,
                 exposure = c(2, 3), strict = TRUE, correct = TRUE))
add("two", list(r1 = 0.1, r2 = 0.2, power = 0.8, group.rate = 0.5,
                 exposure = c(1, 2), alternative = "one.sided"))
add("two", list(n = 15, r1 = 1.5, r2 = 2.5, correct = TRUE))
add("two", list(n = 10, r1 = 0.1, r2 = 0.1, strict = TRUE))
path <- "R/power_rates.R"
output <- list(source = path, source_md5 = unname(tools::md5sum(path)),
               cases = cases)
jsonlite::write_json(output,
  "../ibist-py/tests/reference/rate_power_r.json", auto_unbox = TRUE,
  null = "null", digits = NA, pretty = TRUE)
