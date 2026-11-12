#' Print an ibist Confidence Interval Object
#'
#' @param x An object of class \code{"ibist_ci"}.
#' @param digits Number of significant digits to print.
#' @param ... Reserved for future extensions.
#'
#' @export
print.ibist_ci <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  interval <- if (nrow(x$conf.int) == 1L) "interval" else "intervals"
  cat("\n\tConfidence ", interval, " for ", x$parameter, "\n", sep = "")
  if (!is.null(x$data.name)) cat("\ndata: ", x$data.name, "\n", sep = "")
  if (length(x$method) == 1L) cat("method: ", x$method, "\n", sep = "")
  if (!is.null(x$side) && x$side != "two.sided") {
    cat("bound: ", x$side,
        "\n", sep = "")
  }
  cat("\nsample estimates:\n")
  print(x$estimate, digits = digits)
  cat(
    "\n", 100 * x$conf.level, "% confidence ", interval,
    ":\n", sep = ""
  )
  intervals <- x$conf.int
  attr(intervals, "conf.level") <- NULL
  if (nrow(intervals) == 1L) {
    print(as.numeric(intervals), digits = digits)
  } else {
    print(intervals, digits = digits)
  }
  invisible(x)
}

new_ibist_ci <- function(method, estimate, intervals, conf.level, parameter,
                         data.name = NULL, null.value = NULL,
                         side = "two.sided") {
  conf.int <- cbind(
    lower = vapply(intervals, `[`, numeric(1), 1L),
    upper = vapply(intervals, `[`, numeric(1), 2L)
  )
  rownames(conf.int) <- method
  attr(conf.int, "conf.level") <- conf.level

  structure(
    list(
      conf.int = conf.int,
      estimate = estimate,
      conf.level = conf.level,
      side = side,
      method = method,
      parameter = parameter,
      data.name = data.name,
      null.value = null.value
    ),
    class = "ibist_ci"
  )
}

# NULL distinguishes an omitted side from an explicitly supplied value.
resolve_ci_side <- function(side = NULL, alternative = NULL) {
  choices <- c("two.sided", "lower", "upper")
  if (!is.null(side) && (length(side) != 1L || !is.character(side) ||
                        is.na(side) || !side %in% choices)) {
    stop("'side' must be two.sided, lower, or upper.")
  }
  if (!is.null(alternative)) {
    alternatives <- c("two.sided", "greater", "less")
    if (length(alternative) != 1L || !is.character(alternative) ||
        is.na(alternative) || !alternative %in% alternatives) {
      stop("'alternative' must be two.sided, less, or greater.")
    }
    legacy <- choices[match(alternative, alternatives)]
    if (!is.null(side) && side != legacy) {
      stop("'side' and 'alternative' specify conflicting intervals.")
    }
    side <- legacy
  }
  if (is.null(side)) side <- "two.sided"
  list(side = side, alternative = c("two.sided", "greater", "less")[
    match(side, choices)])
}

ci_tail <- function(conf.level, side) {
  (1 - conf.level) / if (side == "two.sided") 2 else 1
}

ci_bounds <- function(limits, side, domain = c(0, Inf)) {
  if (side == "lower") limits[2L] <- domain[2L]
  if (side == "upper") limits[1L] <- domain[1L]
  limits
}
