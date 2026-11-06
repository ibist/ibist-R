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
  if (!is.null(x$alternative) && x$alternative != "two.sided") {
    cat("bound: ", if (x$alternative == "greater") "lower" else "upper",
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
                         data.name = NULL, null.value = NULL) {
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
      method = method,
      parameter = parameter,
      data.name = data.name,
      null.value = null.value
    ),
    class = "ibist_ci"
  )
}
