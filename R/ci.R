#' Print an ibist Confidence Interval Object
#'
#' @param x An object of class \code{"ibist_ci"}.
#' @param digits Number of significant digits to print.
#' @param ... Reserved for future extensions.
#'
#' @export
print.ibist_ci <- function(x, digits = getOption("digits"), ...) {
  interval <- if (nrow(x$conf.int) == 1L) "interval" else "intervals"
  cat(
    "\n", 100 * x$conf.level, "% confidence ", interval,
    " for ", x$parameter, "\n\n", sep = ""
  )
  if (length(x$estimate) == 1L) {
    cat("Estimate: ", format(x$estimate, digits = digits), "\n\n", sep = "")
  } else {
    cat("Estimates:\n")
    print(x$estimate, digits = digits)
    cat("\n")
  }
  intervals <- x$conf.int
  attr(intervals, "conf.level") <- NULL
  print(intervals, digits = digits)
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
