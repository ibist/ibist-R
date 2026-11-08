#' Tests of Nominal and Ordinal Association
#'
#' Tests zero association using one or more named measures and a common
#' inference method. Existing measure and test functions remain available.
#'
#' @param x A two-dimensional matrix or table of non-negative integer counts.
#'   For ordinal inference, rows and columns must already be ordered.
#' @param measure One or more measure names. The default is \code{"pearson"}.
#'   Nominal measures are \code{"phi"}, \code{"contingency"}, and
#'   \code{"cramer-v"}. Ordinal measures are \code{"gamma"}, \code{"tau-b"},
#'   \code{"tau-c"}, \code{"somers-c|r"}, \code{"somers-r|c"},
#'   \code{"pearson"}, and \code{"spearman"}. Preserve order and repetitions.
#' @param method Inference procedure: \code{"multinomial"} (default),
#'   \code{"pearson"}, or \code{"cmh"}. One method applies to the whole call.
#' @param alternative Alternative hypothesis: \code{"two.sided"} (default),
#'   \code{"greater"}, or \code{"less"}. Nominal Pearson independence tests
#'   require \code{"two.sided"}; directional tests apply to ordinal measures.
#' @param conf.level Confidence level, default \code{0.95}. Used for
#'   multinomial intervals only; validated for every method.
#' @param row.scores,col.scores Optional Pearson category scores. Multinomial
#'   inference follows \code{ordinal.assoc()}. For CMH, supply both or neither;
#'   scores apply to Pearson only in a mixed Pearson/Spearman request.
#'   Scores are errors for nominal inference or a CMH Spearman-only request.
#' @param correct Logical, default \code{FALSE}. Request capped Yates
#'   correction for nominal Pearson tests of 2 x 2 tables only. Nominal
#'   coefficient estimates always use the uncorrected Pearson statistic.
#'
#' @details
#' \code{method = "pearson"} supports only the three nominal measures.
#' They share the independence test from \code{indep.test()}, while their
#' descriptive estimates come from \code{nominal.assoc()}. No nominal
#' coefficient confidence interval is currently supplied. Empty margins are
#' rejected. The null coefficient is zero, equivalently independence.
#'
#' \code{method = "multinomial"} supports all seven ordinal measures and
#' delegates their estimates, standard errors, Wald intervals, and tests
#' to \code{ordinal.assoc()}. The five concordance-based measures share one
#' test statistic; Pearson and Spearman use their own null variances.
#' Two-sided intervals are untruncated; one-sided finite bounds are restricted
#' to [-1, 1]. Intervals use estimation variances and need not invert tests.
#'
#' \code{method = "cmh"} supports Pearson and Spearman only. It uses
#' \code{cmh.ll.test()}, with signed statistic \eqn{Z=r\sqrt{n-1}} and
#' \code{chisq = Z^2}. This is different from multinomial inference for the
#' same coefficient. No CMH confidence interval is supplied. Ordinal methods
#' drop empty margins; see the underlying functions for score definitions.
#'
#' Incompatible measure/method combinations are errors, including mixed
#' nominal and ordinal requests. Multiple results use the same observations;
#' their tests are not independent and p-values are not multiplicity-adjusted.
#'
#' @return One measure returns an \code{htest} with an estimate, zero null,
#'   statistic, p-value, alternative, and method description. Additional fields
#'   \code{measure} and \code{test.method} identify the selected coefficient
#'   and inference procedure. Multinomial results also contain
#'   \code{std.error} and \code{conf.int}; nominal Pearson results retain
#'   independence diagnostics and degrees of freedom; CMH results include
#'   \code{chisq}. Multiple measures return an ordered named list of
#'   \code{htest} objects of class \code{ibist_assoc_tests}.
#'   \code{summary()} and \code{as.data.frame()} return a full-precision
#'   combined table. Printing identifies shared tests and rounds display only.
#'
#' @references
#' Brown, M. B., and Benedetti, J. K. (1977). Sampling behavior of tests for
#' correlation in two-way contingency tables. \emph{Journal of the American
#' Statistical Association}, 72(358), 309--315.
#' \doi{10.1080/01621459.1977.10480995}
#'
#' Mantel, N. (1963). Chi-square tests with one degree of freedom: extensions
#' of the Mantel-Haenszel procedure. \emph{Journal of the American Statistical
#' Association}, 58(303), 690--700.
#' \doi{10.1080/01621459.1963.10500879}
#'
#' SAS Institute Inc. (2016). \emph{SAS/STAT 14.2 User's Guide}.
#' The FREQ Procedure, Chi-Square Tests and Statistics; Measures of Association.
#' \url{https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf}
#'
#' @seealso \code{\link{nominal.assoc}}, \code{\link{ordinal.assoc}},
#'   \code{\link{indep.test}}, \code{\link{cmh.ll.test}}
#'
#' @examples
#' tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), 3, byrow = TRUE)
#' assoc.test(tab)
#' assoc.test(tab, measure = c("gamma", "tau-b", "somers-c|r"))
#' assoc.test(tab, measure = c("phi", "cramer-v"), method = "pearson")
#' assoc.test(tab, measure = c("pearson", "spearman"), method = "cmh")
#'
#' @export
assoc.test <- function(x, measure = "pearson", method = "multinomial",
                      alternative = "two.sided", conf.level = 0.95,
                      row.scores = NULL, col.scores = NULL, correct = FALSE) {
  nominal <- c("phi", "contingency", "cramer-v")
  ordinal <- c("gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c",
               "pearson", "spearman")
  if (!is.character(measure) || !length(measure) || anyNA(measure)) {
    stop("'measure' must contain one or more valid measure names.")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !is.character(alternative) || length(alternative) != 1L ||
      is.na(alternative)) {
    stop("'method' and 'alternative' must each be one valid name.")
  }
  measure <- match.arg(measure, c(nominal, ordinal), several.ok = TRUE)
  method <- match.arg(method, c("multinomial", "pearson", "cmh"))
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  valid <- switch(method, pearson = nominal, multinomial = ordinal,
                  cmh = c("pearson", "spearman"))
  if (!all(measure %in% valid)) {
    stop("Selected measures are incompatible with method = '", method, "'.")
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }
  if (!is.logical(correct) || length(correct) != 1L || is.na(correct)) {
    stop("'correct' must be TRUE or FALSE.")
  }
  if (correct && method != "pearson") {
    stop("'correct = TRUE' is only available for nominal Pearson tests.")
  }
  scores <- !is.null(row.scores) || !is.null(col.scores)
  if (method == "pearson" && (scores || alternative != "two.sided")) {
    stop("Nominal tests require 'two.sided' and do not accept category scores.")
  }
  if (method == "cmh" &&
      (xor(is.null(row.scores), is.null(col.scores)) ||
       (scores && !"pearson" %in% measure))) {
    stop("CMH scores require both axes and a Pearson measure request.")
  }
  data.name <- paste(deparse(substitute(x)), collapse = "\n")
  if (method == "multinomial") {
    fit <- ordinal.assoc(x, measure, conf.level, alternative,
                         row.scores, col.scores)
    results <- lapply(seq_along(measure), function(i) {
      ci <- unname(fit$conf.int[i, ])
      attr(ci, "conf.level") <- conf.level
      structure(list(
        statistic = c(Z = unname(fit$statistic[i])),
        p.value = unname(fit$p.value[i]),
        estimate = stats::setNames(unname(fit$estimate[i]), measure[i]),
        null.value = stats::setNames(0, measure[i]),
        alternative = alternative,
        method = paste("Multinomial test of zero", measure[i], "association"),
        data.name = data.name, conf.int = ci,
        std.error = unname(fit$std.error[i])), class = "htest")
    })
  } else if (method == "pearson") {
    fit <- indep.test(x, method = "pearson", correct = correct)
    estimates <- nominal.assoc(x, measure = unique(measure))$estimate
    results <- lapply(measure, function(name) {
      out <- fit
      out$estimate <- stats::setNames(unname(estimates[name]), name)
      out$null.value <- stats::setNames(0, name)
      out$alternative <- alternative
      out$method <- paste("Pearson independence test:", name)
      if (correct) out$method <- paste0(out$method, "; Yates correction")
      out
    })
  } else {
    results <- lapply(measure, function(name) {
      fit <- cmh.ll.test(
        x, name, alternative,
        row.scores = if (name == "pearson") row.scores else NULL,
        col.scores = if (name == "pearson") col.scores else NULL)
      fit$data.name <- data.name
      names(fit$estimate) <- names(fit$null.value) <- name
      fit$method <- paste("CMH linear-by-linear test:", name)
      fit
    })
  }
  for (i in seq_along(results)) {
    results[[i]]$measure <- measure[i]
    results[[i]]$test.method <- method
    results[[i]]$data.name <- data.name
  }
  if (length(results) == 1L) return(results[[1L]])
  names(results) <- measure
  group <- if (method == "pearson") measure else if (method == "multinomial") {
    measure[measure %in% ordinal[1:5]]
  } else character()
  shared <- if (length(unique(group)) > 1L) list(group) else list()
  structure(results, class = c("ibist_assoc_tests", "list"),
            shared.tests = shared)
}

#' @rdname assoc.test
#' @param object An ordered collection returned by \code{assoc.test()}.
#' @param ... Reserved for future extensions.
#' @export
summary.ibist_assoc_tests <- function(object, ...) {
  scalar <- function(x, name) {
    if (is.null(x[[name]])) NA_real_ else unname(x[[name]][1L])
  }
  rows <- lapply(object, function(x) {
    ci <- if (is.null(x$conf.int)) c(NA_real_, NA_real_) else x$conf.int
    data.frame(measure = x$measure, estimate = scalar(x, "estimate"),
               std.error = scalar(x, "std.error"), lower = ci[1L],
               upper = ci[2L], statistic = scalar(x, "statistic"),
               df = scalar(x, "parameter"), p.value = x$p.value,
               row.names = NULL, check.names = FALSE)
  })
  out <- do.call(rbind, unname(rows))
  rownames(out) <- NULL
  out
}

#' @rdname assoc.test
#' @param row.names Optional row names for the combined data frame.
#' @param optional Passed through for compatibility with \code{as.data.frame}.
#' @usage \method{as.data.frame}{ibist_assoc_tests}(
#'   x, row.names = NULL, optional = FALSE, ...
#' )
#' @export
as.data.frame.ibist_assoc_tests <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  out <- summary(x)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' @rdname assoc.test
#' @param digits Significant digits used for printing.
#' @export
print.ibist_assoc_tests <- function(x, digits = 7L, ...) {
  cat("\nAssociation tests (", x[[1L]]$test.method, ")\n\n", sep = "")
  print(summary(x), digits = digits, row.names = FALSE)
  for (group in attr(x, "shared.tests")) {
    note <- paste("Shared test:", paste(group, collapse = ", "))
    cat("\n", paste(strwrap(note, width = getOption("width")),
                    collapse = "\n"), "\n", sep = "")
  }
  cat("\nTests use the same data; p-values are unadjusted.\n")
  invisible(x)
}
