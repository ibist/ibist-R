#' Diagnostics for Asymptotic Tests and Confidence Intervals
#'
#' Inspects inputs before inference or retained inputs after inference.
#' Criteria are method-specific heuristics, not guarantees of accuracy.
#'
#' @param object A supported inference function, a supported result, or raw
#'   counts with an explicit \code{method}.
#' @param ... Arguments of the inference function for a pre-check. Raw
#'   one-rate counts accept the arguments of \code{rate.1s.ci()}.
#' @param method Raw-data procedure: \code{"pearson"}, \code{"lr"},
#'   \code{"fisher"}, or \code{"mantel-haenszel"}. Rate CI methods
#'   \code{"exact"}, \code{"score"}, \code{"wh"}, \code{"wald"}, and
#'   \code{"log"} select one-rate intervals. For a function pre-check,
#'   this is that function's method argument. Omit for result objects.
#'
#' @details
#' Function pre-checks support \code{indep.test()}, \code{rate.1s.ci()},
#' and \code{stats::mantelhaen.test()} for 2 x 2 x K tables. They do not
#' run the test or interval. Results from \code{indep.test()} and newly
#' computed \code{rate.1s.ci()} objects support post-checks. Arbitrary
#' htest objects without the required inputs are not reconstructed by
#' evaluating their printed data names.
#'
#' Pearson and G-squared checks report minimum expected count >= 1 and
#' fraction of expected counts >= 5 of at least 0.8. These are screening
#' heuristics, not universal adequacy conditions. The G-squared report
#' additionally notes that these criteria were developed for Pearson's
#' approximation. Mantel-Haenszel retains the legacy sum of conditional
#' hypergeometric variances >= 5 rule. Empty strata contribute zero.
#'
#' Approximate rate intervals report observed events and the plug-in
#' Poisson skewness 1/sqrt(x), undefined at zero. No universal cutoff is
#' asserted: their condition.met values are NA. Exact procedures are
#' marked not applicable, including exact fallbacks for zero-event log
#' and uncorrected Wald intervals. Dependence, sampling assumptions, and
#' common-effect homogeneity are not assessed. No p-value or interval
#' decision enters a diagnostic, and no inference method is changed.
#'
#' @return An \code{"ibist_asymcheck"} object with \code{procedure},
#'   \code{method}, \code{applicable}, \code{checks} (a data frame),
#'   \code{details}, \code{notes}, and \code{condition.met}. The last is
#'   NA if no numerical criterion was assessed. It summarizes only the
#'   stated criteria, not all assumptions needed for valid inference.
#' @references R Core Team. Fisher's Exact Test for Count Data,
#'   section Details (Cochran conditions).
#'   R stats documentation, \url{https://stat.ethz.ch/R-manual/}.
#'
#'   NIST/SEMATECH e-Handbook, section 1.3.6.6.19, Poisson Distribution
#'   (skewness). \url{https://www.itl.nist.gov/div898/handbook/}.
#' @seealso \code{\link{indep.test}}, \code{\link{rate.1s.ci}},
#'   \code{\link{mantelhaen.asymcheck}}
#' @examples
#' tab <- matrix(c(12, 8, 5, 15), 2)
#' asymcheck(indep.test, x = tab)
#' fit <- indep.test(tab, asymcheck = TRUE)
#' fit$diagnostics
#' asymcheck(fit)
#' asymcheck(rate.1s.ci, x = 12, exposure = 100, method = "wald")
#' @export
asymcheck <- function(object, ..., method = NULL) {
  args <- list(...)
  if (is.function(object)) {
    if (!is.null(method)) args$method <- method
    if (identical(object, indep.test)) {
      args <- .asymcheck_bind(object, args)
      return(.asymcheck_indep(args$x, args$method, args$correct))
    }
    if (identical(object, rate.1s.ci)) {
      args <- .asymcheck_bind(object, args)
      return(.asymcheck_rate(args))
    }
    if (identical(object, stats::mantelhaen.test)) {
      args <- .asymcheck_bind(object, args)
      dims <- dim(args$x)
      if (!is.null(args$y) || !is.null(args$z) || !is.array(args$x) ||
          length(dims) != 3L || any(dims < 2L)) {
        stop("Mantel-Haenszel pre-check requires a 2 x 2 x K count array")
      }
      return(.asymcheck_mh(args$x, exact = args$exact))
    }
    stop("unsupported inference function")
  }
  if (inherits(object, "htest") || inherits(object, "ibist_ci")) {
    if (length(args) || !is.null(method)) {
      stop("result post-checks do not accept method or additional inputs")
    }
    if (inherits(object, "htest") && !is.null(object$observed) &&
        !is.null(object$expected)) {
      label <- object$method
      selected <- if (grepl("Fisher", label, fixed = TRUE)) "fisher" else
        if (grepl("G-squared", label, fixed = TRUE)) "lr" else
          if (grepl("Pearson", label, fixed = TRUE)) "pearson" else NULL
      if (!is.null(selected)) {
        return(.asymcheck_indep(object$observed, selected,
                               grepl("Yates", label, fixed = TRUE)))
      }
    }
    if (inherits(object, "ibist_ci") &&
        identical(object$parameter, "Poisson rate") &&
        !is.null(object$inputs)) {
      return(.asymcheck_rate(list(x = object$inputs$x,
        exposure = object$inputs$exposure, method = object$method,
        conf.level = object$conf.level, correct = object$correct,
        side = object$side, alternative = NULL)))
    }
    stop("unsupported result or missing retained inputs; use a pre-check")
  }
  if (is.null(method)) stop("raw inputs require an explicit method")
  if (identical(method, "mantel-haenszel")) {
    if (length(args)) stop("raw Mantel-Haenszel checks accept only counts")
    return(.asymcheck_mh(object))
  }
  if (length(method) == 1L && method %in% c("pearson", "lr", "fisher")) {
    args$x <- object
    args$method <- method
    return(do.call(asymcheck, c(list(object = indep.test), args)))
  }
  args$x <- object
  args$method <- method
  do.call(asymcheck, c(list(object = rate.1s.ci), args))
}

.asymcheck_bind <- function(fun, args) {
  call <- as.call(c(list(as.name("precheck")), args))
  matched <- as.list(match.call(fun, call, expand.dots = FALSE))[-1L]
  if (is.null(matched$x)) stop("a pre-check requires 'x'")
  defaults <- formals(fun)
  for (name in setdiff(names(defaults), c(names(matched), "...", "x"))) {
    matched[[name]] <- eval(defaults[[name]], envir = environment(fun))
  }
  # match.arg uses the first default; an explicit rate vector is retained.
  if (is.null(args$method) && !is.null(matched$method)) {
    matched$method <- matched$method[1L]
  }
  if (length(matched$...)) {
    stop("additional inference options are not supported by this pre-check")
  }
  if (!is.null(matched$asymcheck)) {
    .asymcheck_flag(matched$asymcheck, "asymcheck")
  }
  matched
}

.asymcheck_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(paste0("'", name, "' must be TRUE or FALSE"))
  }
}

.asymcheck_table <- function(x, margins = TRUE) {
  tab <- as.matrix(x)
  if (!is.numeric(tab) || length(dim(tab)) != 2L || any(dim(tab) < 2L) ||
      any(!is.finite(tab)) || any(tab < 0) || any(tab != floor(tab)) ||
      sum(tab) > 2^53 - 1) stop("invalid non-negative integer count table")
  if (margins && (any(rowSums(tab) == 0) || any(colSums(tab) == 0))) {
    stop("table must have positive margins")
  }
  tab
}

.asymcheck_new <- function(procedure, method, applicable, checks,
                           details = list(), notes = character()) {
  assessed <- checks$condition.met[!is.na(checks$condition.met)]
  structure(list(procedure = procedure, method = method,
    applicable = applicable, checks = checks, details = details,
    notes = c(notes, "Criteria do not verify all sampling assumptions."),
    condition.met = if (length(assessed)) all(assessed) else NA),
    class = "ibist_asymcheck")
}

.asymcheck_rows <- function(method, criterion, value, threshold = NA_real_,
                            comparison = NA_character_, met = NA) {
  data.frame(method = method, criterion = criterion, value = value,
    threshold = threshold, comparison = comparison, condition.met = met)
}

.asymcheck_indep <- function(x, method, correct = FALSE) {
  method <- match.arg(method, c("pearson", "lr", "fisher"))
  .asymcheck_flag(correct, "correct")
  tab <- .asymcheck_table(x)
  if (correct && (method != "pearson" || any(dim(tab) != 2L))) {
    stop("correction requires a 2 x 2 Pearson test")
  }
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  dimnames(expected) <- dimnames(tab)
  values <- c(min(expected), mean(expected >= 5))
  applicable <- method != "fisher"
  checks <- .asymcheck_rows(method,
    c("minimum expected count", "fraction expected >= 5"), values,
    if (applicable) c(1, 0.8) else NA_real_,
    if (applicable) ">=" else NA_character_,
    if (applicable) values >= c(1, 0.8) else NA)
  notes <- if (!applicable) {
    "Conditional Fisher inference needs no asymptotic check."
  } else "Expected-count rules are approximation screening heuristics."
  if (method == "lr") notes <- c(notes,
    "Pearson screening rules do not establish G-squared accuracy.")
  .asymcheck_new("independence", method, applicable, checks,
    list(expected = expected, correct = correct), notes)
}

.asymcheck_mh <- function(x, exact = FALSE) {
  .asymcheck_flag(exact, "exact")
  dims <- dim(x)
  tables <- if (is.array(x) && length(dims) == 3L &&
                all(dims[1:2] == 2L)) {
    lapply(seq_len(dims[3L]), function(i) x[, , i])
  } else if (is.list(x)) x else stop("use a list or 2 x 2 x K array")
  if (!length(tables)) stop("at least one stratum is required")
  tables <- lapply(tables, function(tab) {
    tab <- .asymcheck_table(tab, margins = FALSE)
    if (any(dim(tab) != 2L)) stop("strata must be 2 x 2 tables")
    tab
  })
  variance <- vapply(tables, mantelhaen_stratum_variance, numeric(1))
  names(variance) <- paste("Stratum", seq_along(variance))
  value <- sum(variance)
  checks <- .asymcheck_rows("mantel-haenszel", "sum of variances", value,
    if (exact) NA_real_ else 5, if (exact) NA_character_ else ">=",
    if (exact) NA else value >= 5)
  .asymcheck_new("stratified independence", "mantel-haenszel", !exact,
    checks, list(variance = variance, variance.sum = value),
    if (exact) "Exact conditional inference needs no asymptotic check." else
      "Legacy variance-sum heuristic; common-effect homogeneity is unchecked.")
}

.asymcheck_rate <- function(args) {
  x <- args$x
  exposure <- args$exposure
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 ||
      x != floor(x) || x > 2^53 - 1) stop("invalid event count")
  for (name in c("exposure", "conf.level")) {
    value <- args[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) stop(paste(name, "must be positive and finite"))
  }
  if (args$conf.level >= 1) stop("conf.level must be in (0, 1)")
  .asymcheck_flag(args$correct, "correct")
  orientation <- resolve_ci_side(args$side, args$alternative)
  methods <- args$method
  if (!is.character(methods) || !length(methods) || anyNA(methods) ||
      any(!methods %in% c("exact", "score", "wh", "wald", "log"))) {
    stop("unsupported one-rate interval method")
  }
  fallback <- x == 0 & (methods == "log" |
    (methods == "wald" & !args$correct))
  approximate <- methods != "exact" & !fallback
  checks <- do.call(rbind, lapply(seq_along(methods), function(i) {
    .asymcheck_rows(methods[i], c("observed events", "plug-in skewness"),
      c(x, if (x > 0 && approximate[i]) 1 / sqrt(x) else NA_real_))
  }))
  .asymcheck_new("one-sample Poisson rate interval", methods,
    any(approximate), checks,
    list(x = x, exposure = exposure, conf.level = args$conf.level,
      side = orientation$side, correct = args$correct,
      approximate = approximate, exact.fallback = fallback),
    c("No universal event-count cutoff is asserted for interval accuracy.",
      if (x == 0) "Zero events: plug-in skewness is undefined.",
      if (any(fallback)) "Zero-event log/uncorrected Wald uses exact fallback.",
      if (!any(approximate)) "No asymptotic interval is being assessed."))
}

#' @rdname asymcheck
#' @param x A diagnostic result to print.
#' @param digits Number of significant digits for display only.
#' @export
print.ibist_asymcheck <- function(x, digits = 4L, ...) {
  cat("\nAsymptotic diagnostics:", x$procedure, "\n")
  print(x$checks, digits = digits, row.names = FALSE)
  cat(paste(x$notes, collapse = "\n"), "\n")
  invisible(x)
}
