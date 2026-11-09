# Approximation criteria from the book. No inference backend is called.
.asymcheck_book_ci <- function(args) {
  family <- args$family
  choices <- switch(family,
    p1 = c("wilson", "wald", "agresti-coull", "clopper-pearson", "blaker"),
    p2 = c("newcombe", "wald", "miettinen-nurminen", "mee", "agresti-caffo",
           "brown-li", "chan-zhang"),
    rr = c("koopman", "katz"), or = c("wald", "gart", "cornfield", "bp"),
    rate2 = c("log", "score", "exact"), stop("unsupported CI family"))
  methods <- match.arg(args$method, choices, several.ok = TRUE)
  orientation <- resolve_ci_side(args$side, args$alternative)
  tol <- if ("tol" %in% names(args)) args$tol else 1e-8
  validate_binomial_ci_settings(args$conf.level, orientation$alternative,
    if (family == "or") 1e-8 else tol)
  if (family %in% c("p1", "p2")) {
    validate_ci_flag(args$correct, "correct")
    correctable <- if (family == "p1") c("wald", "wilson") else
      c("wald", "newcombe")
    if (args$correct && !any(methods %in% correctable)) {
      stop("correct applies to none of the requested methods")
    }
  }
  midp <- if (family %in% c("p1", "or")) args$midp else FALSE
  .asymcheck_flag(midp, "midp")
  if (family == "p1" && midp && !"clopper-pearson" %in% methods) {
    stop("midp applies only to Clopper-Pearson")
  }
  if (family == "p2" && (!is.numeric(args$nuisance.tol) ||
      length(args$nuisance.tol) != 1L || !is.finite(args$nuisance.tol) ||
      args$nuisance.tol <= 0 || args$nuisance.tol >= 1)) {
    stop("invalid nuisance.tol")
  }
  exact <- methods %in% c("clopper-pearson", "blaker", "chan-zhang",
                           "cornfield", "bp", "exact")
  if (family == "or" && any(exact) && (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0)) stop("invalid tol")
  if (family == "rate2") {
    .asymcheck_rate_inputs(args$x, args$exposure, sizes = 2L)
    total <- sum(args$x)
    if (total == 0) stop("at least one event is required for a rate-ratio CI")
    counts <- list(x = args$x[1L], n = total)
  } else if (family == "p1") {
    counts <- validate_binomial_counts(args$x, args$n)
  } else {
    if (family %in% c("rr", "or") && !is.matrix(args$x) &&
        !inherits(args$x, "table")) stop("x must be a 2 x 2 table")
    counts <- independent_binomial_counts(args$x, args$n)
  }
  p <- counts$x/counts$n
  value <- counts$n*p*(1-p)
  rows <- list()
  for (j in seq_along(methods)) {
    log.rate <- family == "rate2" && methods[j] == "log"
    values <- if (log.rate) args$x else value
    assessed <- !exact[j]
    criterion <- if (log.rate) "observed events" else "n * p.hat * (1-p.hat)"
    check <- .asymcheck_rows(methods[j], criterion, values,
      if (assessed) if (log.rate) 10 else 5 else NA_real_,
      if (assessed) if (log.rate) ">" else ">=" else NA_character_,
      if (!assessed) NA else if (log.rate) values > 10 else values >= 5)
    check$input <- if (family == "p1") seq_along(values) else 1L
    check$sample <- if (family == "p1") 1L else seq_along(values)
    rows[[j]] <- check
  }
  checks <- do.call(rbind, rows)
  if (family == "p1") checks <- checks[order(checks$input), ]
  rownames(checks) <- NULL
  chapter <- switch(family, p1 = "11", p2 = "12", rr = "15", or = "15",
                     rate2 = "22 (and 11 for conditional score)")
  .asymcheck_new(paste(switch(family, p1 = "one-proportion", p2 =
    "independent proportion-difference", rr = "risk-ratio", or = "odds-ratio",
    rate2 = "Poisson rate-ratio"), "interval"), methods, any(!exact), checks,
    list(counts = data.frame(x = counts$x, n = counts$n),
      sample.proportions = p, methods = data.frame(method = methods,
        applicable = !exact), conf.level = args$conf.level,
      side = orientation$side, alternative = orientation$alternative,
      correct = args$correct, midp = midp, tol = tol,
      nuisance.tol = args$nuisance.tol, exposure = args$exposure,
      original.events = if (family == "rate2") args$x else NULL,
      source = paste("Book chapter", chapter)),
    c("The book's approximation guidelines are not coverage guarantees.",
      "Criteria use original observations; pseudo-counts are not sample size.",
      if (family == "p1")
        "Ch. 11 binomial-normal guideline; Wald can under-cover even if met.",
      if (family %in% c("p2", "rr", "or"))
        "Check n_i * p.hat_i * (1-p.hat_i) >= 5 in each original sample.",
      if (family == "rate2")
        "Log-Wald: both y_i > 10; score: conditional-binomial CI guideline.",
      if (any(exact)) "Exact methods need no asymptotic approximation check.",
      if (midp) "Mid-p does not guarantee nominal coverage."))
}

.asymcheck_rate_inputs <- function(x, exposure, sizes = c(1L, 2L)) {
  if (!is.numeric(x) || !length(x) %in% sizes || any(!is.finite(x)) ||
      any(x < 0) || any(x != floor(x))) stop("invalid event counts")
  if (!is.numeric(exposure) || length(exposure) != length(x) ||
      any(!is.finite(exposure)) || any(exposure <= 0)) stop("invalid exposure")
}

.asymcheck_rate_test <- function(args) {
  .asymcheck_rate_inputs(args$x, args$exposure)
  settings <- validate_binomial_ci_settings(args$conf.level,
                                            args$alternative, 1e-8)
  .asymcheck_flag(args$correct, "correct")
  if (length(args$x) == 1L) {
    if (!is.numeric(args$r) || length(args$r) != 1L || !is.finite(args$r) ||
        args$r <= 0 || !is.finite(args$r*args$exposure)) {
      stop("invalid null rate")
    }
    value <- args$r*args$exposure
    threshold <- 10
    criterion <- "null expected events"
    p0 <- NULL
  } else {
    p0 <- args$exposure[1L]/sum(args$exposure)
    value <- sum(args$x)*p0*(1-p0)
    threshold <- 5
    criterion <- "n * pi0 * (1-pi0)"
  }
  .asymcheck_new("Poisson rate test", "normal", TRUE,
    .asymcheck_rows("normal", criterion, value, threshold, ">=",
      value >= threshold),
    list(x = args$x, exposure = args$exposure, r = args$r, pi0 = p0,
      correct = args$correct, alternative = settings$alternative,
      source = if (length(args$x) == 1L) "Book Ch. 21" else "Book Ch. 22"),
    c("Book null-distribution guideline; not an accuracy guarantee.",
      "This checks the test approximation, not its accompanying interval."))
}

.asymcheck_prop_test <- function(args, exact = FALSE) {
  x <- args$x
  n <- args[["n", exact = TRUE]]
  if (exact && length(x) == 2L) {
    n <- sum(x)
    x <- x[1L]
  } else if (!exact && is.matrix(x)) {
    if (ncol(x) != 2L || !is.numeric(x) || any(!is.finite(x)) ||
        any(x < 0) || any(x != floor(x))) stop("invalid success/failure table")
    n <- rowSums(x)
    x <- x[, 1L]
  } else if (!exact && is.table(x) && length(dim(x)) == 1L && length(x) == 2L) {
    n <- sum(x)
    x <- x[1L]
  }
  if (length(x) != length(n) || !length(x) %in% c(1L, 2L) ||
      (exact && length(x) != 1L)) stop("supply one or two matching samples")
  counts <- validate_binomial_counts(x, n)
  validate_binomial_ci_settings(args$conf.level, args$alternative, 1e-8)
  p0 <- args$p
  if (length(x) == 2L && !is.null(p0)) {
    stop("two-sample pre-checks support the equality null (p = NULL)")
  }
  if (is.null(p0)) p0 <- if (length(x) == 1L) 0.5 else sum(x)/sum(n)
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) || p0 < 0 ||
      p0 > 1 || (!exact && length(x) == 1L && p0 %in% c(0, 1))) {
    stop("invalid null probability")
  }
  if (!exact) .asymcheck_flag(args$correct, "correct")
  value <- if (length(x) == 1L) n*p0*(1-p0) else
    c(rbind(n*p0, n*(1-p0)))
  criterion <- if (length(x) == 1L) "n * pi0 * (1-pi0)" else
    rep(c("null expected successes", "null expected failures"), 2L)
  method <- if (exact) "exact-binomial" else "score"
  .asymcheck_new("proportion test", method, !exact,
    .asymcheck_rows(method, criterion, value, if (exact) NA_real_ else 5,
      if (exact) NA_character_ else ">=", if (exact) NA else value >= 5),
    list(counts = counts, null.probability = p0,
      source = if (length(x) == 1L) "Book Ch. 11" else "Book Ch. 12"),
    if (exact) "Exact binomial inference needs no asymptotic check." else
      "Book null-distribution guideline; not an accuracy guarantee.")
}
