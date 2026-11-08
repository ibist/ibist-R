#' Goodness-of-Fit Tests for Generalized Linear Models
#'
#' Tests binomial and Poisson GLMs using methods appropriate to their family
#' and link. The fitted model is the first argument; it is never refitted.
#'
#' @param fit A converged, unpenalized \code{glm} object. Supported families
#'   are binomial and Poisson with fixed dispersion one. Prior weights must
#'   represent integer trial counts or integer frequencies, not survey weights.
#'   Retain y=TRUE and either model=TRUE or x=TRUE when fitting.
#' @param method One or more of \code{"pearson"} (default),
#'   \code{"deviance"}, \code{"hosmer-lemeshow"}, \code{"osius-rojek"}, or
#'   \code{"stukel"}. The last three require a binomial logit model.
#'   Order and repetitions are retained; p-values are unadjusted.
#' @param g Requested number of Hosmer-Lemeshow risk groups, integer >= 3.
#'   Tied fitted probabilities are kept together; actual groups can be fewer.
#'   Ties are identified at 14 significant digits for numerical stability.
#' @param asymcheck Logical; attach approximation diagnostics under
#'   \code{diagnostics}. Default \code{FALSE}. Standalone checks are also
#'   available through \code{asymcheck(glm.gof, fit = fit, ...)} and
#'   \code{asymcheck(result)} without computing a test in the pre-check.
#'
#' @details
#' Pearson and deviance aggregate identical model-matrix rows AND offsets.
#' Binomial responses may be binary, weighted proportions, or success/failure
#' matrices. Groups have observed events O, expected events E, trial totals n
#' (or Poisson frequency totals), and model rank r; df is J-r. No automatic
#' pooling is done. Small expected counts can invalidate the chi-square
#' reference. Poisson aggregation assesses grouped means, not all aspects of
#' individual-count dispersion. Offsets and exposure offsets are retained.
#'
#' Hosmer-Lemeshow uses SAS's target-frequency grouping rule, always retaining
#' probability ties. Its statistic uses both expected events and nonevents,
#' with actual group count minus two df. An intercept must be in the model.
#' Grouping may differ from software that splits ties or uses sample quantiles.
#'
#' Osius-Rojek follows the SAS normal standardization: center the profile
#' Pearson statistic at J, not J-r. For profile trial totals n and probabilities
#' p, its variance is 2*sum(1-1/n) plus the squared norm of the residuals from
#' projecting (1-2p)/sqrt(n*p*(1-p)) onto the weighted model space. The reported
#' statistic is the squared standardized statistic, with one df.
#'
#' Stukel uses a score test of eta^2*I(eta>=0) and eta^2*I(eta<0), where eta
#' includes offsets. Weighted projection removes nuisance model terms. Degrees
#' of freedom equal the estimable increment, which can be one instead of two.
#' No augmented model is fitted. Singular or zero-variance tests are errors.
#' Stable projections use column-scaled SVD with relative tolerance 1e-10.
#' Projection noise below 1e-10*max(original column norm, 1) is discarded.
#' Standardizing variances <= 1e-12 are treated as numerically degenerate.
#'
#' All five methods are asymptotic. Diagnostics report expected counts for
#' grouped tests, positive standardizing variance for Osius-Rojek, and
#' estimable added terms for Stukel. They do not certify finite-sample accuracy.
#' Quasi, estimated-dispersion, constrained, and penalized fits are outside
#' this interface. Existing logistic.gof() retains its separate pooling rules.
#'
#' @return One method returns an htest with statistic, parameter (df), p.value,
#'   test.method, family, link, and groups. Multiple methods return an ordered
#'   \code{ibist_glm_gof_tests} list of these objects, with a combined summary.
#'   Inputs needed for diagnostic post-checks are retained independently of
#'   the asymcheck option. No model fit is stored in the result.
#' @references SAS Institute. SAS/STAT User's Guide, LOGISTIC Procedure,
#'   Goodness-of-Fit Tests (pp. 6050-6053, chapter 79).
#'   \url{https://support.sas.com/documentation/}.
#'   Hosmer and Lemeshow (1980), Communications in Statistics A9, 1043-1069.
#'   Osius and Rojek (1992), JASA 87, 1145-1152.
#'   Stukel (1988), JASA 83, 426-431.
#' @seealso \code{\link{logistic.gof}}, \code{\link{asymcheck}}
#' @examples
#' d <- data.frame(age = seq(25, 75, 10), y = c(0, 0, 0, 7, 6, 17),
#'                 n = c(14, 20, 19, 18, 12, 17))
#' fit <- glm(cbind(y, n-y) ~ age, data = d, family = binomial())
#' glm.gof(fit, method = c("pearson", "deviance"), asymcheck = TRUE)
#' @export
glm.gof <- function(fit, method = "pearson", g = 10, asymcheck = FALSE) {
  .asymcheck_flag(asymcheck, "asymcheck")
  data <- .glm_gof_data(fit)
  method <- .glm_gof_methods(method, g, data)
  results <- lapply(method, function(name) {
    parts <- .glm_gof_parts(data, name, g)
    check <- .glm_gof_check(data, name, parts)
    if (identical(check$condition.met, FALSE)) {
      warning(name, ": approximation criterion not met; inspect diagnostics",
              call. = FALSE)
    }
    statistic <- .glm_gof_stat(data, name, parts)
    out <- structure(list(statistic = c(`X-squared` = statistic),
      parameter = c(df = parts$df),
      p.value = stats::pchisq(statistic, parts$df, lower.tail = FALSE),
      method = paste("GLM", name, "goodness-of-fit test"),
      test.method = name, family = data$family, link = data$link,
      data.name = paste(deparse(stats::formula(fit)), collapse = " "),
      groups = parts$groups), class = "htest")
    attr(out, "gof.diagnostics") <- check
    if (asymcheck) out$diagnostics <- check
    out
  })
  if (length(results) == 1L) return(results[[1L]])
  structure(stats::setNames(results, method),
            class = c("ibist_glm_gof_tests", "list"))
}

.glm_gof_basis <- function(x) {
  if (!ncol(x)) return(matrix(0, nrow(x), 0))
  scale <- sqrt(colSums(x^2))
  x <- x[, scale > 0, drop = FALSE]
  if (!ncol(x)) return(matrix(0, nrow(x), 0))
  x <- sweep(x, 2L, scale[scale > 0], "/")
  s <- svd(x, nv = 0L)
  s$u[, s$d > 1e-10 * max(s$d), drop = FALSE]
}

.glm_gof_data <- function(fit) {
  if (!inherits(fit, "glm") || !isTRUE(fit$converged)) {
    stop("fit must be a converged glm object")
  }
  family <- fit$family$family
  if (!family %in% c("binomial", "poisson")) {
    stop("only binomial and Poisson GLMs with dispersion one are supported")
  }
  mu <- fit$fitted.values
  design <- fit[["x", exact = TRUE]]
  if (is.null(design) && is.null(fit[["model", exact = TRUE]])) {
    stop("fit must retain the model frame or design matrix")
  }
  x <- if (!is.null(design)) design else stats::model.matrix(fit)
  n <- fit$prior.weights
  y <- fit$y
  if (is.null(y)) stop("fit must retain the response (y = TRUE)")
  offset <- fit$offset
  if (is.null(offset)) offset <- rep(0, length(y))
  active.mu <- mu[n > 0]
  if (any(!is.finite(c(x, y, n, active.mu, offset))) || any(n < 0) ||
      any(active.mu <= 0) ||
      (family == "binomial" && any(active.mu >= 1))) {
    stop("invalid retained inputs or boundary fitted means")
  }
  observed <- y * n
  counts <- if (family == "binomial") c(n, observed, n-observed) else
    c(y, n)
  if (any(counts < 0) || any(abs(counts - round(counts)) > 1e-7) ||
      sum(n) > 2^53 - 1 || sum(observed) > 2^53 - 1) {
    stop("weights must encode integer counts/frequencies, not analytic weights")
  }
  keep <- n > 0
  x <- x[keep, , drop = FALSE]
  offset <- offset[keep]
  n <- n[keep]
  observed <- observed[keep]
  expected <- n * mu[keep]
  if (!length(n)) stop("no positive-weight observations")
  key <- apply(cbind(x, offset), 1L, function(z)
    paste(formatC(z, digits = 17, format = "g"), collapse = ":"))
  id <- match(key, unique(key))
  aggregate <- function(z) as.vector(rowsum(z, id, reorder = FALSE))
  first <- !duplicated(key)
  out <- list(family = family, link = fit$family$link,
    x = x[first, , drop = FALSE], observed = aggregate(observed),
    expected = aggregate(expected), n = aggregate(n))
  out$mu <- out$expected / out$n
  out$variance <- if (family == "binomial")
    out$n * out$mu * (1-out$mu) else out$expected
  out$q <- .glm_gof_basis(sqrt(out$variance) * out$x)
  out$rank <- ncol(out$q)
  if (out$rank != fit$rank) {
    stop("numerical model rank differs from fitted rank")
  }
  out
}

.glm_gof_methods <- function(method, g, data) {
  allowed <- c("pearson", "deviance", "hosmer-lemeshow", "osius-rojek",
               "stukel")
  if (!is.character(method) || !length(method) || anyNA(method) ||
      any(!method %in% allowed)) stop("unsupported GOF method")
  if (!is.numeric(g) || length(g) != 1L || !is.finite(g) || g < 3 ||
      g != floor(g)) stop("g must be an integer >= 3")
  if (g != 10 && !"hosmer-lemeshow" %in% method) {
    stop("g applies only to Hosmer-Lemeshow")
  }
  if (any(!method %in% c("pearson", "deviance")) &&
      (data$family != "binomial" || data$link != "logit")) {
    stop("requested method requires a binomial logit model")
  }
  method
}

.glm_gof_hl <- function(data, g) {
  # Projecting a constant checks for an intercept despite aliased columns.
  q <- .glm_gof_basis(data$x)
  if (max(abs(1 - q %*% colSums(q))) > 1e-7) {
    stop("Hosmer-Lemeshow requires an intercept")
  }
  order <- order(data$mu)
  p <- signif(data$mu[order], 14)
  id <- match(p, unique(p))
  agg <- function(z) as.vector(rowsum(z[order], id, reorder = FALSE))
  blocks <- data.frame(observed = agg(data$observed),
                       expected = agg(data$expected), n = agg(data$n))
  target <- floor(sum(blocks$n) / g + 0.5)
  group <- integer(nrow(blocks))
  current <- 1L
  size <- 0
  for (i in seq_len(nrow(blocks))) {
    f <- blocks$n[i]
    if (i > 1L && current < g &&
        !(size < target && size + floor(f/2) <= target)) {
      current <- current + 1L
      size <- 0
    }
    group[i] <- current
    size <- size + f
  }
  if (current > 1L && size < sum(blocks$n)/(2*g)) {
    group[group == current] <- current - 1L
  }
  as.data.frame(rowsum(as.matrix(blocks), group, reorder = FALSE))
}

.glm_gof_parts <- function(data, method, g, validate = TRUE) {
  groups <- data.frame(observed = data$observed, expected = data$expected,
                        n = data$n)
  parts <- list(groups = groups, df = nrow(groups) - data$rank)
  if (method == "hosmer-lemeshow") {
    parts$groups <- .glm_gof_hl(data, g)
    parts$df <- nrow(parts$groups) - 2L
  }
  if (method == "osius-rojek") {
    v <- (1-2*data$mu) / sqrt(data$variance)
    residual <- v - data$q %*% crossprod(data$q, v)
    parts$variance <- 2*sum(1-1/data$n) + sum(residual^2)
    if (validate &&
        (!is.finite(parts$variance) || parts$variance <= 1e-12)) {
      stop("Osius-Rojek standardizing variance is zero or unstable")
    }
    parts$df <- 1L
  }
  if (method == "stukel") {
    eta <- stats::qlogis(data$mu)
    z <- cbind(eta^2 * (eta >= 0), eta^2 * (eta < 0))
    a <- sqrt(data$variance) * z
    a <- a - data$q %*% crossprod(data$q, a)
    # Discard numerical projection noise before column scaling.
    original <- sqrt(colSums((sqrt(data$variance) * z)^2))
    small <- sqrt(colSums(a^2)) <= 1e-10 * pmax(original, 1)
    a[, small] <- 0
    parts$q <- .glm_gof_basis(a)
    parts$df <- ncol(parts$q)
  }
  if (validate && parts$df <= 0) {
    stop("no positive residual or added-term degrees of freedom")
  }
  parts
}

.glm_gof_stat <- function(data, method, parts) {
  tab <- parts$groups
  if (method == "stukel") {
    residual <- (data$observed-data$expected) / sqrt(data$variance)
    return(sum(crossprod(parts$q, residual)^2))
  }
  variance <- if (data$family == "binomial")
    tab$expected * (1-tab$expected/tab$n) else tab$expected
  pearson <- sum((tab$observed-tab$expected)^2 / variance)
  if (method == "osius-rojek") {
    return((pearson-nrow(tab))^2 / parts$variance)
  }
  if (method != "deviance") return(pearson)
  contribution <- function(o, e) {
    positive <- o > 0
    sum(o[positive] * log(o[positive]/e[positive]))
  }
  value <- contribution(tab$observed, tab$expected)
  value <- value + if (data$family == "binomial") {
    contribution(tab$n-tab$observed, tab$n-tab$expected)
  } else -sum(tab$observed-tab$expected)
  max(0, 2*value)
}

.glm_gof_check <- function(data, method, parts) {
  tab <- parts$groups
  rows <- .asymcheck_rows(method, "degrees of freedom", parts$df, 0, ">",
                          parts$df > 0)
  if (method %in% c("pearson", "deviance", "hosmer-lemeshow")) {
    values <- min(tab$expected)
    names <- "minimum expected events"
    if (data$family == "binomial") {
      values <- c(values, min(tab$n-tab$expected))
      names <- c(names, "minimum expected nonevents")
    }
    rows <- rbind(rows, .asymcheck_rows(method, names, values, 5, ">=",
                                        values >= 5))
  } else if (method == "osius-rojek") {
    rows <- rbind(rows, .asymcheck_rows(method, "standardizing variance",
      parts$variance, 1e-12, ">",
      is.finite(parts$variance) && parts$variance > 1e-12))
  } else {
    rows <- rbind(rows, .asymcheck_rows(method, "estimable added terms",
      parts$df, 0, ">", parts$df > 0))
  }
  .asymcheck_new("GLM goodness of fit", method, TRUE, rows,
    list(family = data$family, link = data$link, rank = data$rank,
         converged = TRUE, boundary.fitted = FALSE, groups = tab,
         minimum.fitted.margin = min(if (data$family == "binomial")
           pmin(data$mu, 1-data$mu) else data$mu)),
    c("All GOF methods here are asymptotic; passing checks is not a guarantee.",
      if (method == "hosmer-lemeshow")
        "Risk grouping and probability ties affect Hosmer-Lemeshow results.",
      if (method == "stukel" && parts$df == 1)
        "Only one Stukel direction is estimable.",
      if (method == "stukel" && parts$df == 0)
        "No Stukel direction is estimable.",
      if (method %in% c("osius-rojek", "stukel"))
        "Sparse-cell cutoffs do not assess this test's asymptotic accuracy."))
}

.glm_gof_precheck <- function(fit, method = "pearson", g = 10,
                               asymcheck = FALSE) {
  .asymcheck_flag(asymcheck, "asymcheck")
  data <- .glm_gof_data(fit)
  methods <- .glm_gof_methods(method, g, data)
  checks <- lapply(methods, function(name)
    .glm_gof_check(data, name, .glm_gof_parts(data, name, g, validate = FALSE)))
  .glm_gof_combine(checks)
}

.glm_gof_combine <- function(checks) {
  checks <- unname(checks)
  if (length(checks) == 1L) return(checks[[1L]])
  methods <- vapply(checks, `[[`, character(1), "method")
  notes <- setdiff(unique(unlist(lapply(checks, `[[`, "notes"))),
                   "Criteria do not verify all sampling assumptions.")
  rows <- do.call(rbind, lapply(checks, `[[`, "checks"))
  rownames(rows) <- NULL
  .asymcheck_new("GLM goodness of fit", methods, TRUE, rows,
    list(methods = stats::setNames(lapply(checks, `[[`, "details"), methods)),
    notes)
}

#' @rdname glm.gof
#' @param x An ordered collection of GLM goodness-of-fit tests.
#' @param ... Reserved for future extensions.
#' @export
as.data.frame.ibist_glm_gof_tests <- function(x, ...) {
  data.frame(method = names(x),
    statistic = vapply(x, function(z) unname(z$statistic), numeric(1)),
    df = vapply(x, function(z) unname(z$parameter), numeric(1)),
    p.value = vapply(x, `[[`, numeric(1), "p.value"), row.names = NULL)
}

#' @rdname glm.gof
#' @export
print.ibist_glm_gof_tests <- function(x, ...) {
  print(as.data.frame(x), row.names = FALSE)
  cat("Tests use the same fit; p-values are unadjusted.\n")
  invisible(x)
}
