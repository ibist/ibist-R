#' Epidemiological Impact Measures for a 2 x 2 Table
#'
#' Computes unadjusted exposed-group and population attributable measures.
#'
#' @param x A 2 x 2 matrix or table of non-negative integer counts. Rows
#'   are exposed and unexposed; columns are events and non-events.
#' @param design \code{"cohort"} (default) or \code{"case-control"}.
#' @param conf.level Confidence level between zero and one, default 0.95.
#' @param exposure Known population exposure proportion for cohort data.
#' @param prevalence Known population outcome prevalence for case-control
#'   data, strictly between zero and one. It is used to reconstruct risks
#'   by Bayes' rule. Both population arguments are treated as fixed values;
#'   uncertainty in external estimates is not included.
#' @param rare Logical. With case-control data and no outcome prevalence,
#'   TRUE explicitly selects the rare-outcome odds-ratio approximation.
#' @param rr.method Risk-ratio interval used for the cohort exposed fraction,
#'   \code{"katz"} (default) or \code{"koopman"}. Other measures use Wald
#'   intervals on their documented scales.
#' @param measure Optional character vector selecting measures and their
#'   order. Choices are \code{"risk-exposed"}, \code{"risk-population"},
#'   \code{"fraction-exposed"}, and \code{"fraction-population"}. By default,
#'   all four are reported for cohort data and case-control data with known
#'   outcome prevalence. Otherwise the two case-control fractions are reported.
#' @param alternative Compatibility alias for \code{side}, using
#'   \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' @param tol Numerical tolerance for Koopman interval inversion.
#' @param side \code{"two.sided"} (default), \code{"lower"}, or
#'   \code{"upper"}. Conflicting side and alternative values are errors.
#'
#' @details
#' Let the table rows be \eqn{(a,b)} and \eqn{(c,d)}, and define
#' \eqn{p_1=a/(a+b)}, \eqn{p_0=c/(c+d)}, and population exposure
#' proportion \eqn{e}. The measures are
#' \deqn{AR_E=p_1-p_0,\quad AF_E=1-p_0/p_1,}
#' \deqn{PAR=e(p_1-p_0),\quad PAF=1-p_0/[ep_1+(1-e)p_0].}
#' Results use proportions, rather than percentages. Negative values are
#' retained for protective associations. Causal interpretations require
#' comparability and appropriate control of confounding.
#'
#' The exposed risk uses the uncorrected Wald interval from
#' \code{\link{prop.2s.ci}}. Exposed fractions transform risk-ratio limits
#' \eqn{[L,U]} to \eqn{[1-1/L,1-1/U]}. Katz uses log-risk-ratio limits;
#' Koopman uses \code{\link{rr.2x2.ci}} score limits. Katz requires positive
#' event counts. With Koopman, a zero exposed event count gives an undefined
#' point estimate (NA) and extended confidence limits; no pseudo-counts are
#' added. Both-zero event counts give the full interval \eqn{[-\infty,1]}.
#'
#' Population-risk intervals use a multivariate delta-method Wald variance.
#' Population-fraction intervals use a delta-method variance on
#' \eqn{\log(1-PAF)} and transform back. Without a supplied exposure
#' proportion, these population calculations assume a representative sample
#' with four-cell multinomial sampling. The covariance between the estimated
#' overall and unexposed risks is included. Fixed exposure proportions use
#' independent binomial variances conditional on the two row totals.
#' A cohort sampled separately within exposure groups needs a population
#' exposure proportion for population inference. Its enrollment proportions
#' alone do not identify the population exposure distribution.
#'
#' Case-control counts do not identify absolute risk differences. With
#' \code{rare = TRUE}, the exposed fraction uses \eqn{1-1/OR} and transforms
#' the log-Wald OR limits. The population fraction uses
#' \eqn{[a/(a+c)](1-1/OR)}. Its log-complement variance is
#' \eqn{a/[c(a+c)]+b/[d(b+d)]}, accounting for both case and control samples.
#' This approximation requires a rare outcome in the source population;
#' its frequency in the sampled table does not assess that assumption.
#'
#' Population-fraction log-Wald intervals require positive unexposed event
#' counts; case-control calculations also require unexposed controls.
#' Plug-in Wald variances can be degenerate at boundaries. These intervals
#' are asymptotic and do not guarantee finite-sample coverage. Risk intervals
#' are restricted to [-1,1]; fraction intervals have domain [-Inf,1].
#'
#' Known case-control outcome prevalence replaces the rare-outcome
#' approximation. Let \eqn{q_1=a/(a+c)}, \eqn{q_0=b/(b+d)}, and let
#' \eqn{r} be the supplied outcome prevalence. Then
#' \eqn{e=rq_1+(1-r)q_0}, \eqn{p_1=rq_1/e}, and
#' \eqn{p_0=r(1-q_1)/(1-e)}. Delta-method variances propagate both independent
#' case/control binomial variances through these equations, including their
#' induced covariance. Absolute-risk intervals use the risk scale; fraction
#' intervals use the log-complement scale.
#'
#' If a valid table does not permit a requested interval, the summary retains
#' its point estimate where defined, reports NA limits, records the reason
#' under \code{notes}, and issues a warning. Invalid inputs and unsupported
#' study-design combinations are errors. These are computability notes;
#' no asymptotic-adequacy pass is implied.
#'
#' @return An \code{"ibist_impact"} object with named \code{estimate},
#'   \code{conf.int} (lower/upper matrix), \code{method} per measure,
#'   \code{measure}, \code{design}, \code{conf.level}, \code{side}, and
#'   \code{observed}. Population values and whether they were estimated are
#'   retained as \code{exposure}, \code{prevalence},
#'   \code{exposure.estimated}, and \code{prevalence.estimated}.
#'   Approximation assumptions and interval failures are recorded under
#'   \code{approximation} and \code{notes}. \code{as.data.frame()} returns
#'   the summary table. All estimates use proportions, without rounding.
#' @references
#' Hanley, J. A. (2001). A heuristic approach to the formulas for population
#' attributable fraction. Journal of Epidemiology and Community Health,
#' 55, 508--514. \doi{10.1136/jech.55.7.508}.
#'
#' Katz, D., Baptista, J., Azen, S. P., and Pike, M. C. (1978). Obtaining
#' confidence intervals for the risk ratio in cohort studies. Biometrics,
#' 34, 469--474.
#'
#' Koopman, P. A. R. (1984). Confidence intervals for the ratio of two binomial
#' proportions. Biometrics, 40, 513--517.
#'
#' Jewell, N. P. (2004). Statistics for Epidemiology. Chapman & Hall/CRC,
#' Section 7.4 (case-control population fraction, p. 90).
#'
#' Pirikahu, S., Jones, G., and Hazelton, M. L. (2021). Bayesian Inference
#' for Population Attributable Measures from Under-identified Models.
#' Sections 1--3 (definitions and sampling-design identifiability).
#' \url{https://arxiv.org/abs/2105.12901}.
#' @seealso \code{\link{rr.2x2.ci}}, \code{\link{or.2x2.ci}},
#'   \code{\link{prop.2s.ci}}
#' @examples
#' vf <- matrix(c(38, 62, 267, 717), 2, byrow = TRUE)
#' impact.2x2(vf, design = "cohort")
#' impact.2x2(vf, design = "cohort", exposure = 0.1)
#' coffee <- matrix(c(265, 483, 112, 280), 2, byrow = TRUE)
#' impact.2x2(coffee, design = "case-control", rare = TRUE)
#' impact.2x2(coffee, design = "case-control", prevalence = 0.01)
#' @export
impact.2x2 <- function(
  x, design = c("cohort", "case-control"), conf.level = 0.95,
  exposure = NULL, prevalence = NULL, rare = FALSE,
  rr.method = c("katz", "koopman"), measure = NULL,
  alternative = NULL, tol = 1e-8, side = NULL
) {
  data.name <- deparse(substitute(x))
  design <- match.arg(design)
  rr.method <- match.arg(rr.method)
  orientation <- resolve_ci_side(side, alternative)
  side <- orientation$side
  validate_binomial_ci_settings(conf.level, orientation$alternative, tol)
  validate_ci_flag(rare, "rare")
  if (!is.matrix(x) && !inherits(x, "table")) {
    stop("'x' must be a 2 x 2 table of non-negative integer counts.")
  }
  counts <- independent_binomial_counts(x)
  tab <- as.matrix(x)
  storage.mode(tab) <- "double"
  if (sum(tab) > 2^53-1) stop("The table total exceeds exact integer range.")
  for (arg in c("exposure", "prevalence")) {
    value <- get(arg)
    if (!is.null(value) && (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 0 || value > 1)) {
      stop("'", arg, "' must be a known proportion between zero and one.")
    }
  }
  if (design == "cohort" && !is.null(prevalence)) {
    stop("'prevalence' is supported only for case-control data.")
  }
  if (design == "case-control") {
    if (!is.null(exposure)) {
      stop("'exposure' is supported only for cohort data.")
    }
    if (any(colSums(tab) == 0)) {
      stop("Both case/control totals must be positive.")
    }
    if (!is.null(prevalence) && (prevalence == 0 || prevalence == 1)) {
      stop("'prevalence' must be strictly between zero and one.")
    }
    if (rr.method != "katz") {
      stop("'rr.method' applies only to the cohort exposed fraction.")
    }
    if (is.null(prevalence) && !rare) {
      stop("Supply outcome prevalence or explicitly set 'rare = TRUE'.")
    }
  }
  choices <- c("risk-exposed", "risk-population",
               "fraction-exposed", "fraction-population")
  available <- if (design == "case-control" && is.null(prevalence)) {
    choices[3:4]
  } else choices
  if (is.null(measure)) measure <- available
  if (!is.character(measure) || !length(measure) || anyNA(measure) ||
      any(!measure %in% choices)) stop("Unknown or empty 'measure'.")
  if (any(!measure %in% available)) {
    stop("Case-control absolute risks require population outcome prevalence.")
  }
  known <- design == "case-control" && !is.null(prevalence)
  model <- if (known) impact_case_control(tab, prevalence) else NULL
  p <- counts$x/counts$n
  e <- if (is.null(exposure)) counts$n[1L]/sum(tab) else exposure
  if (design == "cohort") {
    risk <- e*p[1L] + (1-e)*p[2L]
    estimates <- c(p[1L]-p[2L], e*(p[1L]-p[2L]),
                   if (p[1L] == 0) NA_real_ else 1-p[2L]/p[1L],
                   if (risk == 0) NA_real_ else 1-p[2L]/risk)
  } else if (known) {
    estimates <- model$estimate
  } else {
    a <- tab[1L, 1L]; b <- tab[1L, 2L]
    c <- tab[2L, 1L]; d <- tab[2L, 2L]
    estimates <- c(NA_real_, NA_real_,
      if (a*d == 0) NA_real_ else 1-b*c/(a*d),
      if (d == 0) NA_real_ else 1-(c/(a+c))/(d/(b+d)))
  }
  names(estimates) <- choices
  notes <- rep(NA_character_, length(measure))
  methods <- vapply(measure, function(m) {
    if (m == "fraction-exposed" && design == "cohort") {
      paste(rr.method, "RR transformation")
    } else if (m == "fraction-exposed" && !known) {
      "Wald OR transformation (rare outcome)"
    } else if (startsWith(m, "fraction")) "log-complement Wald" else "Wald"
  }, character(1))
  intervals <- lapply(seq_along(measure), function(i) {
    m <- measure[i]
    tryCatch({
      if (known) {
        j <- match(m, choices)
        impact_wald(model$center[j], model$variance[j],
                    conf.level, side, startsWith(m, "fraction"))
      } else if (m == "risk-exposed") {
        as.numeric(prop.2s.ci(tab, conf.level = conf.level, method = "wald",
                             correct = FALSE, side = side)$conf.int)
      } else if (m == "fraction-exposed") {
        fit <- if (design == "cohort") {
          rr.2x2.ci(tab, conf.level = conf.level, method = rr.method,
                    side = side, tol = tol)
        } else or.2x2.ci(tab, conf.level = conf.level, method = "wald",
                         side = side)
        as.numeric(1-1/fit$conf.int)
      } else {
        moments <- impact_population(tab, m, design, exposure)
        impact_wald(moments$center, moments$variance, conf.level, side,
                    startsWith(m, "fraction"))
      }
    }, error = function(err) {
      notes[i] <<- conditionMessage(err)
      c(NA_real_, NA_real_)
    })
  })
  names(notes) <- measure
  if (any(!is.na(notes))) {
    warning("Some impact intervals are unavailable; see 'notes'.",
            call. = FALSE)
  }
  limits <- do.call(rbind, intervals)
  dimnames(limits) <- list(measure, c("lower", "upper"))
  attr(limits, "conf.level") <- conf.level
  structure(list(estimate = estimates[measure], conf.int = limits,
    method = unname(methods), measure = measure, design = design,
    conf.level = conf.level, side = side,
    alternative = orientation$alternative, observed = tab,
    data.name = data.name, rr.method = rr.method,
    exposure = if (known) model$exposure else {
      if (design == "cohort") unname(e) else NULL
    },
    prevalence = if (design == "cohort") unname(risk) else prevalence,
    exposure.estimated = if (known) TRUE else {
      if (design == "cohort") is.null(exposure) else NA
    },
    prevalence.estimated = if (design == "cohort") TRUE else {
      if (known) FALSE else NA
    },
    approximation = if (design == "case-control" && !known) {
      "Rare-outcome approximation: sample odds ratio approximates risk ratio."
    } else NULL,
    notes = notes), class = "ibist_impact")
}

impact_wald <- function(center, variance, conf.level, side, fraction) {
  if (!is.finite(center) || !is.finite(variance) || variance < 0) {
    stop("The log-complement or its variance is undefined at this boundary.")
  }
  z <- stats::qnorm(ci_tail(conf.level, side), lower.tail = FALSE)
  limits <- if (fraction) {
    -expm1(center + c(1, -1)*z*sqrt(variance))
  } else center + c(-1, 1)*z*sqrt(variance)
  domain <- if (fraction) c(-Inf, 1) else c(-1, 1)
  ci_bounds(pmax(domain[1L], pmin(domain[2L], limits)), side, domain)
}

impact_case_control <- function(tab, prevalence) {
  totals <- colSums(tab)
  q <- tab[1L, ]/totals
  r <- prevalence
  e <- r*q[1L] + (1-r)*q[2L]
  p1 <- r*q[1L]/e
  p0 <- r*(1-q[1L])/(1-e)
  dp1 <- r*(1-r)*c(q[2L], -q[1L])/e^2
  dp0 <- r*(1-r)*c(-(1-q[2L]), 1-q[1L])/(1-e)^2
  gradients <- rbind(dp1-dp0, -dp0, dp0/p0-dp1/p1, dp0/p0)
  variance <- as.numeric(gradients^2 %*% (q*(1-q)/totals))
  estimate <- c(p1-p0, r-p0,
                if (p1 == 0) NA_real_ else 1-p0/p1, 1-p0/r)
  center <- c(estimate[1:2], log(p0)-log(p1), log(p0)-log(r))
  list(estimate = unname(estimate), center = unname(center),
       variance = variance, exposure = unname(e))
}
impact_population <- function(tab, measure, design, exposure) {
  a <- tab[1L, 1L]
  b <- tab[1L, 2L]
  c <- tab[2L, 1L]
  d <- tab[2L, 2L]
  n <- sum(tab)
  n0 <- c+d
  events <- a+c
  fraction <- measure == "fraction-population"
  if (design == "case-control") {
    if (c == 0 || d == 0) {
      stop("population fraction intervals need unexposed cases and controls")
    }
    center <- log(c/events)-log(d/(b+d))
    return(list(estimate = -expm1(center), center = center,
                variance = a/(c*events) + b/(d*(b+d))))
  }
  n1 <- a+b
  p <- c(a/n1, c/n0)
  e <- if (is.null(exposure)) n1/n else exposure
  total.risk <- e*p[1L] + (1-e)*p[2L]
  if (fraction && (c == 0 || total.risk == 0)) {
    stop("population fraction log-Wald intervals need unexposed events")
  }
  estimate <- e*(p[1L]-p[2L])
  center <- if (fraction) log(p[2L])-log(total.risk) else estimate
  if (!is.null(exposure)) {
    gradient <- if (fraction) c(-e/total.risk,
      1/p[2L]-(1-e)/total.risk) else c(e, -e)
    variance <- sum(gradient^2 * p*(1-p)/c(n1, n0))
  } else {
    # These gradients are with respect to the four multinomial cell counts.
    gradient <- if (fraction) {
      c(1/n-1/events, 1/n, 1/n+1/c-1/n0-1/events, 1/n-1/n0)
    } else c(1/n-events/n^2, -events/n^2,
             1/n-events/n^2-d/n0^2, -events/n^2+c/n0^2)
    # Homogeneity of degree zero makes the multinomial covariance term zero.
    variance <- sum(as.vector(t(tab)) * gradient^2)
  }
  list(estimate = if (fraction) -expm1(center) else estimate,
       center = center, variance = variance)
}

#' @rdname impact.2x2
#' @param row.names Optional row names for the summary data frame.
#' @param optional Logical, passed to \code{as.data.frame()}.
#' @param ... Further arguments passed to the print or data-frame method.
#' @export
as.data.frame.ibist_impact <- function(x, row.names = NULL, optional = FALSE,
                                      ...) {
  data.frame(measure = x$measure, estimate = unname(x$estimate),
             lower = x$conf.int[, 1L], upper = x$conf.int[, 2L],
             method = x$method,
             row.names = if (is.null(row.names)) seq_along(x$measure) else {
               row.names
             }, check.names = !optional, ...)
}

#' @rdname impact.2x2
#' @param digits Number of significant digits to print.
#' @export
print.ibist_impact <- function(x, digits = max(3L, getOption("digits")-3L),
                               ...) {
  cat("\nEpidemiological impact measures (", x$design, ")\n", sep = "")
  cat(100*x$conf.level, "% confidence intervals; ", x$side, "\n", sep = "")
  print(as.data.frame(x), digits = digits, row.names = FALSE, ...)
  if (!is.null(x$approximation)) cat("\n", x$approximation, "\n", sep = "")
  for (i in which(!is.na(x$notes))) {
    cat(x$measure[i], ": ", x$notes[i], "\n", sep = "")
  }
  invisible(x)
}
