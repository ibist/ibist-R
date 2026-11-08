#' Power Calculation for a One-Sample Poisson Rate Test
#'
#' Computes normal-approximation power or solves one missing design quantity
#' for the score test of a Poisson rate against a positive reference rate.
#'
#' @param n Sample size; total exposure is \code{n * exposure}. Supplied and
#'   solved sample sizes remain continuous. Round up for subject recruitment.
#' @param r0 Positive null rate per unit exposure.
#' @param r1 Positive alternative rate per unit exposure.
#' @param power Target power in (0, 1); omit to compute power.
#' @param sig.level Significance level in (0, 1). Set to \code{NULL} to solve.
#' @param exposure Positive exposure per subject, in the units of the rates.
#' @param alternative One of \code{"two.sided"}, \code{"less"}, or
#'   \code{"greater"}.
#' @param correct Logical; apply a half-count continuity correction.
#'   Default \code{FALSE}, as in the proportion power functions.
#' @param strict Logical; include both tails for a two-sided calculation.
#'   Default \code{FALSE} uses the dominant tail, as in proportion power.
#' @param tol Positive root tolerance on the log or logit scale.
#'
#' @details
#' Exactly one of \code{n}, \code{r0}, \code{r1}, \code{power}, and
#' \code{sig.level} must be \code{NULL}. Exposure is a fixed design input.
#' For total exposure \eqn{T}, the rate estimator has approximate mean
#' \eqn{r_1} and variance \eqn{r_1/T}. The score rejection boundaries use
#' \eqn{r_0 \pm \{z\sqrt{r_0/T}+c\}}, where \eqn{c=0.5/T} with
#' correction and zero otherwise. The critical value uses half the level
#' for two-sided tests. Power uses the alternative normal distribution.
#'
#' For two-sided rate inversion, seek \code{r1 > r0} or \code{r0 < r1}
#' first, then the opposite branch. Directional inversions prefer the branch
#' matching the alternative. Positive roots are bracketed within 100
#' doublings or halvings of 1 (for n) or the other rate. Unattainable targets
#' raise an error. Equal rates cannot determine a sample size.
#' This is a large-sample planning approximation to \code{rate.test()},
#' not exact Poisson power; accuracy can be poor with few expected events.
#'
#' @return An object of class \code{"power.htest"} with \code{n}, rates,
#'   \code{sig.level}, \code{power}, \code{exposure}, \code{alternative},
#'   \code{method}, and \code{note}. Inverse calls retain the target power.
#' @references Minitab, LLC. Example of Power and Sample Size for 1-Sample
#'   Poisson Rate; Methods and Formulas for the same procedure. Minitab
#'   Support, accessed 2026-11-08. \url{https://support.minitab.com}.
#'   Published powers 0.842947 and 0.898200 are reproduced for n = 25 and
#'   30, r0 = 15, r1 = 13, level .05, and a less-than alternative.
#'
#'   R Core Team. The Poisson Distribution. R stats documentation.
#'   \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Poisson.html}
#'   Documents the Poisson mean and variance used in the derivation above.
#' @seealso \code{\link{rate.test}}, \code{\link{power.r2s.test}},
#'   \code{\link{power.p1s.test}}
#' @examples
#' power.r1s.test(n = 100, r0 = 0.1, r1 = 0.2, exposure = 2)
#' power.r1s.test(r0 = 0.1, r1 = 0.2, power = 0.8, exposure = 2)
#' @export
power.r1s.test <- function(n = NULL, r0 = NULL, r1 = NULL, power = NULL,
                          sig.level = 0.05, exposure = 1,
                          alternative = c("two.sided", "less", "greater"),
                          correct = FALSE, strict = FALSE,
                          tol = .Machine$double.eps^0.5) {
  alternative <- match.arg(alternative)
  values <- list(n = n, r0 = r0, r1 = r1, power = power,
                 sig.level = sig.level)
  .rate_power_validate(values, exposure, correct, strict, tol, 1L)
  evaluate <- function(v) {
    total <- v$n * exposure
    delta <- v$r1 - v$r0
    if (alternative == "less") delta <- -delta
    if (alternative == "two.sided") delta <- abs(delta)
    boundary <- qnorm(v$sig.level /
      if (alternative == "two.sided") 2 else 1, lower.tail = FALSE) *
      sqrt(v$r0 / total) + if (correct) 0.5 / total else 0
    se <- sqrt(v$r1 / total)
    answer <- pnorm((delta - boundary) / se)
    if (strict && alternative == "two.sided") {
      answer <- answer + pnorm((-delta - boundary) / se)
    }
    answer
  }
  missing <- names(values)[vapply(values, is.null, logical(1))]
  other <- if (missing == "r0") r1 else if (missing == "r1") r0 else 1
  prefer <- if (missing == "r0") -1 else 1
  if (alternative == "less") prefer <- -prefer
  if (missing == "n" && r0 == r1) {
    stop("equal rates cannot determine a sample size")
  }
  values <- .rate_power_solve(values, evaluate, other, prefer, tol)
  structure(c(values, list(exposure = exposure, alternative = alternative,
    method = "One-sample Poisson rate power calculation",
    note = .rate_power_note(correct))), class = "power.htest")
}

#' Power Calculation for Two Independent Poisson Rates
#'
#' Computes normal-approximation power or solves one missing design quantity
#' for a score test of equal Poisson rates with unequal sizes and exposures.
#'
#' @inheritParams power.r1s.test
#' @param r1 Positive rate per unit exposure in the first group.
#' @param r2 Positive rate per unit exposure in the second group.
#' @param group.rate Positive sample-size ratio \eqn{n_2/n_1}, default 1.
#' @param exposure Positive exposure per subject: a scalar for both groups
#'   or a vector of length two for separate follow-up durations.
#' @param alternative Either \code{"two.sided"} or \code{"one.sided"}.
#'   One-sided power uses the direction of the supplied rate difference,
#'   matching \code{power.p2s.test()}.
#'
#' @details
#' Exactly one of \code{n}, \code{r1}, \code{r2}, \code{power}, and
#' \code{sig.level} must be \code{NULL}. The first group has size n and
#' the second \code{n * group.rate}. Total exposures are
#' \eqn{T_1=n e_1} and \eqn{T_2=n\,group.rate\,e_2}.
#' For \eqn{d=|r_1-r_2|}, use the expected pooled null rate
#' \eqn{\bar r=(r_1T_1+r_2T_2)/(T_1+T_2)}, null standard error
#' \eqn{s_0=\sqrt{\bar r(1/T_1+1/T_2)}}, and alternative standard error
#' \eqn{s_1=\sqrt{r_1/T_1+r_2/T_2}}. The dominant-tail power is
#' \eqn{\Phi((d-zs_0-c)/s_1)}, with
#' \eqn{c=0.5(1/T_1+1/T_2)} when corrected. Strict two-sided power adds
#' \eqn{\Phi((-d-zs_0-c)/s_1)}. This replaces the random pooled rate
#' in \code{rate.test()} by its expectation for large-sample planning.
#'
#' Rate inversion first seeks \code{r1 < r2} or \code{r2 > r1}, then the
#' opposite branch. Numerical bounds and limitations follow
#' \code{power.r1s.test()}. This is not exact conditional Poisson power.
#'
#' @return A \code{"power.htest"} object as for \code{power.r1s.test()},
#'   with both rates, a length-two \code{exposure}, and \code{group.rate}.
#' @references Minitab, LLC. Example of Power and Sample Size for 2-Sample
#'   Poisson Rate; Methods and Formulas for the same procedure. Minitab
#'   Support, accessed 2026-11-08. \url{https://support.minitab.com}.
#'   With strict = TRUE, baseline rate 35, level .05 and unit exposures,
#'   comparison rates 32 and 38 give powers 0.902793 and 0.902550 at
#'   n = 79 and 86 per group. Target power .9 gives these rounded sizes.
#'
#'   R Core Team. The Poisson Distribution. R stats documentation.
#'   \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Poisson.html}
#'   Documents the Poisson mean and variance used in the derivation above.
#' @seealso \code{\link{rate.test}}, \code{\link{power.r1s.test}},
#'   \code{\link{power.p2s.test}}
#' @examples
#' power.r2s.test(n = 100, r1 = 0.1, r2 = 0.2, group.rate = 2)
#' power.r2s.test(r1 = 0.1, r2 = 0.2, power = 0.8,
#'                group.rate = 2, exposure = c(2, 3))
#' @export
power.r2s.test <- function(n = NULL, r1 = NULL, r2 = NULL,
                          sig.level = 0.05, power = NULL, group.rate = 1,
                          exposure = 1,
                          alternative = c("two.sided", "one.sided"),
                          correct = FALSE, strict = FALSE,
                          tol = .Machine$double.eps^0.5) {
  alternative <- match.arg(alternative)
  values <- list(n = n, r1 = r1, r2 = r2, power = power,
                 sig.level = sig.level)
  .rate_power_validate(values, exposure, correct, strict, tol, c(1L, 2L))
  .rate_power_positive(group.rate, "group.rate")
  exposure <- rep(exposure, length.out = 2L)
  evaluate <- function(v) {
    totals <- v$n * c(1, group.rate) * exposure
    pooled <- sum(c(v$r1, v$r2) * (totals / sum(totals)))
    se0 <- sqrt(pooled * sum(1 / totals))
    se1 <- sqrt(sum(c(v$r1, v$r2) / totals))
    delta <- abs(v$r1 - v$r2)
    boundary <- qnorm(v$sig.level /
      if (alternative == "two.sided") 2 else 1, lower.tail = FALSE) * se0 +
      if (correct) 0.5 * sum(1 / totals) else 0
    answer <- pnorm((delta - boundary) / se1)
    if (strict && alternative == "two.sided") {
      answer <- answer + pnorm((-delta - boundary) / se1)
    }
    answer
  }
  missing <- names(values)[vapply(values, is.null, logical(1))]
  other <- if (missing == "r1") r2 else if (missing == "r2") r1 else 1
  prefer <- if (missing == "r1") -1 else 1
  if (missing == "n" && r1 == r2) {
    stop("equal rates cannot determine a sample size")
  }
  values <- .rate_power_solve(values, evaluate, other, prefer, tol)
  structure(c(values, list(exposure = exposure, group.rate = group.rate,
    alternative = alternative,
    method = "Two-sample Poisson rates power calculation",
    note = paste("n is number in the 1st group;",
                 .rate_power_note(correct)))), class = "power.htest")
}

.rate_power_positive <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0) {
    stop(paste0("'", name, "' must be a single positive finite number"))
  }
}

.rate_power_validate <- function(values, exposure, correct, strict, tol,
                                 lengths) {
  if (sum(vapply(values, is.null, logical(1))) != 1L) {
    stop("exactly one design parameter must be NULL")
  }
  for (name in names(values)) {
    value <- values[[name]]
    if (!is.null(value)) {
      .rate_power_positive(value, name)
      if (name %in% c("power", "sig.level") && value >= 1) {
        stop(paste0("'", name, "' must be in (0, 1)"))
      }
    }
  }
  if (!is.numeric(exposure) || !length(exposure) %in% lengths ||
      any(!is.finite(exposure)) || any(exposure <= 0)) {
    stop("'exposure' must have the required length and positive values")
  }
  for (name in c("correct", "strict")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(paste0("'", name, "' must be TRUE or FALSE"))
    }
  }
  .rate_power_positive(tol, "tol")
}

.rate_power_note <- function(correct) {
  if (correct) "continuity correction applied" else "no continuity correction"
}

.rate_power_solve <- function(values, evaluate, reference, prefer, tol) {
  missing <- names(values)[vapply(values, is.null, logical(1))]
  if (missing == "power") {
    values$power <- evaluate(values)
    return(values)
  }
  transform <- if (missing == "sig.level") stats::plogis else exp
  objective <- function(x) {
    v <- values
    v[[missing]] <- transform(x)
    evaluate(v) - values$power
  }
  center <- if (missing == "sig.level") 0 else log(reference)
  start <- objective(center)
  if (is.finite(start) && start == 0) {
    values[[missing]] <- transform(center)
    return(values)
  }
  for (direction in c(prefer, -prefer)) {
    previous <- center
    previous.value <- start
    for (step in seq_len(100L)) {
      endpoint <- center + direction * step * log(2)
      if (missing == "sig.level" && abs(endpoint) > 36) break
      value <- objective(endpoint)
      if (is.finite(value) && is.finite(previous.value) &&
          (value == 0 || sign(value) != sign(previous.value))) {
        root <- uniroot(objective, sort(c(previous, endpoint)), tol = tol)
        values[[missing]] <- transform(root$root)
        return(values)
      }
      previous <- endpoint
      previous.value <- value
    }
  }
  stop(paste0("could not solve for '", missing, "' with supplied parameters"))
}
