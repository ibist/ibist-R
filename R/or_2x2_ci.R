#' Confidence Intervals for a 2 x 2 Odds Ratio
#'
#' Computes confidence intervals for the odds ratio in a 2 x 2 table.
#'
#' @param x A 2 x 2 table of non-negative integer counts.
#' @param conf.level Confidence level for the interval.
#' @param method Method or methods for confidence interval. One or more of
#'   \code{"wald"}, \code{"gart"}, \code{"cornfield"}, or
#'   \code{"bp"}.
#' @param midp Logical; should the conditional exact method be modified using
#'   mid-p probabilities? The default is \code{FALSE}.
#' @param ... Reserved for future extensions.
#'
#' @details
#' For a table with entries \eqn{n_{11}}, \eqn{n_{12}}, \eqn{n_{21}}, and
#' \eqn{n_{22}}, the sample odds ratio is
#' \eqn{n_{11} n_{22} / (n_{12} n_{21})}.
#'
#' The \code{"wald"} method uses the usual large-sample interval on the
#' log odds-ratio scale. The addition of 0.5 to every cell is the
#' Haldane--Anscombe correction, introduced independently by Haldane (1956)
#' and Anscombe (1956). The \code{"gart"} method applies this correction
#' using the interval construction in Gart (1966). The \code{"cornfield"} method
#' inverts two one-sided conditional tests, while the \code{"bp"} method
#' inverts a two-sided conditional test ordered
#' by null probabilities. Both use the noncentral hypergeometric distribution.
#' Setting \code{midp = TRUE} subtracts half the observed-table probability
#' from the corresponding exact p-value. Exact intervals are conservative;
#' mid-p intervals are not guaranteed to attain at least the nominal coverage
#' level.
#'
#' @return An object of class \code{"ibist_ci"} containing the common estimate
#'   and a confidence-limit matrix with one row per method.
#'
#' @references
#' Anscombe, F. J. (1956). On estimating binomial response relations.
#' \emph{Biometrika}, 43(3--4), 461--464.
#' \doi{10.1093/biomet/43.3-4.461}
#'
#' Baptista, J., and Pike, M. C. (1977). Algorithm AS 115: Exact two-sided
#' confidence limits for the odds ratio in a 2 x 2 table.
#' \emph{Journal of the Royal Statistical Society. Series C}, 26(2), 214--220.
#'
#' Fagerland, M. W., Lydersen, S., and Laake, P. (2015). Recommended confidence
#' intervals for two independent binomial proportions.
#' \emph{Statistical Methods in Medical Research}, 24(2), 224--254.
#'
#' Gart, J. J. (1966). Alternative analyses of contingency tables.
#' \emph{Journal of the Royal Statistical Society: Series B}, 28(1), 164--179.
#'
#' Haldane, J. B. S. (1956). The estimation and significance of the logarithm
#' of a ratio of frequencies. \emph{Annals of Human Genetics}, 20(4), 309--311.
#' \doi{10.1111/j.1469-1809.1955.tb01285.x}
#'
#' @examples
#' tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)
#' or.2x2.ci(tab)
#' or.2x2.ci(tab, method = "gart")
#' or.2x2.ci(tab, method = "cornfield")
#' or.2x2.ci(tab, method = "bp", midp = TRUE)
#'
#' @importFrom stats qnorm uniroot
#' @export
or.2x2.ci <- function(
  x,
  conf.level = 0.95,
  method = c("wald", "gart", "cornfield", "bp"),
  midp = FALSE,
  ...
) {
  tab <- validate_2x2_table(x)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }

  methods <- if (missing(method)) method[1L] else {
    match.arg(method, several.ok = TRUE)
  }

  if (!is.logical(midp) || length(midp) != 1L || is.na(midp)) {
    stop("'midp' must be TRUE or FALSE.")
  }

  ci_methods <- list(
    wald = or_ci_wald,
    gart = or_ci_gart
  )

  estimate <- odds_ratio_estimate(tab)
  intervals <- lapply(
    methods,
    function(method) {
      if (method %in% names(ci_methods)) {
        return(ci_methods[[method]](tab, conf.level, ...))
      }
      if (method == "cornfield") {
        return(or_ci_cornfield(tab, conf.level, midp = midp, ...))
      }
      or_ci_baptista_pike(tab, conf.level, midp = midp, ...)
    }
  )

  method_labels <- paste(
    methods,
    if (midp) "midp" else "exact",
    sep = "-"
  )
  method_labels[methods == "wald"] <- "wald"
  method_labels[methods == "gart"] <- "gart"

  new_ibist_ci(
    method = method_labels,
    estimate = c("odds ratio" = estimate),
    intervals = intervals,
    conf.level = conf.level,
    parameter = "odds ratio",
    data.name = deparse(substitute(x)),
    null.value = c("odds ratio" = 1)
  )
}

validate_2x2_table <- function(x) {
  tab <- as.matrix(x)
  if (!is.numeric(tab) || !identical(dim(tab), c(2L, 2L))) {
    stop("'x' must be a 2 x 2 numeric table.")
  }
  if (any(!is.finite(tab)) || any(tab < 0) || any(tab != floor(tab))) {
    stop("'x' must contain non-negative integer counts.")
  }
  storage.mode(tab) <- "integer"
  tab
}

odds_ratio_estimate <- function(tab) {
  a <- tab[1L, 1L]
  b <- tab[1L, 2L]
  c <- tab[2L, 1L]
  d <- tab[2L, 2L]
  (a * d) / (b * c)
}

or_ci_wald <- function(tab, conf.level, ...) {
  if (any(tab == 0)) {
    stop("Wald confidence interval requires all cell counts to be positive.")
  }

  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  estimate <- odds_ratio_estimate(tab)
  se <- sqrt(sum(1 / tab))

  exp(log(estimate) + c(-1, 1) * z * se)
}

or_ci_gart <- function(tab, conf.level, ...) {
  tab <- tab + 0.5
  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  estimate <- odds_ratio_estimate(tab)
  se <- sqrt(sum(1 / tab))

  exp(log(estimate) + c(-1, 1) * z * se)
}

or_ci_cornfield <- function(tab, conf.level, midp, tol = 1e-8, ...) {
  alpha <- 1 - conf.level
  estimate <- odds_ratio_estimate(tab)

  lower <- if (estimate == 0) {
    0
  } else {
    find_cornfield_limit(
      tab, alpha, lower = TRUE, estimate = estimate, midp = midp, tol = tol
    )
  }

  upper <- if (is.infinite(estimate)) {
    Inf
  } else {
    find_cornfield_limit(
      tab, alpha, lower = FALSE, estimate = estimate, midp = midp, tol = tol
    )
  }

  c(lower, upper)
}

find_cornfield_limit <- function(tab, alpha, lower, estimate, midp, tol) {
  objective <- function(theta) {
    cornfield_tail(tab, theta, lower, midp) - alpha / 2
  }

  if (lower) {
    high <- if (is.finite(estimate)) max(estimate, .Machine$double.eps) else 1
    while (objective(high) < 0 && high < .Machine$double.xmax / 10) {
      high <- high * 10
    }
    low <- high
    while (objective(low) > 0 && low > .Machine$double.xmin * 10) {
      low <- low / 10
    }
    if (objective(low) > 0) {
      return(0)
    }
  } else {
    low <- if (is.finite(estimate)) max(estimate, .Machine$double.eps) else 1
    high <- low
    while (objective(high) > 0 && high < .Machine$double.xmax / 10) {
      high <- high * 10
    }
    if (objective(high) > 0) {
      return(Inf)
    }
  }

  stats::uniroot(objective, interval = c(low, high), tol = tol)$root
}

cornfield_tail <- function(tab, theta, lower, midp) {
  x_obs <- tab[1L, 1L]
  row1 <- sum(tab[1L, ])
  row2 <- sum(tab[2L, ])
  col1 <- sum(tab[, 1L])
  support <- seq.int(max(0L, col1 - row2), min(row1, col1))
  prob <- noncentral_hypergeom_prob(support, row1, row2, col1, theta)
  prob_obs <- noncentral_hypergeom_prob(x_obs, row1, row2, col1, theta)

  tail <- if (lower) {
    sum(prob[support >= x_obs])
  } else {
    sum(prob[support <= x_obs])
  }
  if (midp) tail - 0.5 * prob_obs else tail
}

or_ci_baptista_pike <- function(tab, conf.level, midp, tol = 1e-8, ...) {
  alpha <- 1 - conf.level
  estimate <- odds_ratio_estimate(tab)

  lower <- if (estimate == 0) {
    0
  } else {
    find_bp_limit(
      tab, alpha, lower = TRUE, estimate = estimate, midp = midp, tol = tol
    )
  }

  upper <- if (is.infinite(estimate)) {
    Inf
  } else {
    find_bp_limit(
      tab, alpha, lower = FALSE, estimate = estimate, midp = midp, tol = tol
    )
  }

  c(lower, upper)
}

find_bp_limit <- function(tab, alpha, lower, estimate, midp, tol) {
  objective <- function(theta) {
    bp_value(tab, theta, midp = midp) - alpha
  }

  if (lower) {
    high <- if (is.finite(estimate)) max(estimate, .Machine$double.eps) else 1
    while (objective(high) < 0 && high < .Machine$double.xmax / 10) {
      high <- high * 10
    }
    if (objective(high) < 0) {
      return(0)
    }
    low <- high
    while (objective(low) > 0 && low > .Machine$double.xmin * 10) {
      low <- low / 10
    }
    if (objective(low) > 0)
      return(0)
    interval <- c(low, high)
  } else {
    low <- max(estimate, .Machine$double.eps)
    while (objective(low) < 0 && low > .Machine$double.xmin * 10) {
      low <- low / 10
    }
    high <- low
    while (objective(high) > 0 && high < .Machine$double.xmax / 10) {
      high <- high * 10
    }
    if (objective(high) > 0)
      return(Inf)
    interval <- c(low, high)
  }

  stats::uniroot(objective, interval = interval, tol = tol)$root
}

bp_value <- function(tab, theta, midp) {
  x_obs <- tab[1L, 1L]
  row1 <- sum(tab[1L, ])
  row2 <- sum(tab[2L, ])
  col1 <- sum(tab[, 1L])

  support <- seq.int(max(0L, col1 - row2), min(row1, col1))
  prob <- noncentral_hypergeom_prob(support, row1, row2, col1, theta)
  prob_obs <- noncentral_hypergeom_prob(x_obs, row1, row2, col1, theta)

  tolerance <- max(1e-14 * prob_obs, .Machine$double.eps)
  equal <- abs(prob - prob_obs) <= tolerance
  if (midp) {
    sum(prob[prob < prob_obs & !equal]) + 0.5 * sum(prob[equal])
  } else {
    sum(prob[prob < prob_obs & !equal]) + sum(prob[equal])
  }
}

noncentral_hypergeom_prob <- function(x, row1, row2, col1, theta) {
  support <- seq.int(max(0L, col1 - row2), min(row1, col1))

  log_weights <- lchoose(row1, support) +
    lchoose(row2, col1 - support) +
    support * log(theta)
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  probs <- weights / sum(weights)

  probs[match(x, support)]
}
