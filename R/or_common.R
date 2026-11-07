#' Inference for Stratified 2 x 2 Tables
#'
#' Computes exact conditional tests and confidence intervals for a common
#' odds ratio, and exact or asymptotic tests of homogeneity, across a
#' collection of 2 x 2 tables.
#'
#' @param x A list of 2 x 2 tables, or a 2 x 2 x K array of non-negative
#'   integer counts. Rows are the two groups and columns are the two outcome
#'   categories.
#' @param theta The null common odds ratio.
#' @param side Interval form: \code{"two.sided"} (default),
#'   \code{"lower"}, or \code{"upper"}. An omitted (\code{NULL}) value
#'   selects two-sided limits. One-sided bounds use the full error probability.
#' @param alternative For \code{or.common.test()}, the alternative hypothesis:
#'   \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#'   For \code{or.common.ci()}, a compatibility alias: \code{"greater"}
#'   selects a lower bound and \code{"less"} an upper bound. Conflicting
#'   \code{side} and \code{alternative} values are errors.
#' @param conf.level Confidence level for the interval.
#' @param algorithm Computational algorithm for the exact confidence interval.
#'   The default, \code{"network"}, uses dynamic programming. The
#'   \code{"convolution"} option is retained as an independent check.
#' @param method Homogeneity test: \code{"zelen"} (default, exact),
#'   \code{"breslow-day"}, or \code{"woolf"}. Select one method per call.
#' @param correct Logical. Apply Tarone's adjustment to the Breslow--Day
#'   statistic. Default is \code{FALSE}. Setting this to \code{TRUE} for
#'   another method is an error. This is not Yates' continuity correction.
#' @param ... Reserved for future extensions.
#'
#' @details
#' Conditional on the margins of every stratum, the upper-left cell in each
#' table has a noncentral hypergeometric distribution. The test combines the
#' strata through the sum of those upper-left cells. The two-sided test uses
#' the Jung et al. (2014) rule
#' \eqn{2\min(p_\mathrm{less},p_\mathrm{greater})},
#' truncated at one. The confidence interval is obtained by inverting these
#' one-sided exact tests.
#'
#' \code{or.homogeneity.test()} tests
#' \eqn{H_0: \theta_1 = \cdots = \theta_K}. Its default, Zelen's exact test,
#' conditions on the observed sum
#' of the upper-left cells, which removes the unspecified common odds ratio
#' from the null distribution. The two-sided p-value uses probability
#' ordering.
#'
#' The Breslow--Day method estimates the common odds ratio by the
#' Mantel--Haenszel estimator. In each stratum it fits cell counts with the
#' observed margins and that odds ratio. If \eqn{a_k} and \eqn{A_k} are the
#' observed and fitted upper-left counts, and \eqn{V_k} is the reciprocal
#' of the sum of reciprocals of the four fitted counts, the statistic is
#' \eqn{\sum_k(a_k-A_k)^2/V_k}. Tarone's adjustment subtracts
#' \eqn{\{\sum_k(a_k-A_k)\}^2/\sum_k V_k}.
#'
#' The Woolf method uses stratum log odds ratios \eqn{l_k} with weights
#' \eqn{w_k=1/\sum_{ij}(1/n_{kij})}. Its statistic is
#' \eqn{\sum_k w_k(l_k-\bar l_w)^2}. To match \code{DescTools::WoolfTest()},
#' if any observed cell is zero, 0.5 is added to every cell in every stratum
#' before computing the log odds ratios and weights. This automatic
#' zero-cell adjustment is separate from the \code{correct} argument.
#'
#' Both asymptotic methods use the upper tail of a chi-squared distribution
#' with \eqn{K-1} degrees of freedom and require at least two strata with
#' positive row and column totals. Breslow--Day additionally requires a
#' finite positive Mantel--Haenszel odds-ratio estimate. A warning is issued
#' if a Breslow--Day fitted count, or a Woolf adjusted observed count, is
#' below 5. These approximations may be unreliable for sparse tables;
#' consider the exact Zelen method. Empty margins are rejected, not dropped.
#'
#' The common-odds-ratio test and interval remain exact conditional
#' procedures. The homogeneity tests assess equality of stratum-specific
#' odds ratios without specifying their common value.
#'
#' @return \code{or.common.test()} returns an object of class \code{"htest"}.
#'   \code{or.common.ci()} returns an object of class \code{"ibist_ci"}.
#'   \code{or.homogeneity.test()} returns an object of class \code{"htest"}.
#'   For the homogeneity test, \code{statistic} is the total upper-left count
#'   and \code{parameter} is the number of strata for Zelen's method; for
#'   asymptotic methods they are the chi-squared statistic and degrees of
#'   freedom. The alternative is always \code{"two.sided"}.
#'
#' @references
#' SAS Institute Inc. (2016). \emph{SAS/STAT 14.2 User's Guide}.
#' The FREQ Procedure, pp. 2819--2820: Breslow-Day Test for Homogeneity
#' of the Odds Ratios. Cary, NC: SAS Institute Inc.
#' \url{https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf}
#'
#' Breslow, N. E., and Day, N. E. (1980). \emph{Statistical Methods in Cancer
#' Research. Volume I: The Analysis of Case-Control Studies}.
#' IARC Scientific Publications No. 32. International Agency for Research
#' on Cancer, Lyon.
#'
#' Tarone, R. E. (1985). On heterogeneity tests based on efficient scores.
#' \emph{Biometrika}, 72, 91--95.
#'
#' Woolf, B. (1955). On estimating the relation between blood group and
#' disease. \emph{Annals of Human Genetics}, 19(4), 251--253.
#' \doi{10.1111/j.1469-1809.1955.tb01348.x}
#'
#' Jung, S. H., Biggerstaff, B. J., and Kim, C. (2014). A class of exact
#' tests for a common odds ratio in stratified 2 x 2 tables.
#' \emph{Biometrical Journal}, 56(1), 97--109.
#' \doi{10.1002/bimj.201300048}
#'
#' Mehta, C. R., Patel, N. R., and Gray, R. (1985). On computing an exact
#' confidence interval for the common odds ratio in several 2 x 2 contingency
#' tables. \emph{Journal of the American Statistical Association}, 80(392),
#' 969--973.
#' \doi{10.1080/01621459.1985.10478212}
#'
#' Zelen, M. (1971). The analysis of several 2 x 2 contingency tables.
#' \emph{Biometrika}, 58(1), 129--137.
#' \doi{10.1093/biomet/58.1.129}
#'
#' @examples
#' tabs <- list(
#'   matrix(c(7, 3, 2, 8), nrow = 2, byrow = TRUE),
#'   matrix(c(4, 6, 1, 9), nrow = 2, byrow = TRUE)
#' )
#' or.common.test(tabs)
#' or.common.ci(tabs)
#' or.homogeneity.test(tabs)
#' or.homogeneity.test(tabs, method = "breslow-day")
#' or.homogeneity.test(tabs, method = "breslow-day", correct = TRUE)
#' or.homogeneity.test(tabs, method = "woolf")
#'
#' @name or.common
NULL

#' @rdname or.common
#' @export
or.common.test <- function(
  x,
  theta = 1,
  alternative = c("two.sided", "greater", "less"),
  ...
) {
  tables <- validate_stratified_tables(x)
  theta <- validate_theta(theta)
  alternative <- match.arg(alternative)
  observed <- sum(vapply(
    tables, function(tab) tab[1L, 1L], numeric(1)
  ))
  distribution <- common_or_distribution(tables, theta)
  totals <- distribution$min.total + seq_along(distribution$prob) - 1L
  p_less <- sum(distribution$prob[totals <= observed])
  p_greater <- sum(distribution$prob[totals >= observed])
  p_value <- switch(
    alternative,
    two.sided = min(1, 2 * min(p_less, p_greater)),
    greater = p_greater,
    less = p_less
  )

  structure(
    list(
      statistic = c("sum of upper-left cells" = observed),
      parameter = c("common odds ratio" = theta),
      p.value = p_value,
      null.value = c("common odds ratio" = theta),
      alternative = alternative,
      method = "Exact conditional test for a common odds ratio",
      data.name = deparse(substitute(x)),
      estimate = c(
        "common odds ratio" = conditional_or_mle(tables, observed)
      )
    ),
    class = "htest"
  )
}

#' @rdname or.common
#' @export
or.common.ci <- function(
  x,
  conf.level = 0.95,
  algorithm = c("network", "convolution"),
  side = NULL, alternative = NULL, ...
) {
  side <- resolve_ci_side(side, alternative)$side
  tables <- validate_stratified_tables(x)
  algorithm <- match.arg(algorithm)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }

  observed <- sum(vapply(
    tables, function(tab) tab[1L, 1L], numeric(1)
  ))
  support <- common_or_support(tables)
  alpha <- ci_tail(conf.level, side)
  lower <- if (side == "upper") 0 else common_or_ci_limit(
    tables, observed, support, alpha, lower = TRUE, algorithm = algorithm
  )
  upper <- if (side == "lower") Inf else common_or_ci_limit(
    tables, observed, support, alpha, lower = FALSE, algorithm = algorithm
  )

  new_ibist_ci(
    method = "exact",
    estimate = c(
      "common odds ratio" = conditional_or_mle(tables, observed)
    ),
    intervals = list(c(lower, upper)),
    conf.level = conf.level,
    parameter = "common odds ratio",
    side = side,
    data.name = deparse(substitute(x)),
    null.value = c("common odds ratio" = 1)
  )
}

#' @rdname or.common
#' @export
or.homogeneity.test <- function(
  x,
  method = c("zelen", "breslow-day", "woolf"),
  correct = FALSE,
  ...
) {
  method <- match.arg(method)
  if (!is.logical(correct) || length(correct) != 1L || is.na(correct)) {
    stop("'correct' must be TRUE or FALSE.", call. = FALSE)
  }
  if (correct && method != "breslow-day") {
    stop("'correct = TRUE' is only available for the Breslow-Day method.",
         call. = FALSE)
  }
  tables <- validate_stratified_tables(x)
  if (method != "zelen") {
    return(asymptotic_or_homogeneity(
      tables, method, correct, paste(deparse(substitute(x)), collapse = "\n")
    ))
  }
  observed <- sum(vapply(
    tables, function(tab) tab[1L, 1L], numeric(1)
  ))
  p_value <- zelen_p_value(tables, observed)

  structure(
    list(
      statistic = c("sum of upper-left cells" = observed),
      parameter = c("number of strata" = length(tables)),
      p.value = p_value,
      alternative = "two.sided",
      method = "Zelen's exact test of homogeneity of odds ratios",
      data.name = deparse(substitute(x))
    ),
    class = "htest"
  )
}

asymptotic_or_homogeneity <- function(tables, method, correct, data.name) {
  if (length(tables) < 2L) {
    stop("Asymptotic homogeneity tests require at least two strata.",
         call. = FALSE)
  }
  counts <- vapply(tables, as.numeric, numeric(4))
  valid <- vapply(tables, function(tab) {
    all(is.finite(tab)) && is.finite(sum(tab)) &&
      all(rowSums(tab) > 0) && all(colSums(tab) > 0)
  }, logical(1))
  if (!all(valid)) {
    stop("Each stratum must have finite counts and positive row and ",
         "column totals.",
         call. = FALSE)
  }
  if (method == "woolf") {
    if (any(counts == 0)) counts <- counts + 0.5
    log.or <- log(counts[1L, ]) + log(counts[4L, ]) -
      log(counts[2L, ]) - log(counts[3L, ])
    weights <- 1 / colSums(1 / counts)
    pooled <- sum(weights * log.or) / sum(weights)
    statistic <- sum(weights * (log.or - pooled)^2)
    sparse <- any(counts < 5)
    description <- "Woolf test of homogeneity of odds ratios"
  } else {
    totals <- colSums(counts)
    theta <- sum(counts[1L, ] * counts[4L, ] / totals) /
      sum(counts[2L, ] * counts[3L, ] / totals)
    if (!is.finite(theta) || theta <= 0) {
      stop("Breslow-Day requires a finite positive Mantel-Haenszel odds ratio;",
           " consider method = 'zelen'.", call. = FALSE)
    }
    fitted <- vapply(seq_along(tables), function(i) {
      # Solve the fixed-margin odds-ratio equation on a probability scale.
      r <- (counts[1L, i] + counts[3L, i]) / totals[i]
      c <- (counts[1L, i] + counts[2L, i]) / totals[i]
      if (theta == 1) {
        a <- r * c
      } else {
        b <- 1 + (theta - 1) * (r + c)
        discriminant <- max(0, b^2 + 4 * (1 - theta) * theta * r * c)
        signed.root <- if (b >= 0) sqrt(discriminant) else -sqrt(discriminant)
        q <- -0.5 * (b + signed.root)
        roots <- c(q / (1 - theta), -theta * r * c / q)
        lower <- max(0, r + c - 1)
        upper <- min(r, c)
        roots <- roots[is.finite(roots) & roots > lower & roots < upper]
        if (length(roots) != 1L) {
          stop("Could not fit positive Breslow-Day cell counts.", call. = FALSE)
        }
        a <- roots[[1L]]
      }
      c(a, c - a, r - a, 1 - r - c + a) * totals[i]
    }, numeric(4))
    if (any(!is.finite(fitted)) || any(fitted <= 0)) {
      stop("Could not fit positive Breslow-Day cell counts.", call. = FALSE)
    }
    variance <- 1 / colSums(1 / fitted)
    difference <- counts[1L, ] - fitted[1L, ]
    statistic <- sum(difference^2 / variance)
    if (correct) statistic <- statistic - sum(difference)^2 / sum(variance)
    sparse <- any(fitted < 5)
    description <- paste0("Breslow-Day test of homogeneity of odds ratios",
                          if (correct) " with Tarone adjustment" else "")
  }
  if (sparse) {
    warning("Chi-squared approximation may be incorrect: ",
            "some cell counts are below 5.",
            call. = FALSE)
  }
  statistic <- max(0, statistic)
  df <- length(tables) - 1L
  structure(list(statistic = c("X-squared" = statistic), parameter = c(df = df),
                 p.value = stats::pchisq(statistic, df, lower.tail = FALSE),
                 alternative = "two.sided", method = description,
                 data.name = data.name), class = "htest")
}

validate_theta <- function(theta) {
  invalid <- !is.numeric(theta) || length(theta) != 1L ||
    !is.finite(theta) || theta <= 0
  if (invalid) {
    stop("'theta' must be a single finite positive number.")
  }
  theta
}

validate_stratified_tables <- function(x) {
  dimensions <- dim(x)
  is_table_array <- is.array(x) && length(dimensions) == 3L &&
    identical(dimensions[1:2], c(2L, 2L))
  tables <- if (is_table_array) {
    lapply(seq_len(dimensions[3L]), function(i) x[, , i])
  } else if (is.list(x)) {
    x
  } else {
    stop("'x' must be a list of 2 x 2 tables or a 2 x 2 x K array.")
  }

  if (!length(tables)) {
    stop("'x' must contain at least one 2 x 2 table.")
  }
  lapply(tables, validate_2x2_table)
}

common_or_support <- function(tables) {
  limits <- lapply(tables, function(tab) {
    row1 <- sum(tab[1L, ])
    row2 <- sum(tab[2L, ])
    col1 <- sum(tab[, 1L])
    c(max(0L, col1 - row2), min(row1, col1))
  })
  c(
    sum(vapply(limits, `[`, numeric(1), 1L)),
    sum(vapply(limits, `[`, numeric(1), 2L))
  )
}

stratum_or_distribution <- function(tab, eta) {
  row1 <- sum(tab[1L, ])
  row2 <- sum(tab[2L, ])
  col1 <- sum(tab[, 1L])
  support <- seq.int(max(0L, col1 - row2), min(row1, col1))
  log_weights <- lchoose(row1, support) + lchoose(row2, col1 - support) +
    support * eta
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  list(min = support[1L], prob = weights / sum(weights))
}

common_or_distribution <- function(
  tables,
  theta,
  algorithm = c("network", "convolution")
) {
  algorithm <- match.arg(algorithm)
  if (algorithm == "network") {
    return(common_or_distribution_network(tables, theta))
  }

  distributions <- lapply(tables, stratum_or_distribution, eta = log(theta))
  result <- list(min.total = 0, prob = 1)
  for (distribution in distributions) {
    result$prob <- as.vector(stats::convolve(result$prob,
      rev(distribution$prob), type = "open"))
    result$min.total <- result$min.total + distribution$min
  }
  result
}

common_or_distribution_network <- function(tables, theta) {
  distributions <- lapply(tables, stratum_or_log_weights, eta = log(theta))
  result <- list(min.total = 0, log.weight = 0)

  for (distribution in distributions) {
    new_log_weight <- rep(
      -Inf, length(result$log.weight) + length(distribution$log.weight) - 1L
    )
    for (i in seq_along(result$log.weight)) {
      for (j in seq_along(distribution$log.weight)) {
        index <- i + j - 1L
        new_log_weight[index] <- log_space_add(
          new_log_weight[index],
          result$log.weight[i] + distribution$log.weight[j]
        )
      }
    }
    result$log.weight <- new_log_weight
    result$min.total <- result$min.total + distribution$min
  }

  result$prob <- exp(result$log.weight - log_sum_exp(result$log.weight))
  result$log.weight <- NULL
  result
}

stratum_or_log_weights <- function(tab, eta) {
  row1 <- sum(tab[1L, ])
  row2 <- sum(tab[2L, ])
  col1 <- sum(tab[, 1L])
  support <- seq.int(max(0L, col1 - row2), min(row1, col1))
  log_weight <- lchoose(row1, support) + lchoose(row2, col1 - support) +
    support * eta
  list(min = support[1L], log.weight = log_weight)
}

log_space_add <- function(x, y) {
  if (is.infinite(x)) {
    return(y)
  }
  if (is.infinite(y)) {
    return(x)
  }
  maximum <- max(x, y)
  maximum + log(exp(x - maximum) + exp(y - maximum))
}

conditional_or_mle <- function(tables, observed) {
  support <- common_or_support(tables)
  if (observed <= support[1L]) {
    return(0)
  }
  if (observed >= support[2L]) {
    return(Inf)
  }

  expected_minus_observed <- function(eta) {
    distributions <- lapply(tables, stratum_or_distribution, eta = eta)
    expected <- sum(vapply(distributions, function(distribution) {
      offsets <- seq_along(distribution$prob) - 1L
      distribution$min + sum(offsets * distribution$prob)
    }, numeric(1)))
    expected - observed
  }
  root <- stats::uniroot(
    expected_minus_observed, c(-700, 700), tol = 1e-10
  )$root
  exp(root)
}

common_or_ci_limit <- function(
  tables, observed, support, tail, lower, algorithm
) {
  if (lower && observed <= support[1L]) {
    return(0)
  }
  if (!lower && observed >= support[2L]) {
    return(Inf)
  }

  tail_probability <- function(eta) {
    distribution <- common_or_distribution(tables, exp(eta), algorithm)
    totals <- distribution$min.total + seq_along(distribution$prob) - 1L
    if (lower) {
      sum(distribution$prob[totals >= observed])
    } else {
      sum(distribution$prob[totals <= observed])
    }
  }
  root <- stats::uniroot(
    function(eta) tail_probability(eta) - tail,
    c(-700, 700),
    tol = 1e-10
  )$root
  exp(root)
}

zelen_p_value <- function(tables, observed) {
  supports <- lapply(tables, function(tab) {
    row1 <- sum(tab[1L, ])
    row2 <- sum(tab[2L, ])
    col1 <- sum(tab[, 1L])
    seq.int(max(0L, col1 - row2), min(row1, col1))
  })
  log_coefficients <- lapply(seq_along(tables), function(i) {
    tab <- tables[[i]]
    row1 <- sum(tab[1L, ])
    row2 <- sum(tab[2L, ])
    col1 <- sum(tab[, 1L])
    support <- supports[[i]]
    lchoose(row1, support) + lchoose(row2, col1 - support)
  })
  observed_log_weight <- sum(vapply(seq_along(tables), function(i) {
    match_observed <- tables[[i]][1L, 1L]
    log_coefficients[[i]][match(match_observed, supports[[i]])]
  }, numeric(1)))

  all_log_weights <- zelen_log_weights(
    supports, log_coefficients, observed
  )
  log_normalizing_constant <- log_sum_exp(all_log_weights)
  included <- all_log_weights <= observed_log_weight + 1e-12
  sum(exp(all_log_weights[included] - log_normalizing_constant))
}

zelen_log_weights <- function(supports, log_coefficients, target) {
  n_strata <- length(supports)
  min_remaining <- rev(cumsum(rev(vapply(supports, min, numeric(1)))))
  max_remaining <- rev(cumsum(rev(vapply(supports, max, numeric(1)))))
  weights <- list()

  visit <- function(stratum, remaining, log_weight) {
    if (stratum > n_strata) {
      if (remaining == 0) {
        weights[[length(weights) + 1L]] <<- log_weight
      }
      return(invisible(NULL))
    }

    support <- supports[[stratum]]
    coefficient <- log_coefficients[[stratum]]
    for (i in seq_along(support)) {
      value <- support[i]
      new_remaining <- remaining - value
      if (stratum == n_strata || (
        new_remaining >= min_remaining[stratum + 1L] &&
          new_remaining <= max_remaining[stratum + 1L]
      )) {
        visit(stratum + 1L, new_remaining, log_weight + coefficient[i])
      }
    }
    invisible(NULL)
  }

  visit(1L, target, 0)
  unlist(weights, use.names = FALSE)
}

log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}
