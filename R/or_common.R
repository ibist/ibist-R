#' Exact Inference for Stratified 2 x 2 Tables
#'
#' Computes an exact conditional test or confidence interval for a common
#' odds ratio across a collection of 2 x 2 tables.
#'
#' @param x A list of 2 x 2 tables, or a 2 x 2 x K array of non-negative
#'   integer counts. Rows are the two groups and columns are the two outcome
#'   categories.
#' @param theta The null common odds ratio.
#' @param alternative A character string specifying the alternative
#'   hypothesis: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' @param conf.level Confidence level for the interval.
#' @param algorithm Computational algorithm for the exact confidence interval.
#'   The default, \code{"network"}, uses dynamic programming. The
#'   \code{"convolution"} option is retained as an independent check.
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
#' \code{or.homogeneity.test()} performs Zelen's exact test of
#' \eqn{H_0: \theta_1 = \cdots = \theta_K}. It conditions on the observed sum
#' of the upper-left cells, which removes the unspecified common odds ratio
#' from the null distribution. The two-sided p-value uses probability
#' ordering.
#'
#' All three procedures are exact conditional procedures; none is an
#' asymptotic alternative. The first two procedures test or estimate a common
#' odds ratio. The homogeneity procedure tests whether the stratum-specific
#' odds ratios are
#' equal, without specifying their common value.
#'
#' @return \code{or.common.test()} returns an object of class \code{"htest"}.
#'   \code{or.common.ci()} returns an object of class \code{"ibist_ci"}.
#'   \code{or.homogeneity.test()} returns an object of class \code{"htest"}.
#'
#' @references
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
  ...
) {
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
  alpha <- 1 - conf.level
  lower <- common_or_ci_limit(
    tables, observed, support, alpha / 2, lower = TRUE, algorithm = algorithm
  )
  upper <- common_or_ci_limit(
    tables, observed, support, alpha / 2, lower = FALSE, algorithm = algorithm
  )

  new_ibist_ci(
    method = "exact",
    estimate = c(
      "common odds ratio" = conditional_or_mle(tables, observed)
    ),
    intervals = list(c(lower, upper)),
    conf.level = conf.level,
    parameter = "common odds ratio",
    data.name = deparse(substitute(x)),
    null.value = c("common odds ratio" = 1)
  )
}

#' @rdname or.common
#' @export
or.homogeneity.test <- function(x, ...) {
  tables <- validate_stratified_tables(x)
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
