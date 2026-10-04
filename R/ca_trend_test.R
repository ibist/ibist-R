#' Cochran--Armitage Trend Test
#'
#' Tests for a linear trend in a binary outcome across ordered groups using
#' the Cochran--Armitage statistic. The test can use a normal approximation
#' or an exact conditional distribution.
#'
#' @param x An \eqn{R \times 2} matrix or table of non-negative integer
#'   counts. Rows are ordered groups and columns are outcome categories. The
#'   first column is treated as the event of interest.
#' @param scores Optional numeric scores for the ordered groups. Defaults to
#'   \code{1:R}; scores must be finite and strictly increasing.
#' @param alternative Alternative hypothesis: \code{"two.sided"},
#'   \code{"greater"} (increasing trend), or \code{"less"} (decreasing
#'   trend).
#' @param exact Logical; if \code{TRUE}, compute the exact conditional
#'   p-value given both margins. Otherwise, use the asymptotic normal
#'   approximation.
#'
#' @details
#' A positive signed statistic indicates that the event proportion tends to
#' increase with the group scores. The exact test conditions on the row totals
#' and the total number of events and uses the conditional distribution of
#' \eqn{T = \sum_i x_i n_{i1}}. For a two-sided exact test, the p-value is
#' twice the smaller one-sided tail probability, truncated at one. The exact
#' calculation is performed by dynamic programming and can become slow when
#' the table has many rows or large margins.
#'
#' @return An object of class \code{"htest"}. The signed Cochran--Armitage
#'   statistic is returned along with its square, the score statistic
#'   \eqn{T}, the p-value, and the alternative hypothesis.
#'
#' @references
#' Cochran, W. G. (1954). Some methods for strengthening the common
#' chi-squared tests. \emph{Biometrics}, 10(4), 417--451.
#' \doi{10.2307/3001616}
#'
#' Armitage, P. (1955). Tests for linear trends in proportions and
#' frequencies. \emph{Biometrics}, 11(3), 375--386.
#' \doi{10.2307/3001775}
#'
#' Portier, C., and Hoel, D. G. (1984). Type 1 error of trend tests in
#' proportions and the design of cancer screens.
#' \emph{Communications in Statistics - Theory and Methods}, 13(1), 1--14.
#'
#' @examples
#' tab <- matrix(c(106, 109, 128, 130, 99, 95, 76, 74), nrow = 4)
#' ca.trend.test(tab)
#' ca.trend.test(tab, scores = c(1, 2, 4, 8), alternative = "greater")
#' small.tab <- matrix(c(2, 0, 1, 1, 0, 2), nrow = 3, byrow = TRUE)
#' ca.trend.test(small.tab, exact = TRUE)
#'
#' @export
ca.trend.test <- function(
  x,
  scores = NULL,
  alternative = c("two.sided", "greater", "less"),
  exact = FALSE
) {
  alternative <- match.arg(alternative)

  tab <- as.matrix(x)
  if (length(dim(tab)) != 2L || nrow(tab) < 2L || ncol(tab) != 2L ||
      !is.numeric(tab) || anyNA(tab) || any(!is.finite(tab)) ||
      any(tab < 0) || any(tab != floor(tab))) {
    stop("'x' must be an R-by-2 table of non-negative integer counts.",
         call. = FALSE)
  }
  if (!is.logical(exact) || length(exact) != 1L || is.na(exact)) {
    stop("'exact' must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(scores)) {
    scores <- seq_len(nrow(tab))
  } else if (!is.numeric(scores) || length(scores) != nrow(tab) ||
             anyNA(scores) || any(!is.finite(scores)) ||
             any(diff(scores) <= 0)) {
    stop("'scores' must be finite, strictly increasing numeric values,",
         " one for each row of 'x'.", call. = FALSE)
  }

  keep <- rowSums(tab) > 0
  tab <- tab[keep, , drop = FALSE]
  scores <- scores[keep]
  if (nrow(tab) < 2L) {
    stop("'x' must contain at least two non-empty rows.", call. = FALSE)
  }

  row.totals <- rowSums(tab)
  n <- sum(row.totals)
  events <- tab[, 1L]
  if (sum(events) == 0 || sum(events) == n) {
    stop("Both outcome categories must have positive column totals.",
         call. = FALSE)
  }
  pbar <- sum(events) / n
  score.mean <- sum(row.totals * scores) / n
  numerator <- sum(scores * (events - row.totals * pbar))
  denominator <- sqrt(
    pbar * (1 - pbar) *
      sum(row.totals * (scores - score.mean)^2)
  )
  z <- if (denominator == 0) 0 else numerator / denominator
  score.statistic <- sum(scores * events)

  if (exact) {
    distribution <- ca_trend_distribution(
      row.totals,
      scores - scores[1L],
      sum(events)
    )
    conditional.score <- sum((scores - scores[1L]) * events)
    lower <- sum(
      distribution$prob[
        distribution$score <= conditional.score + distribution$tolerance
      ]
    )
    upper <- sum(
      distribution$prob[
        distribution$score >= conditional.score - distribution$tolerance
      ]
    )
    p.value <- switch(
      alternative,
      two.sided = min(1, 2 * min(lower, upper)),
      greater = upper,
      less = lower
    )
  } else {
    p.value <- switch(
      alternative,
      two.sided = 2 * stats::pnorm(-abs(z)),
      greater = stats::pnorm(z, lower.tail = FALSE),
      less = stats::pnorm(z)
    )
  }

  structure(
    list(
      statistic = c(`Z` = z),
      score.statistic = c(`T` = score.statistic),
      chisq = c(`CA X-squared` = z^2),
      p.value = p.value,
      alternative = alternative,
      method = if (exact) {
        "Exact conditional Cochran-Armitage trend test"
      } else {
        "Asymptotic Cochran-Armitage trend test"
      },
      data.name = deparse(substitute(x))
    ),
    class = "htest"
  )
}

ca_trend_distribution <- function(row.totals, scores, events) {
  states <- data.frame(cases = 0, score = 0, log.weight = 0)

  for (i in seq_along(row.totals)) {
    later.total <- if (i < length(row.totals)) {
      sum(row.totals[(i + 1L):length(row.totals)])
    } else {
      0
    }
    values <- seq.int(0, min(row.totals[i], events))
    updates <- lapply(values, function(value) {
      cases <- states$cases + value
      keep <- cases <= events & cases + later.total >= events
      if (!any(keep)) return(NULL)
      data.frame(
        cases = cases[keep],
        score = states$score[keep] + scores[i] * value,
        log.weight = states$log.weight[keep] +
          lchoose(row.totals[i], value)
      )
    })
    states <- do.call(rbind, Filter(Negate(is.null), updates))
    states$score.key <- format(
      signif(states$score, 15),
      scientific = TRUE,
      trim = TRUE
    )
    states <- aggregate(
      log.weight ~ cases + score.key,
      data = states,
      FUN = log_sum_exp
    )
    names(states)[names(states) == "log.weight"] <- "log.weight"
    states$score <- as.numeric(states$score.key)
    states$score.key <- NULL
  }

  states <- states[states$cases == events, , drop = FALSE]
  log.normalizer <- lchoose(sum(row.totals), events)
  log.prob <- states$log.weight - log.normalizer
  probability <- exp(log.prob)
  probability <- probability / sum(probability)
  tolerance <- 100 * .Machine$double.eps * max(1, abs(states$score))
  list(score = states$score, prob = probability, tolerance = tolerance)
}
