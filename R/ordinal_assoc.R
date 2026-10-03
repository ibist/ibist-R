#' Ordinal Association Measures with Asymptotic Confidence Intervals and Tests
#'
#' Computes ordinal association measures and their asymptotic confidence
#' intervals from an RxC table. It also tests the null hypothesis that each
#' association measure is zero.
#'
#' @param x A two-dimensional matrix or table of non-negative integer counts.
#'   Rows and columns must be ordered from the lowest to highest category.
#' @param measure One or more measures: \code{"gamma"}, \code{"tau-b"},
#'   \code{"tau-c"}, \code{"somers-c|r"}, \code{"somers-r|c"},
#'   \code{"pearson"}, and \code{"spearman"}.
#'   By default, all seven are computed.
#' @param conf.level Confidence level for the two-sided asymptotic confidence
#'   intervals.
#' @param alternative Alternative hypothesis for the tests:
#'   \code{"two.sided"}, \code{"greater"}, or \code{"less"}. Positive
#'   association means that higher row categories tend to occur with higher
#'   column categories.
#' @param row.scores,col.scores Optional numeric scores for row and column
#'   categories when computing Pearson correlation. By default, numeric
#'   category labels are used; otherwise, category order numbers are used.
#'   Spearman correlation always uses marginal midrank scores.
#'
#' @details
#' The procedure computes the concordant and discordant pair counts,
#' \eqn{n_C} and \eqn{n_D}, and sets \eqn{S=n_C-n_D}. Its null variance is
#' \deqn{\mathrm{Var}_0(S) = \sum_{ij} n_{ij}(A_{ij}-D_{ij})^2 - 4S^2/n,}
#' where \eqn{A_{ij}} and \eqn{D_{ij}} are the concordant and discordant
#' contributions for cell \eqn{(i,j)}. The test statistic is
#' \eqn{Z=S/\sqrt{\mathrm{Var}_0(S)}}.
#'
#' The concordance-based measures share a test statistic, but Pearson and
#' Spearman correlations use their own multinomial null variances. Standard
#' errors and confidence intervals are computed separately for each measure.
#' Confidence intervals are two-sided regardless of the test alternative.
#' Rows and columns with zero marginal totals are dropped before computation.
#'
#' The \code{"somers-c|r"} measure is Somers' \eqn{\Delta(C|R)}, treating the
#' row variable as explanatory; \code{"somers-r|c"} is
#' \eqn{\Delta(R|C)}, treating the column variable as explanatory.
#'
#' Pearson correlation uses the supplied category scores or the default
#' category scores described above (numeric labels such as 1, 2, and 10 are
#' used as-is). Spearman correlation uses marginal
#' midranks, so observations tied within a category receive the same rank.
#' The correlation standard errors and null variances use multinomial
#' contingency-table formulas. The confidence intervals are Wald intervals.
#'
#' @return An object of class \code{"ibist_ordinal_assoc"} containing:
#' \describe{
#'   \item{estimate}{Named vector of selected association estimates.}
#'   \item{std.error}{Named vector of asymptotic standard errors.}
#'   \item{conf.int}{Matrix of lower and upper confidence limits.}
#'   \item{statistic}{Named asymptotic test statistics, one per selected
#'     measure.}
#'   \item{p.value}{Named p-values, one per selected measure.}
#'   \item{S, var.S}{The observed concordance contrast and its null variance.}
#'   \item{concordant, discordant}{Concordant and discordant pair counts.}
#' }
#'
#' @references
#' Brown, M. B., and Benedetti, J. K. (1977). Sampling behavior of tests for
#' correlation in two-way contingency tables. \emph{Journal of the American
#' Statistical Association}, 72(358), 309--315.
#'
#' Goodman, L. A., and Kruskal, W. H. (1954). Measures of association for
#' cross classifications. \emph{Journal of the American Statistical
#' Association}, 49(268), 732--764.
#'
#' Kendall, M. G. (1970). \emph{Rank Correlation Methods} (4th ed.). Griffin.
#'
#' Somers, R. H. (1962). A new asymmetric measure of association for ordinal
#' variables. \emph{American Sociological Review}, 27(6), 799--811.
#'
#' Stuart, A. (1953). The estimation and comparison of strengths of association
#' in contingency tables. \emph{Biometrika}, 40(1--2), 105--110.
#'
#' @examples
#' tab <- matrix(c(12, 3, 1, 4, 9, 2, 1, 5, 13), nrow = 3, byrow = TRUE)
#' ordinal.assoc(tab)
#' ordinal.assoc(
#'   tab,
#'   measure = c("tau-b", "somers-c|r", "pearson", "spearman")
#' )
#'
#' @export
ordinal.assoc <- function(
  x,
  measure = c(
    "gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c",
    "pearson", "spearman"
  ),
  conf.level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  row.scores = NULL,
  col.scores = NULL
) {
  all.measures <- c(
    "gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c",
    "pearson", "spearman"
  )
  measures <- if (missing(measure)) {
    all.measures
  } else {
    match.arg(
      measure,
      choices = all.measures,
      several.ok = TRUE
    )
  }
  alternative <- match.arg(alternative)

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }

  tab <- as.matrix(x)
  if (length(dim(tab)) != 2L || any(dim(tab) < 2L) || !is.numeric(tab) ||
      anyNA(tab) || any(!is.finite(tab)) || any(tab < 0) ||
      any(tab != floor(tab))) {
    stop("'x' must be a two-dimensional table of non-negative integer counts.")
  }

  keep.rows <- rowSums(tab) > 0
  keep.cols <- colSums(tab) > 0
  if (!is.null(row.scores)) {
    row.scores <- validate_ordinal_scores(row.scores, nrow(tab), "row.scores")
    row.scores <- row.scores[keep.rows]
  }
  if (!is.null(col.scores)) {
    col.scores <- validate_ordinal_scores(col.scores, ncol(tab), "col.scores")
    col.scores <- col.scores[keep.cols]
  }
  tab <- tab[keep.rows, keep.cols, drop = FALSE]
  if (any(dim(tab) < 2L)) {
    stop("'x' must have at least two non-empty rows and columns.")
  }

  if (is.null(row.scores)) {
    row.scores <- default_ordinal_scores(rownames(tab), nrow(tab))
  }
  if (is.null(col.scores)) {
    col.scores <- default_ordinal_scores(colnames(tab), ncol(tab))
  }

  components <- ordinal_assoc_components(tab)
  n <- sum(tab)
  margin <- n * (n - 1) / 2
  u_row <- margin - sum(choose(rowSums(tab), 2))
  u_col <- margin - sum(choose(colSums(tab), 2))
  u <- sqrt(u_row * u_col)
  d <- min(dim(tab))

  all.estimates <- c(
    gamma = if (components$concordant + components$discordant > 0) {
      components$S / (components$concordant + components$discordant)
    } else {
      NA_real_
    },
    `tau-b` = components$S / u,
    `tau-c` = 2 * d * components$S / (n^2 * (d - 1)),
    `somers-c|r` = components$S / u_row,
    `somers-r|c` = components$S / u_col,
    pearson = NA_real_,
    spearman = NA_real_
  )

  A <- components$concordant.by.cell
  D <- components$discordant.by.cell
  d_cell <- A - D
  v_cell <- outer(rowSums(tab), rep(1, ncol(tab))) * u_col +
    outer(rep(1, nrow(tab)), colSums(tab)) * u_row
  var0S <- components$var.S

  all.variances <- c(
    gamma = if (components$concordant + components$discordant > 0) {
      4 * sum(tab * (components$discordant * A -
                       components$concordant * D)^2) /
        (components$concordant + components$discordant)^4
    } else {
      NA_real_
    },
    `tau-b` = (
      sum(tab * (2 * u * d_cell + all.estimates[["tau-b"]] * v_cell)^2) -
        n^3 * all.estimates[["tau-b"]]^2 * (u_row + u_col)^2
    ) / (4 * u^4),
    `tau-c` = 4 * d^2 * var0S / ((d - 1)^2 * n^4),
    `somers-c|r` = sum(
      tab * (u_row * d_cell - components$S *
               outer(n - rowSums(tab), rep(1, ncol(tab))))^2
    ) / u_row^4,
    `somers-r|c` = sum(
      tab * (u_col * d_cell - components$S *
               outer(rep(1, nrow(tab)), n - colSums(tab)))^2
    ) / u_col^4,
    pearson = NA_real_,
    spearman = NA_real_
  )

  pearson <- ordinal_assoc_score_correlation(tab, row.scores, col.scores)
  spearman <- ordinal_assoc_spearman(tab)
  all.estimates[c("pearson", "spearman")] <- c(
    pearson$estimate, spearman$estimate
  )
  all.variances[c("pearson", "spearman")] <- c(
    pearson$variance, spearman$variance
  )

  estimates <- all.estimates[measures]
  variances <- all.variances[measures]

  variances[variances < 0 & variances > -1e-10] <- 0
  if (any(variances < 0, na.rm = TRUE)) {
    stop("A negative sampling variance was computed; check the table and formulas.")
  }
  std.error <- sqrt(variances)
  z_critical <- stats::qnorm(1 - (1 - conf.level) / 2)
  conf.int <- cbind(
    lower = estimates - z_critical * std.error,
    upper = estimates + z_critical * std.error
  )

  null.variances <- stats::setNames(rep(var0S, length(measures)), measures)
  null.variances["pearson"] <- pearson$null.variance
  null.variances["spearman"] <- spearman$null.variance
  numerators <- estimates
  is.correlation <- names(numerators) %in% c("pearson", "spearman")
  numerators[!is.correlation] <- components$S
  statistic <- numerators / sqrt(null.variances[measures])
  unusable.null.variance <- !is.finite(null.variances[measures]) |
    null.variances[measures] <= 0
  statistic[unusable.null.variance] <- NA_real_
  p.value <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(-abs(statistic)),
    greater = stats::pnorm(statistic, lower.tail = FALSE),
    less = stats::pnorm(statistic)
  )

  structure(
    list(
      estimate = estimates,
      std.error = std.error,
      conf.int = conf.int,
      conf.level = conf.level,
      statistic = statistic,
      p.value = p.value,
      alternative = alternative,
      method = "Asymptotic ordinal association measures and tests",
      measure = measures,
      S = components$S,
      var.S = var0S,
      concordant = components$concordant,
      discordant = components$discordant,
      n = n,
      call = match.call()
    ),
    class = "ibist_ordinal_assoc"
  )
}

ordinal_assoc_components <- function(tab) {
  nr <- nrow(tab)
  nc <- ncol(tab)
  concordant.by.cell <- matrix(0, nrow = nr, ncol = nc)
  discordant.by.cell <- matrix(0, nrow = nr, ncol = nc)

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (i > 1L && j > 1L) {
        concordant.by.cell[i, j] <- concordant.by.cell[i, j] +
          sum(tab[seq_len(i - 1L), seq_len(j - 1L), drop = FALSE])
      }
      if (i < nr && j < nc) {
        concordant.by.cell[i, j] <- concordant.by.cell[i, j] +
          sum(tab[seq.int(i + 1L, nr), seq.int(j + 1L, nc), drop = FALSE])
      }
      if (i > 1L && j < nc) {
        discordant.by.cell[i, j] <- discordant.by.cell[i, j] +
          sum(tab[seq_len(i - 1L), seq.int(j + 1L, nc), drop = FALSE])
      }
      if (i < nr && j > 1L) {
        discordant.by.cell[i, j] <- discordant.by.cell[i, j] +
          sum(tab[seq.int(i + 1L, nr), seq_len(j - 1L), drop = FALSE])
      }
    }
  }

  concordant <- sum(tab * concordant.by.cell) / 2
  discordant <- sum(tab * discordant.by.cell) / 2
  S <- concordant - discordant
  n <- sum(tab)
  var.S <- sum(tab * (concordant.by.cell - discordant.by.cell)^2) -
    4 * S^2 / n

  structure(
    list(
      concordant = concordant,
      discordant = discordant,
      S = S,
      var.S = var.S,
      concordant.by.cell = concordant.by.cell,
      discordant.by.cell = discordant.by.cell
    )
  )
}

# Helpers for score-based correlation measures.
validate_ordinal_scores <- function(x, n, arg) {
  if (!is.numeric(x) || length(x) != n || anyNA(x) ||
      any(!is.finite(x)) || length(unique(x)) < 2L) {
    stop("'", arg, "' must contain ", n,
         " finite numeric scores with at least two distinct values.")
  }
  as.numeric(x)
}

default_ordinal_scores <- function(labels, n) {
  if (is.null(labels)) {
    return(seq_len(n))
  }
  scores <- suppressWarnings(as.numeric(labels))
  if (all(!is.na(scores)) && all(is.finite(scores))) scores else seq_len(n)
}

ordinal_assoc_score_correlation <- function(tab, row.scores, col.scores) {
  n <- sum(tab)
  row.margin <- rowSums(tab)
  col.margin <- colSums(tab)
  row.mean <- sum(row.margin * row.scores) / n
  col.mean <- sum(col.margin * col.scores) / n
  row.centered <- row.scores - row.mean
  col.centered <- col.scores - col.mean
  row.ss <- sum(row.margin * row.centered^2)
  col.ss <- sum(col.margin * col.centered^2)
  cross.ss <- sum(tab * outer(row.centered, col.centered))
  scale <- sqrt(row.ss * col.ss)
  if (!is.finite(scale) || scale <= 0) {
    stop("The category scores must have positive marginal variance.")
  }
  estimate <- cross.ss / scale
  bij <- outer(row.centered^2, rep(1, ncol(tab))) * col.ss +
    outer(rep(1, nrow(tab)), col.centered^2) * row.ss
  influence <- scale * outer(row.centered, col.centered) -
    bij * cross.ss / (2 * scale)
  variance <- sum(tab * influence^2) / scale^4
  null.variance <- (
    sum(tab * outer(row.centered^2, col.centered^2)) - cross.ss^2 / n
  ) / (row.ss * col.ss)

  list(
    estimate = estimate,
    variance = variance,
    null.variance = null.variance
  )
}

ordinal_assoc_spearman <- function(tab) {
  n <- sum(tab)
  row.margin <- rowSums(tab)
  col.margin <- colSums(tab)
  row.rank <- cumsum(row.margin) - row.margin / 2
  col.rank <- cumsum(col.margin) - col.margin / 2
  row.centered <- row.rank - n / 2
  col.centered <- col.rank - n / 2
  F <- n^3 - sum(row.margin^3)
  G <- n^3 - sum(col.margin^3)
  scale <- sqrt(F * G) / 12
  if (!is.finite(scale) || scale <= 0) {
    stop("Spearman correlation requires at least two non-empty categories.")
  }
  v <- sum(tab * outer(row.centered, col.centered))
  estimate <- v / scale

  vij <- matrix(0, nrow(tab), ncol(tab))
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      row.tail <- if (i < nrow(tab)) {
        sum(tab[seq.int(i + 1L, nrow(tab)), , drop = FALSE] %*% col.centered)
      } else {
        0
      }
      col.tail <- if (j < ncol(tab)) {
        sum(row.centered %*% tab[, seq.int(j + 1L, ncol(tab)), drop = FALSE])
      } else {
        0
      }
      vij[i, j] <- n * (
        row.centered[i] * col.centered[j] +
          sum(tab[i, ] * col.centered) / 2 +
          sum(tab[, j] * row.centered) / 2 + row.tail + col.tail
      )
    }
  }
  wij <- -n / (96 * scale) * (
    outer(row.margin^2, rep(G, ncol(tab))) +
      outer(rep(F, nrow(tab)), col.margin^2)
  )
  zij <- scale * vij - v * wij
  z.mean <- sum(tab * zij) / n
  variance <- sum(tab * (zij - z.mean)^2) / (n^2 * scale^4)
  v.mean <- sum(tab * vij) / n
  null.variance <- sum(tab * (vij - v.mean)^2) / (n^2 * scale^2)

  list(
    estimate = estimate,
    variance = variance,
    null.variance = null.variance
  )
}

#' @export
print.ibist_ordinal_assoc <- function(
  x,
  digits = 7L,
  ...
) {
  cat("\nAsymptotic ordinal association measures\n")
  cat("Confidence level:", format(100 * x$conf.level), "%\n\n")

  result <- cbind(
    estimate = x$estimate,
    `std. error` = x$std.error,
    lower = x$conf.int[, "lower"],
    upper = x$conf.int[, "upper"]
  )
  print(round(result, digits = digits))

  cat("\nTests of zero association\n")
  concordance.measures <- intersect(
    x$measure,
    c("gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c")
  )
  if (length(concordance.measures) > 1L) {
    cat(
      "  Note: Concordance-based measures share the same test; their ",
      "identical test results are not separate tests.\n",
      sep = ""
    )
  }
  cat("  Alternative:", x$alternative, "\n")
  tests <- cbind(
    statistic = x$statistic,
    p.value = x$p.value
  )
  print(tests, digits = digits)

  invisible(x)
}
