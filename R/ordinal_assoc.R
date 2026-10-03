#' Ordinal Association Measures with Asymptotic Confidence Intervals and Test
#'
#' Computes concordance-based ordinal association measures and their
#' asymptotic confidence intervals from an RxC table. It also tests the common
#' null hypothesis that the association measure is zero.
#'
#' @param x A two-dimensional matrix or table of non-negative integer counts.
#'   Rows and columns must be ordered from the lowest to highest category.
#' @param measure One or more measures: \code{"gamma"}, \code{"tau-b"},
#'   \code{"tau-c"}, \code{"somers-c|r"}, and \code{"somers-r|c"}. By
#'   default, all five are computed.
#' @param conf.level Confidence level for the two-sided asymptotic confidence
#'   intervals.
#' @param alternative Alternative hypothesis for the shared test:
#'   \code{"two.sided"}, \code{"greater"}, or \code{"less"}. Positive
#'   association means that higher row categories tend to occur with higher
#'   column categories.
#'
#' @details
#' The procedure computes the concordant and discordant pair counts,
#' \eqn{n_C} and \eqn{n_D}, and sets \eqn{S=n_C-n_D}. Its null variance is
#' \deqn{\mathrm{Var}_0(S) = \sum_{ij} n_{ij}(A_{ij}-D_{ij})^2 - 4S^2/n,}
#' where \eqn{A_{ij}} and \eqn{D_{ij}} are the concordant and discordant
#' contributions for cell \eqn{(i,j)}. The test statistic is
#' \eqn{Z=S/\sqrt{\mathrm{Var}_0(S)}}.
#'
#' The selected measures have different sampling variances, so their standard
#' errors and confidence intervals differ. Their tests of a zero association
#' nevertheless reduce to the same \eqn{Z} statistic and p-value. Confidence
#' intervals are two-sided regardless of the test alternative. Rows and columns
#' with zero marginal totals are dropped before computation.
#'
#' The \code{"somers-c|r"} measure is Somers' \eqn{\Delta(C|R)}, treating the
#' row variable as explanatory; \code{"somers-r|c"} is
#' \eqn{\Delta(R|C)}, treating the column variable as explanatory.
#'
#' @return An object of class \code{"ibist_ordinal_assoc"} containing:
#' \describe{
#'   \item{estimate}{Named vector of selected association estimates.}
#'   \item{std.error}{Named vector of asymptotic standard errors.}
#'   \item{conf.int}{Matrix of lower and upper confidence limits.}
#'   \item{statistic}{The shared asymptotic \eqn{Z} statistic.}
#'   \item{p.value}{P-value for the selected alternative.}
#'   \item{S, var.S}{The observed concordance contrast and its null variance.}
#'   \item{concordant, discordant}{Concordant and discordant pair counts.}
#' }
#'
#' @references
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
#' ordinal.assoc(tab, measure = c("tau-b", "somers-c|r"))
#'
#' @export
ordinal.assoc <- function(
  x,
  measure = c("gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c"),
  conf.level = 0.95,
  alternative = c("two.sided", "greater", "less")
) {
  measures <- if (missing(measure)) {
    c("gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c")
  } else {
    match.arg(
      measure,
      choices = c("gamma", "tau-b", "tau-c", "somers-c|r", "somers-r|c"),
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

  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (any(dim(tab) < 2L)) {
    stop("'x' must have at least two non-empty rows and columns.")
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
    `somers-r|c` = components$S / u_col
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
    ) / u_col^4
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

  statistic <- if (var0S > 0) {
    components$S / sqrt(var0S)
  } else {
    NA_real_
  }
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
      statistic = c(Z = statistic),
      p.value = p.value,
      alternative = alternative,
      method = "Asymptotic ordinal association measures and test",
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

  cat("\nShared test of zero association\n")
  cat("  Alternative:", x$alternative, "\n")
  cat("  S =", format(x$S),
      "; Var0(S) =", format(x$var.S),
      "; Z =", format(unname(x$statistic)), "\n")
  cat("  p-value =", format.pval(x$p.value, digits = digits), "\n")

  invisible(x)
}
