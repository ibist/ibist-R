#' Pearson Chi-squared Goodness-of-Fit Test for Logistic Regression
#'
#' Performs the Pearson chi-squared goodness-of-fit (GoF) test for a
#' logistic regression model fitted via \code{glm(..., family = "binomial")}.
#' The test assesses calibration by comparing observed and expected
#' counts across covariate patterns.
#'
#' The Pearson GoF statistic is defined as
#' \deqn{
#'   T_P = \sum_{\ell=1}^J \frac{(O_\ell - E_\ell)^2}{V_\ell},
#' }
#' where \eqn{O_\ell} is the observed number of events, \eqn{E_\ell =
#' n_\ell \hat{\pi}_\ell} is the expected number of events, and
#' \eqn{V_\ell = n_\ell \hat{\pi}_\ell (1 - \hat{\pi}_\ell)} is the
#' variance for the \eqn{\ell}-th covariate pattern. Under the null
#' hypothesis of correct model specification, \eqn{T_P} approximately
#' follows a chi-squared distribution with \eqn{J - (k + 1)} degrees of
#' freedom, where \eqn{J} is the number of groups and \eqn{k + 1} is the
#' number of model parameters.
#'
#' This test requires replicated covariate patterns or grouped binomial
#' data. When groups are small, optional pooling can be applied to ensure
#' stable chi-squared approximation.
#'
#' @param fit A fitted \code{glm} object with \code{family = "binomial"}.
#' @param min_n Minimum total count per group for pooling. Default is 5.
#' @param min_expected Minimum expected count per group for pooling.
#'   Default is 5.
#' @param pool Logical; if \code{TRUE}, adjacent groups (ordered by fitted
#'   probabilities) are pooled to meet minimum thresholds. Default is
#'   \code{TRUE}.
#'
#' @return An object of class \code{"htest"} with components:
#' \item{statistic}{The Pearson chi-squared statistic.}
#' \item{parameter}{Degrees of freedom.}
#' \item{p.value}{p-value of the test.}
#' \item{method}{Description of the test.}
#' \item{data.name}{Description of the data.}
#'
#' Additional components (not printed) include:
#' \item{observed}{Observed counts by group.}
#' \item{expected}{Expected counts by group.}
#' \item{group_n}{Group sizes.}
#' \item{groups}{Data frame of grouped summaries.}
#'
#' @details
#' The function supports three types of binomial responses:
#' \itemize{
#'   \item Binary outcomes (0/1)
#'   \item Proportions with weights (grouped binomial)
#'   \item Two-column matrix response via \code{cbind(success, failure)}
#' }
#'
#' When no replication exists (i.e., each observation has a unique
#' covariate pattern), the Pearson GoF test is not appropriate. In such
#' cases, consider alternative methods such as the Hosmer--Lemeshow test.
#'
#' Pooling stabilizes the test by merging adjacent groups with small
#' counts, trading a small bias for reduced variance in the chi-squared
#' approximation.
#'
#' @examples
#' ## Example: Surfactant use and birthweight (RDS data)
#' data(rds, package = "ibist")
#'
#' fit <- glm(death ~ surf + bwt,
#'   data = rds, weights = count,
#'   family = "binomial"
#' )
#'
#' ## Pearson GoF test with pooling (default)
#' (res <- logistic_gof(fit))
#'
#' ## Without pooling (may be unstable if small groups exist)
#' logistic_gof(fit, pool = FALSE)
#'
#' ## Inspect grouped diagnostics
#' res$groups
#'
#' @seealso \code{\link{glm}}, \code{\link{chisq.test}}
#'
#' @importFrom stats family
#' @importFrom stats model.frame model.response model.weights
#' @importFrom stats fitted model.matrix aggregate coef pchisq formula
#' @export
logistic_gof <- function(fit, min_n = 5, min_expected = 5,
                         pool = TRUE) {
  ## --- Sanity checks ---
  if (!inherits(fit, "glm")) {
    stop("fit must be a glm object.")
  }

  if (family(fit)$family != "binomial") {
    stop("Model must be binomial.")
  }

  if (!is.logical(pool) || length(pool) != 1L || is.na(pool)) {
    stop("'pool' must be TRUE or FALSE.")
  }

  if (
    !is.numeric(min_n) || length(min_n) != 1L ||
      !is.finite(min_n) || min_n <= 0
  ) {
    stop("'min_n' must be a single positive number.")
  }

  if (
    !is.numeric(min_expected) || length(min_expected) != 1L ||
      !is.finite(min_expected) || min_expected < 0
  ) {
    stop("'min_expected' must be a single non-negative number.")
  }

  mf <- model.frame(fit)
  y <- model.response(mf)
  w <- model.weights(mf)

  response <- binomial_response_counts(y, w)
  observed_i <- response$observed
  n_i <- response$n

  pi_hat <- fitted(fit)

  if (any(pi_hat <= 0 | pi_hat >= 1)) {
    warning("Fitted probabilities at boundary; variance may be unstable.")
  }

  ## --- Construct groups ---
  model_mat <- model.matrix(fit)
  key <- apply(model_mat, 1, paste, collapse = ":")

  df <- data.frame(key = key, O = observed_i, n = n_i, pi_hat = pi_hat)

  agg <- aggregate(
    cbind(O = O, n = n, piw = pi_hat * n) ~ key,
    data = df,
    FUN = sum
  )

  observed <- agg$O
  n <- agg$n
  pi_g <- agg$piw / n

  df_group <- data.frame(O = observed, n = n, pi_hat = pi_g)

  n_patterns <- nrow(df_group)
  k <- fit$rank
  if (is.null(k)) {
    k <- sum(!is.na(coef(fit)))
  }

  if (n_patterns <= k) {
    stop("Not enough unique covariate patterns: J <= k.")
  }

  if (n_patterns == NROW(y)) {
    warning("No replication detected; Pearson GoF not appropriate.")
  }

  ## --- Pool if requested ---
  if (pool) {
    df_group <- pool_groups(df_group,
      min_n = min_n,
      min_expected = min_expected
    )
  }

  ## --- Check post-pooling ---
  n_groups <- nrow(df_group)
  df_chi <- n_groups - k

  if (df_chi <= 0) {
    stop("Degrees of freedom <= 0 after pooling.")
  }

  ## --- Compute statistic ---
  expected <- df_group$n * df_group$pi_hat
  variance <- df_group$n * df_group$pi_hat * (1 - df_group$pi_hat)

  if (any(variance == 0)) {
    stop("Zero variance encountered; adjust pooling thresholds.")
  }

  statistic <- sum((df_group$O - expected)^2 / variance)
  pval <- 1 - pchisq(statistic, df_chi)

  ## --- Return htest ---
  res <- list(
    statistic = c(X2 = statistic),
    parameter = c(df = df_chi),
    p.value = pval,
    method = paste(
      "Pearson Chi-squared GoF for logistic regression",
      if (pool) "(with pooling)" else ""
    ),
    data.name = deparse(formula(fit))
  )

  # attach diagnostics
  res$observed <- df_group$O
  res$expected <- expected
  res$group_n <- df_group$n
  res$groups <- df_group

  class(res) <- "htest"
  res
}

binomial_response_counts <- function(y, w = NULL) {
  if (is.matrix(y)) {
    if (ncol(y) != 2L) {
      stop("Matrix response must have two columns: successes and failures.")
    }
    if (any(!is.finite(y)) || any(y < 0)) {
      stop("Matrix response counts must be finite and non-negative.")
    }
    n <- rowSums(y)
    if (any(n <= 0)) {
      stop("Matrix response rows must have positive totals.")
    }
    if (is.null(w)) {
      w <- rep(1, nrow(y))
    } else if (any(!is.finite(w)) || any(w <= 0)) {
      stop("weights must be positive.")
    }
    return(list(observed = y[, 1] * w, n = n * w))
  }

  if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0 | y > 1)) {
    stop(
      paste(
        "Response must be binary, a proportion with weights,",
        "or a two-column matrix."
      )
    )
  }

  has_weights <- !is.null(w)
  if (is.null(w)) w <- rep(1, length(y))
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("weights must be positive.")
  }

  if (!all(y %in% c(0, 1)) && !has_weights) {
    stop("Proportion responses require binomial weights.")
  }

  list(observed = y * w, n = w)
}

pool_groups <- function(df, min_n = 5, min_expected = 5) {
  if (!all(c("O", "n", "pi_hat") %in% names(df))) {
    stop("df must contain O, n, pi_hat.")
  }

  # order by fitted probability
  df <- df[order(df$pi_hat), ]

  pooled <- list()
  current_observed <- 0
  cur_n <- 0
  cur_piw <- 0

  for (i in seq_len(nrow(df))) {
    current_observed <- current_observed + df$O[i]
    cur_n <- cur_n + df$n[i]
    cur_piw <- cur_piw + df$n[i] * df$pi_hat[i]

    pi_curr <- cur_piw / cur_n
    current_expected <- cur_n * pi_curr

    if (cur_n >= min_n && current_expected >= min_expected) {
      pooled[[length(pooled) + 1]] <-
        data.frame(O = current_observed, n = cur_n, pi_hat = pi_curr)
      current_observed <- 0
      cur_n <- 0
      cur_piw <- 0
    }
  }

  # Merge leftover rows into the last group.
  if (cur_n > 0) {
    if (length(pooled) == 0) {
      pooled[[1]] <- data.frame(
        O = current_observed,
        n = cur_n,
        pi_hat = cur_piw / cur_n
      )
    } else {
      last <- pooled[[length(pooled)]]
      new_n <- last$n + cur_n
      new_pi <- (last$n * last$pi_hat + cur_piw) / new_n

      pooled[[length(pooled)]] <- data.frame(
        O = last$O + current_observed,
        n = new_n,
        pi_hat = new_pi
      )
    }
  }

  do.call(rbind, pooled)
}
