## Chan--Zhang tail inversion. All calculations are independent of Exact.
## Search intervals are discarded only when their entire acceptance set
## is excluded, rather than assuming that the exact p-value is monotone.
independent_exact_limits <- function(x, n, tail, alternative, tol,
                                     nuisance.tol) {
  if (prod(n + 1) > 1e6) {
    stop("Chan-Zhang enumeration requires at most 1,000,000 possible tables.")
  }
  lower <- function(x, n) {
    estimate <- x[1L] / n[1L] - x[2L] / n[2L]
    if (estimate == -1) return(-1)
    support <- cbind(rep(0:n[1L], n[2L] + 1),
                     rep(0:n[2L], each = n[1L] + 1))
    observed <- x[1L] + 1 + x[2L] * (n[1L] + 1)
    differences <- support[, 1L] / n[1L] - support[, 2L] / n[2L]
    probabilities <- function(delta) {
      left <- rep(max(0, -delta), nrow(support))
      right <- rep(min(1, 1 - delta), nrow(support))
      ## Concavity of the constrained log likelihood brackets every MLE.
      for (iteration in seq_len(52L)) {
        q <- (left + right) / 2
        p <- cbind(q + delta, q)
        p <- pmax(pmin(p, 1 - .Machine$double.eps),
                  .Machine$double.xmin)
        slope <- rowSums(support / p - (rep(n, each = nrow(support)) -
                                         support) / (1 - p))
        positive <- slope > 0
        left[positive] <- q[positive]
        right[!positive] <- q[!positive]
      }
      list(low = cbind(left + delta, left),
           high = cbind(right + delta, right))
    }
    score.bounds <- function(left, right, pl, pr) {
      ## Each constrained p1 increases with delta, while p2 decreases.
      ## Enclose p(1-p) over those ranges before dividing the numerator.
      low <- cbind(pl$low[, 1L], pr$low[, 2L])
      high <- cbind(pr$high[, 1L], pl$high[, 2L])
      low <- pmax(pmin(low, 1), 0)
      high <- pmax(pmin(high, 1), 0)
      vlo <- pmin(low * (1 - low), high * (1 - high))
      vhi <- pmax(low * (1 - low), high * (1 - high))
      vhi[low <= 0.5 & high >= 0.5] <- 0.25
      divisor <- rep(n, each = nrow(support))
      vlo <- sqrt(rowSums(vlo / divisor))
      vhi <- sqrt(rowSums(vhi / divisor))
      divide <- function(a, b) {
        value <- a / b
        value[a == 0 & b == 0] <- 0
        value
      }
      numerator.low <- differences - right
      numerator.high <- differences - left
      low.score <- divide(numerator.low,
                           ifelse(numerator.low >= 0, vhi, vlo))
      high.score <- divide(numerator.high,
                            ifelse(numerator.high >= 0, vlo, vhi))
      ## At the two constant-response corners both numerator and variance
      ## vanish at delta=0. Their analytic limits avoid an unbounded 0/0
      ## enclosure that would otherwise retain zero for every interval.
      corner <- function(delta, sizes) {
        if (delta == 0) return(0)
        if (delta < 0) sqrt(sizes[2L] * (-delta) / (1 + delta)) else {
          -sqrt(sizes[1L] * delta / (1 - delta))
        }
      }
      for (index in c(1L, nrow(support))) {
        sizes <- if (index == 1L) n else rev(n)
        low.score[index] <- corner(right, sizes)
        high.score[index] <- corner(left, sizes)
      }
      ## Include numerical ties conservatively.
      pad <- 256 * .Machine$double.eps * pmax(1, abs(low.score[observed]))
      if (!is.finite(pad)) pad <- 0
      included <- matrix(high.score >= low.score[observed] - pad,
                          nrow = n[1L] + 1)
      ## Enclose the union by an event increasing in x1 and decreasing in x2.
      cutoff <- apply(included, 1L, function(row) {
        positions <- which(row)
        if (length(positions)) max(positions) - 1L else -1L
      })
      cummax(cutoff)
    }
    possible <- function(left, right, pl, pr) {
      cutoff <- score.bounds(left, right, pl, pr)
      range <- c(max(0, -right), min(1, 1 - left))
      ## Increasing delta increases the probability of this enclosing event.
      ## Above q=1-right, p1 is capped at one; split at this kink.
      breaks <- sort(unique(c(range, max(range[1L], min(range[2L],
                                                       1 - right)))))
      probability <- function(q) {
        p1 <- pmax(0, pmin(1, q + right))
        cdf <- stats::pbinom(cutoff, n[2L], q)
        mass <- stats::dbinom(0:n[1L], n[1L], p1)
        value <- sum(mass * cdf)
        slope <- -n[2L] * sum(mass *
                                stats::dbinom(cutoff, n[2L] - 1, q))
        if (q + right < 1) {
          slope <- slope + n[1L] * sum(
            stats::dbinom(0:(n[1L] - 1), n[1L] - 1, p1) * diff(cdf)
          )
        }
        c(value, slope)
      }
      if (any(vapply(breaks, function(q) probability(q)[1L], 0.0) > tail)) {
        return(TRUE)
      }
      l <- breaks[-length(breaks)]
      r <- breaks[-1L]
      evaluations <- 0L
      repeat {
        mid <- (l + r) / 2
        values <- vapply(mid, probability, c(0.0, 0.0))
        if (any(values[1L, ] > tail)) return(TRUE)
        h <- (r - l) / 2
        ## For a Boolean event in N independent Bernoulli trials,
        ## |f''(q)| <= 2*N*(N-1). Taylor's bound covers every local maximum.
        upper <- pmin(1, values[1L, ] + abs(values[2L, ]) * h +
                        sum(n) * (sum(n) - 1) * h^2 +
                        256 * .Machine$double.eps)
        if (all(upper <= tail)) return(FALSE)
        if (max(upper) - max(values[1L, ]) <= nuisance.tol) return(TRUE)
        keep <- upper > tail
        l <- c(l[keep], mid[keep])
        r <- c(mid[keep], r[keep])
        evaluations <- evaluations + length(mid)
        if (evaluations > 1e6) {
          stop("Nuisance maximization did not converge; ",
               "increase 'nuisance.tol'.")
        }
      }
    }
    ## Explore the full parameter space from its lower boundary upward.
    ## Cached endpoint MLE brackets avoid repeating the constrained fits.
    search <- function(left, right, pl, pr) {
      if (!possible(left, right, pl, pr)) return(NULL)
      if (right - left <= tol) return(left)
      middle <- (left + right) / 2
      if (middle == left || middle == right) return(left)
      pm <- probabilities(middle)
      answer <- search(left, middle, pl, pm)
      if (!is.null(answer)) return(answer)
      search(middle, right, pm, pr)
    }
    answer <- search(-1, 1, probabilities(-1), probabilities(1))
    if (is.null(answer)) 1 else answer
  }
  lo <- if (alternative == "less") -1 else lower(x, n)
  hi <- if (alternative == "greater") 1 else -lower(rev(x), rev(n))
  c(lo, hi)
}
