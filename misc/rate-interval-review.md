# Review: negative Wilson--Hilferty rate limit

Status: pending author review. The R WH formula has not been changed.
Observed while preparing the Python rate-function port on 2026-11-07.
R source revision: `5534927c6bd27bc3f949a85aba179ebe49682105`.

## Reproduction

From the R package root:

```r
source("R/ci.R")
source("R/rate_1s_ci.R")
rate.1s.ci(0, exposure = 10, conf.level = 0.95,
           method = "wh", correct = TRUE)$conf.int
```

The limits are:

```text
lower = -0.00015611287295404339
upper =  0.24639362099673709
```

With `correct = FALSE`, the same data give:

```text
lower = 0
upper = 0.36680118321737509
```

## Cause

`ci_wh()` in `R/rate_1s_ci.R` uses `x + 0.5` for both endpoints when
`correct = TRUE`. At zero events the lower-limit expression is

```text
0.5 * (1 - 1 / (9 * 0.5) - qnorm(0.975) / (3 * sqrt(0.5)))^3 / exposure
```

The quantity being cubed is negative. The current guard checks only whether
the adjusted count is nonpositive; the adjusted count is 0.5 in this case,
so the negative limit is returned. A Poisson mean or rate cannot be negative.

## Proposed handling

Restrict the computed lower limit to zero with `max(0, lower)`, retaining
the current formula and upper limit. For the reproduction above, the proposed
interval is `[0, 0.24639362099673709]`.

This restricts the interval to the parameter space; it does not establish
coverage or validate the midpoint correction itself. Reviewing that
correction against the book's intended construction is a separate decision.

The Python port will apply this lower-limit restriction and document it as
an intentional difference from the current R output. Any later R change
should add a zero-count regression test and regenerate comparison cases.

## Scope

R rate functions now accept the common `side` argument; this separate API
change preserves the two-sided WH formula above. Book sources remain
unchanged. Review zero counts across confidence levels and exposure scales
before changing the R WH formula.
