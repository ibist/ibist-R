# Binomial interval specification

Agreed naming and interface principles: 2026-11-06. R is implemented first;
the corresponding Python functions remain to be implemented.

## Interfaces and results

| Estimand | R | Python | Default method |
| --- | --- | --- | --- |
| One probability | `prop.1s.ci(x, n, ...)` | `prop_1s_ci(x, n, ...)` | `wilson` |
| Difference of independent probabilities | `prop.2s.ci(x, n, ...)` | `prop_2s_ci(x, n, ...)` | `newcombe` |
| Risk ratio | `rr.2x2.ci(x, ...)` | `rr_2x2_ci(x, ...)` | `koopman` |

All default to 95% two-sided intervals. Explicit method vectors request
several methods. One-sided alternatives `greater` and `less` request a lower
or upper bound, respectively, with the other limit at the parameter-space
boundary. Defaults are `correct = FALSE` and `midp = FALSE` where supported.
Python argument names will use snake_case.

One-sample counts may be vectors; a scalar trial count is recycled, otherwise
the trial-count vector must have the same length. Two-sample counts and trial
counts are vectors of length two. The difference function also accepts a
2x2 table instead of those vectors. The risk-ratio function accepts a 2x2
table. In both cases rows are independent samples, the first column contains
successes/events, and the second contains failures/non-events. The difference
is row 1 minus row 2; the ratio is row 1 divided by row 2.

R results use the existing `ibist_ci` class. `estimate` is the named observed
estimate, separate from the `conf.int` matrix, which has `lower` and `upper`
columns. Other fields include `conf.level`, `method`, and `alternative`.
One-sample vector results have one row per input and method, ordered by input
and then method; `inputs` explicitly identifies the corresponding counts and
trial counts. Difference and ratio results include `sample.proportions`.

Method labels identify corrections, for example `wilson-cc` and
`clopper-pearson-midp`. Uncorrected and corrected analyses use separate calls.
Printing displays a single pair of limits or a method comparison table,
with base R-style data and sample-estimate fields. `digits` controls printed
significant digits, not stored precision. There is no fabricated hypothesis
test, statistic, or p-value in an interval-only result.

## Statistical definitions

For one sample, write `p = x/n` and let `z` be the standard-normal quantile
at `1 - alpha/2`, or `1 - alpha` for a one-sided bound.

- **Wald:** `p +/- z * sqrt(p*(1-p)/n)`. Correction adds `1/(2n)` to the
  margin. Restrict the endpoints to `[0,1]`.
- **Wilson:** invert the score equation. For correction, replace `p` by
  `(x-0.5)/n` for the lower limit and `(x+0.5)/n` for the upper limit,
  clipping those adjusted proportions to `[0,1]`. The half-count is not
  capped using a null probability. The boundary limits are zero for `x=0`
  and one for `x=n`.
- **Agresti–Coull:** use `n+z^2` and `x+z^2/2` in the Wald formula, then
  restrict to `[0,1]`. The reported estimate remains `x/n`.
- **Clopper–Pearson:** invert the inclusive binomial tails using beta
  quantiles. Mid-p subtracts half the probability of the observed count
  from each inclusive tail before inversion. Mid-p is not coverage-exact.
- **Blaker:** take the interval hull of the two-sided acceptance set based
  on smaller inclusive-tail ordering. Partition the search at opposite-tail
  ordering changes, so discontinuities and possible gaps are respected.
  Use exact binomial tail inversion for a one-sided bound.

For independent proportions:

- **Newcombe:** combine the distances from each observed proportion to its
  Wilson limits in quadrature. Continuity correction uses corrected Wilson
  limits. This differs from inverting the pooled score test for equality.
- **Wald:** use the unpooled binomial variance. Correction adds
  `1/(2n1)+1/(2n2)` to the margin.
- **Agresti–Caffo:** add one success and one failure to each sample and
  compute the Wald interval using those adjusted samples.
- **Miettinen–Nurminen:** maximize the independent binomial likelihood under
  `p1-p2=delta`, then invert the score statistic using the constrained
  variance multiplied by `(n1+n2)/(n1+n2-1)`. Always include this variance
  factor; it is distinct from continuity correction.
- **Mee:** use the same constrained likelihood estimates and score inversion
  as Miettinen–Nurminen, with no `(n1+n2)/(n1+n2-1)` variance factor. Select
  `method = "mee"`; `correct` remains a continuity-correction argument for
  Wald and Newcombe and does not switch between these score methods.
- **Brown–Li:** set `pbar_i=(x_i+0.5)/(n_i+1)` and form a normal-approximation
  interval centered at `pbar_1-pbar_2`, with variance
  `pbar_1*(1-pbar_1)/n1 + pbar_2*(1-pbar_2)/n2`. The denominators are the
  original sample sizes. This is not a posterior-quantile interval; retain
  the observed difference as the reported estimate. Continuity correction
  does not apply. Restrict the limits to `[-1,1]`.

- **Chan–Zhang:** select `method = "chan-zhang"` to invert inclusive
  one-sided unconditional score tests at `alpha/2` for two-sided intervals,
  or `alpha` for a one-sided bound. Enumerate both binomial counts and
  maximize the tail probability over the feasible nuisance probability.
  The common Miettinen–Nurminen variance factor cancels in the ordering.
  This is the tail/central construction, not the Agresti–Min square interval.
  Neither continuity correction nor mid-p applies; the default stays Newcombe.
  Constrained-MLE brackets enclose the score ordering throughout each interval
  of candidate differences. Writing `p2=q`, implicit differentiation of the
  constrained likelihood equation gives `dq/d(delta)=-I1/(I1+I2)` for an
  interior maximum, where each `Ii` is the nonnegative observed binomial
  information. Hence fitted `p1` increases and fitted `p2` decreases with
  `delta`; boundary fits have the same monotonicity. These ranges bound each
  variance and score. The outer search discards only wholly rejected
  intervals, so discontinuities and disconnected accepted components are
  retained. A nuisance search uses first derivatives and the bound
  `abs(f''(q)) <= 2*N*(N-1)`, where `N=n1+n2`, to cover all local maxima.
  This follows by writing the event probability as an expectation of a Boolean
  function of `N` Bernoulli variables: every mixed second finite difference
  has absolute value at most two. Split the search at a capped probability's
  kink. The endpoint tolerance is `tol`; absolute probability tolerance
  `nuisance.tol` defaults to `1e-10`. Unresolved differences within that
  probability tolerance are accepted; endpoints are rounded outward.
  These are conservative floating-point enclosures, not symbolic arithmetic.
  Computation is intended for small samples and errors if
  `(n1+1)*(n2+1) > 1e6`. Large samples should use approximate intervals.

Restrict difference limits to `[-1,1]`. For the ratio:

- **Katz:** exponentiate a Wald interval for `log(p1/p2)`, with variance
  `1/x1-1/n1+1/x2-1/n2`. Require positive event counts; do not add pseudo-counts.
- **Koopman:** maximize the binomial likelihood under `p1=r*p2` and invert
  the constrained Pearson statistic. No finite-sample factor or continuity
  correction is used. Solve on the log-ratio scale, swapping rows and taking
  reciprocals to preserve symmetry and avoid large-ratio arithmetic.

Koopman permits a zero event count, giving a boundary estimate and limit.
With zero events in both rows, report an undefined estimate (`NA` in R) and
`[0,Inf]`. Both sample sizes must still be positive.

Flags that affect none of the requested methods are errors. In a mixed
request, correction affects only Wald/Wilson or Wald/Newcombe, and mid-p
affects only Clopper–Pearson. Unsupported combinations must not be silently
renamed as other statistical methods.

## Validation and migration

The implementation is independently written from the statistical definitions;
DescTools is used for comparison only. It is not an ibist runtime dependency.

- Published SAS/STAT 14.2 FREQ Output 41.4.1 (p. 2857) supplies one-proportion
  Agresti–Coull, Wilson, and Clopper–Pearson limits for 341/762 at 90%.
- Output 41.5.3 (p. 2861) supplies Katz risk-ratio limits for rows `(11,4)`
  and `(2,6)`, for both choices of event column.
- Brown–Li reference values from the search-extracted WUSS table are
  `(56/70,48/80)`, with 95% limits
  `(0.0540,0.3400)`, and `(9/10,3/10)`, with limits `(0.1869,0.9040)`.
  For `(10/10,0/20)`, the published lower limit is `0.7854` and the raw upper
  limit exceeds one; ibist truncates it to one. See Table 2, p. 9, in the
  [WUSS comparison paper](https://www.lexjansen.com/wuss/2016/127_Final_Paper_PDF.pdf).
  The PDF has not been downloaded or visually checked, so these values still
  need confirmation against the original document. They supplement the
  independently executed formula and DescTools checks.
  Formula tests distinguish the original sample-size denominators from
  adjusted sample sizes and posterior variances. Hidden DescTools comparisons
  use `method = "blj"` with absolute tolerance `1e-10`.
- Chan–Zhang matches the same search-extracted WUSS Table 2 in all three
  scenarios: `(0.0382,0.3453)`, `(0.1608,0.8795)`, and `(0.6915,1)`.
  The original-PDF qualification above applies to these values too.
  Independent coverage enumeration uses every outcome for sample sizes
  `(1,1)`, `(2,3)`, and `(3,2)`, at 80% and 95%, on an 11-by-11 probability
  grid including boundaries. Reversal, one-sided equivalence, constant
  responses, and a disconnected acceptance component are tested.
  The latter uses the data in Hawila and Berg (2023), Figure 1; direct
  likelihood optimization and full table summation verify the component.
  Direct comparisons with `Exact::exact.test(method = "z-pooled",
  tsmethod = "central")` agree within `1e-4` in the `(9/10,3/10)` and
  `(7/84,2/92)` examples, allowing for Exact's numerical inversion.
  Measured native R times in this workspace were about 1 second for
  `(9/10,3/10)`, 10 seconds for `(56/70,48/80)`, and 35 seconds for
  `(7/84,2/92)` at default tolerances. These are illustrative timings,
  not performance guarantees. Book comparisons retain Exact in hidden code.
- Base R checks cover Wilson and Clopper–Pearson intervals and Wald
  difference intervals. Account explicitly for `prop.test()`'s null-dependent
  cap on continuity correction when comparing corrected intervals.
- Independent likelihood optimization verifies the defining score equations
  at Miettinen–Nurminen, Mee, and Koopman endpoints. Closed-form zero-event
  cases distinguish the Mee and Miettinen–Nurminen variance factors.
  Small-sample enumeration checks
  Blaker acceptance sets and nominal coverage. Tests include zero counts,
  all successes, sample reversal, one-sided bounds, and vector/method ordering.
- Book examples retain executable `include=FALSE` DescTools comparisons.
  Use absolute tolerances: `1e-10` for closed-form calculations, `1e-6` for
  Miettinen–Nurminen, `2e-5` for Blaker/Koopman, and `1e-4` for Mee and mid-p
  at the current DescTools defaults. An 81-case Mee comparison across three
  confidence levels, both one-sided directions, and interior/boundary counts
  had maximum endpoint difference `3.013e-5`; independent score-equation
  checks support the tighter native inversion. These comparison tolerances do not
  replace the tighter defining-equation and acceptance-set tests.

The SAS source is the
[FREQ procedure manual](https://support.sas.com/documentation/onlinedoc/stat/current/freq.pdf).
The statistical papers are cited in each function's help page. Retain the
book's comparison code until the authors approve removing it; DescTools must
remain in the shared book/manual environment in the meantime.
