# R, Python, and Julia development principles

Agreed: 2026-11-06. Last updated: 2026-11-11.

This document records the agreed approach to the Python and Julia sister
packages for ibist. It is the authoritative policy for all three implementations.
The Python repository is [ibist-py](https://github.com/ibist/ibist-py).
Read [development lessons](DEVELOPMENT_LESSONS.md) before changing statistical
diagnostics or translating them; it records failures and required safeguards.
Its initial package structure uses `ibist` for both the distribution and import
names. The first implementation provides dataset loading and `or_2x2_ci()`.
These names describe the development package and do not imply a reserved
PyPI name.

## Shared statistical specification

Maintain one specification of each method and three native interfaces. The
specification defines the estimand, sampling assumptions, null hypothesis,
statistic, variance formulas, confidence interval construction, and numerical
algorithm. It also defines defaults, alternatives, table orientation,
reference categories, category order, continuity corrections, exact versus
mid-p calculations, boundary behavior, and handling of invalid inputs.

R is the reference implementation because it contains the existing methods
and validation examples. Develop new features in R first. Research and settle
the statistical definitions, interface, output, documentation, and focused
validation before translating the stabilized design to Python and Julia.
This sequence avoids revising all implementations while decisions remain open.
An implementation's behavior alone does not establish statistical correctness.

## Native interfaces

The Python package should work without an R installation. Use snake_case for
Python functions and arguments, for example `or_2x2_ci()`,
`agreement_test()`, `conf_level`, and `test_variance`. These examples illustrate
naming conventions; individual signatures still require design review.
Preserve the current R interface until its separately planned naming migration.

### Julia sister package

On 2026-11-11 the authors approved a native Julia sister package. Its module
and package name is `IBist`, using `using IBist`; the local repository is
`IBist.jl`, published at https://github.com/ibist/ibist-jl. The mixed-case
module spelling is deliberate. Registry availability does not reserve the
name; General registration remains a separate step.
Use snake_case functions and keyword arguments corresponding to Python,
including `or_2x2_ci()` and `conf_level`, with CamelCase result types.
Julia users must not require R or Python at runtime.

The first milestone is the seven shared teaching datasets and the four
2x2 odds-ratio interval methods. Preserve defaults, method order and
repetitions, fixed-margin conditioning, exact/mid-p definitions, one-sided
bounds, and unadjusted common estimates. `CIResult` provides full-precision
estimates and limits, readable printing, and Tables.jl/DataFrame views.
Interval-only results contain no fabricated hypothesis test.

DataFrames is a regular dependency so the teaching examples work after
standard installation. Julia factors use categorical arrays with the
exported levels and order; all missing values use `missing`. Numeric value
labels remain metadata. This differs deliberately from Python's numeric NaN
and nullable-integer representations. Each load and explicit DataFrame copy
is independent. Dataset metadata retain definitions and provenance.

`data-raw/export-julia.R` invokes the existing portable exporters with
Julia destinations; it neither duplicates cleaning nor changes Python's
exports. Its additional one-sided reference calls supplement the shared
two-sided OR fixtures. Julia's conditional limits use log-scale numerical
inversion, allowing `2e-6` relative or `2e-8` absolute endpoint differences
from R's direct odds-ratio inversion. Reuse regression cases alongside
independent formulas, published results, exact coverage enumeration, and
boundary checks. Approximation diagnostics and the remaining families are
deferred rather than translating unsettled behavior.

Use Julia-native numerical code for ordinary methods. Algorithms designated
as shared cores, including Wang's, must retain their authoritative C++
implementation, with a thin binding and portable binary artifacts when added.
Do not independently rewrite those algorithms in Julia. Tests cover Julia
1.10 and current stable on Linux, macOS, and Windows; build executable
documentation without automatic registration or publishing.

On 2026-11-06, the authors agreed to use `prop.1s.ci()`, `prop.2s.ci()`,
and `rr.2x2.ci()` in R, with `prop_1s_ci()`, `prop_2s_ci()`, and
`rr_2x2_ci()` as their Python counterparts. The new R proportion interfaces
follow base R's `x` and `n` convention. Both implementations should use the
same descriptive method names, including `newcombe`, `miettinen-nurminen`,
`agresti-caffo`, `katz`, and `koopman`, with continuity correction and mid-p
modification specified separately where applicable. See
[the binomial interval specification](BINOMIAL_INTERVALS.md) for defaults,
result fields, table orientation, and validation requirements. These interval
families were implemented in R first. `prop_1s_ci()` is now implemented in
Python, including vector inputs and all one-proportion methods.
`prop_2s_ci()` now implements the six approximate independent-difference
methods with scalar estimates and separate group metadata. `rr_2x2_ci()`
now provides Katz and Koopman intervals. Basic one- and two-sample Poisson
rate intervals and normal-approximation tests are also implemented in Python.
Chan--Zhang exact independent-difference intervals are also implemented.
Paired proportion, paired odds-ratio, and common odds-ratio intervals now
have Python counterparts with the same method names and CI-side interface.
Python uses zero-based input metadata and a DataFrame MultiIndex for vector
intervals, while scalar interval output retains its existing layout. See the
shared specification for numerical tolerances and deliberate differences.

### Function documentation

Keep Python function help comparable to the corresponding R help page in
statistical content and detail. Follow PEP 257 and numpydoc, with Parameters,
Returns, Raises, See Also, Notes, References, and Examples sections for public
statistical functions; include Warns where relevant.
Document every argument, effective default, table orientation, method,
boundary case, and relevant numerical limitation. State intentional interface
or result-field differences explicitly. Examples must run against the public
API. Dataset helpers use the applicable sections without adding irrelevant
statistical references. Reuse verified references from the shared method
specification and R help; verify new sources before citing them.

### Result structures

Use R's statistical result structures as the design basis for the Python
interfaces. `CIResult` corresponds to R's `ibist_ci`, and `TestResult` to
`htest`. These are ibist's package-wide conventions, with consistent field
semantics and documented method-specific information.

Translate dotted R field names to Python names, such as `conf.int` to
`conf_int`, `conf.level` to `conf_level`, `p.value` to `p_value`, and
`null.value` to `null_value`. Preserve the distinctions between parameter
estimates, sample metadata, test statistics, null hypotheses, and intervals.
State when an interval and a hypothesis test use different constructions.

Return structured results at full numerical precision. Summaries, interval
arrays, and DataFrame exports are views of those results. Keep method names,
defaults, and supported alternatives consistent across functions where their
statistical meaning is the same, and document intentional differences from R.
Evaluate other Python packages' conventions by statistical clarity and
consistency with this design.

### Risk ratios and basic Poisson rates

Risk-ratio methods follow `BINOMIAL_INTERVALS.md`: `rr_2x2_ci()` returns a
scalar estimate and copied group counts/proportions, supports one-sided
bounds, and preserves method order and repetitions. Katz requires positive
event counts; Koopman permits zero events, with an undefined estimate and
`[0, Inf]` when both counts are zero. Python uses complementary constrained
probabilities and failure counts to retain precision near complete success.

The basic rate interfaces are `rate_1s_ci()`, `rate_2s_ci()`, and
`rate_test()`, corresponding to R's dotted names. The one-rate default is
Garwood exact; other methods are score, Wilson--Hilferty (`wh`), Wald, and
log. `correct=True` is the one-rate default, matching R, and affects score,
WH, and Wald only. The two-rate default is log-Wald, with conditional Wilson
score and Clopper--Pearson exact alternatives. At least one event is needed
for conditional two-rate intervals; log-Wald requires both counts positive.
Counts and their total must not exceed `2**53-1`; exposures are finite and
positive and must match the count vector length, without scalar recycling.

Rate intervals use a scalar `CIResult.estimate`, with `samples` containing
`x` and `exposure`, and copied `sample_rates` arrays. This deliberately
separates the rate-ratio target from the group rates; R's two-rate estimate
vector contains all three. Scalar DataFrame output retains its method index,
and sample metadata attributes are immutable tuples. Method repetitions and
order are retained. Summaries round for presentation only.

All CI entry points use `side`: `"two.sided"` (default), `"lower"`, or
`"upper"`. Both packages use an omitted value (`NULL`/`None`) to select the
default while detecting conflicts with the compatibility alias `alternative`.
Every rate-CI method now honors the requested side in R and Python. Tests
retain the hypothesis argument `alternative`. The complete contract and
translation table are in [CI_INTERFACE.md](CI_INTERFACE.md).

`rate_test()` preserves R's normal approximation and correction rules.
One sample tests a positive null unit rate `r` with a signed Z statistic,
uncapped directional half-count correction, and score interval. Two samples
test equal rates with the conditional binomial Pearson statistic and capped
Yates correction; their rate-ratio interval is log-Wald, unaffected by
correction. Its interval and test are not exact duals. Both counts must be
positive for that interval. A two-sample small-expected-count warning follows
R. The existing `TestResult` holds the scalar estimate, null, confidence
limits, and native test fields; `estimate_name` labels rates without calling
them conditional odds ratios. `samples` and `sample_rates` return copies.

One deliberate boundary difference is recorded in
`misc/rate-interval-review.md`: Python restricts negative WH limits to zero,
while R's corrected WH lower limit at zero events can be negative. The
midpoint correction itself remains subject to author review. The new R
`side` argument preserves the existing two-sided WH formula; book sources
remain unchanged. Clipping is not a coverage proof.

Independent validation includes the published SAS Katz example, constrained
likelihood fits for Koopman endpoints, closed-form boundary limits, Poisson
tail equations, score equations, conditional binomial bounds, exposure
scaling, reversal, and correction rules. The exporter
`data-raw/export-python-risk-rate-reference.R` records R source hashes and
revision and writes 336 risk-ratio and 123 rate cases. These supplement
independent tests and preserve raw R limits, including the WH finding.
Python regression tests account explicitly for clipping; risk-ratio limits
allow `3e-8` relative or `2e-8` absolute difference, and rate limits use
`2e-10` relative and absolute tolerance.

### Attributable measures

The R addition `impact.2x2()` summarizes all four attributable measures
by default for cohort data: `risk-exposed`, `risk-population`,
`fraction-exposed`, and `fraction-population`. Rows are exposed/unexposed;
columns are event/nonevent. Optional `measure` selects rows in the supplied
order, including repetitions. The `ibist_impact` result retains named
`estimate`, a lower/upper `conf.int` matrix, per-row `method`, study design,
observed counts, and population-information metadata. `as.data.frame()`
exports the summary; print methods round only for display. Retain side,
alternative, confidence-level, and numerical-tolerance conventions.

Exposed risk uses the existing uncorrected Wald risk difference. Cohort
exposed fractions transform RR limits selected with `rr.method`: `katz`
(default) or `koopman`. Other measures use Wald intervals on their stated
scales. Negative estimates are retained, including protective fractions.
Risk limits use [-1,1]; fraction limits use [-Inf,1]. Undefined log-Wald
intervals retain identified estimates and NA limits, with a warning and
measure-specific `notes`. Other computable rows remain available. Koopman
permits extended limits with zero events. None of these computability
notes establish asymptotic adequacy.

Population cohort measures assume a representative joint multinomial
sample when exposure prevalence is estimated from the table. Include the
shared-data covariance in their delta-method variances. Optional known
population `exposure` instead uses conditional binomial variances and treats
that supplied proportion as fixed. Separately sampled exposure groups need
population exposure information for population inference.

Case-control counts alone do not identify absolute risks. Without external
outcome prevalence, require `rare=TRUE` and report only the two fractions.
The exposed fraction transforms sample-OR Wald limits; the population
fraction uses the exposed fraction of cases and the OR approximation, with
the log-complement variance from Jewell (2004), p. 90. The sampled case
fraction does not assess the rare-outcome assumption.

Known population outcome `prevalence` reconstructs risks by Bayes' rule
and permits all four case-control measures without the rare-outcome
approximation. With case/control exposure proportions q1/q0 and supplied
outcome prevalence r, e = r*q1+(1-r)*q0, p1 = r*q1/e, and
p0 = r*(1-q1)/(1-e). Delta-method variances propagate the independent
case/control binomial variances through these equations, including induced
covariance. Fractions use the log-complement scale. Supplied population
values are fixed: uncertainty in external estimates is not included.

Validation includes Hanley's population-fraction identities, Jewell's
published case-control example (95% limits rounding to 0.36 and 0.75), and
independent delta calculations using exposure prevalence and the two risks.
Known-prevalence case-control intervals are checked against numerical
Bayes-rule derivatives for a population with exposure .25 and risks .3/.1;
its outcome prevalence is .15, distinct from the sampled case fraction.
The VF cohort exposed fraction matches epiR 2.0.93 `AFRisk.strata.wald`.
Population risk matches its covariance-aware `PARisk.strata.piri`.
Population fractions match `PAFRisk.strata.wald` and `PAFest.strata.wald`
when obtained as individual strata, rather than the single-table crude
display. For example, the VF population fraction is 0.035626 with limits
(0.001924, 0.068189), whereas epiR's single-table display gives
(0.028706, 0.043138). Case-control exposed fractions use sample-OR Wald
limits; epiR's per-stratum `AFest.strata.wald` uses conditional-OR limits.
These constructions must be identified when comparing results.

R documentation cites the definitions and identifiability references.
Original book concepts and the epiR dependency are retained during review.
Python translation follows R stabilization.

### Contingency-table independence tests

On 2026-11-06, the authors chose `indep.test()` in R and `indep_test()` in
Python to provide a common interface to independence tests. Both
implementations support `pearson` (default), `lr`, and `fisher` methods.
The asymptotic methods prepare replacements for the book's use of DescTools for the
likelihood-ratio test and expected cell counts;
manuscript migration is deferred until after the final exam.

- Input is a two-dimensional table of finite, non-negative integer counts,
  with at least two rows and two columns, a finite total, and positive row
  and column totals. Reject empty margins rather than dropping them.
- Both methods test independence of row and column variables. Expected
  counts are `E[i,j] = row_total[i] * column_total[j] / total`. Rows and
  columns may be exchanged without changing the test.
- Pearson uses `sum((observed - expected)^2 / expected)`. The likelihood
  ratio uses `2 * sum(observed * log(observed / expected))`, treating
  zero observed counts as zero contributions. Both use the upper tail of
  a chi-squared distribution with `(nrow - 1) * (ncol - 1)` degrees of freedom.
- `correct = FALSE` is the common default, unlike R's `chisq.test()`.
  `correct = TRUE` is supported only for Pearson's test of a 2 x 2 table;
  other combinations are errors. Use R's capped Yates correction:
  subtract `min(0.5, min(abs(observed - expected)))` from each absolute
  deviation before squaring. Warn when any expected count is below 5.
- The R Pearson method wraps `stats::chisq.test()`; its expected counts and
  Pearson residual diagnostics are also reused by the likelihood-ratio
  method, whose statistic is computed directly. Python may use SciPy where
  its behavior matches this specification; custom calculations are needed
  for any differences, including the cap on the Yates correction.
- Return the statistic, degrees of freedom, asymptotic p-value, test
  description, observed and expected counts. R returns an `htest` with
  `statistic`, `parameter`, `p.value`, `method`, `data.name`, `observed`,
  `expected`, `residuals`, and `stdres`. The residual fields are Pearson
  diagnostics for either method, not deviance residuals. Preserve dimension
  labels. Python should expose corresponding native fields, including
  `statistic`, `df`, `p_value`, `method`, `observed`, and `expected`.
- Validate against direct calculations and an independent Poisson
  log-linear model as well as R's Pearson implementation. Include zero
  cells, independent tables, sparse-count warnings, transposition and
  permutation invariance, invalid inputs, and the capped Yates case.

Fisher's method wraps `stats::fisher.test()` for 2 x 2 and general RxC
tables. It uses the same input validation and adds observed and expected
counts to the returned `htest`; it does not call `chisq.test()` or issue
an asymptotic small-count warning. The null hypothesis is independence
(odds ratio one for 2 x 2 tables), with a two-sided exact test by default.

- Named R options forwarded through `...` are `alternative`, `conf.int`,
  `conf.level`, `workspace`, `simulate.p.value`, and `B`. Reject other names,
  duplicates, unnamed options, and Fisher options on an asymptotic method.
  Reject continuity correction for Fisher. Hybrid approximations and a
  non-unit null odds ratio are outside the independence interface.
- For 2 x 2 tables, support `two.sided`, `less`, and `greater` alternatives.
  Retain R's conditional maximum likelihood odds-ratio estimate and exact
  confidence interval (default confidence level 0.95); `conf.int = FALSE`
  omits the interval. The estimate is not the sample cross-product ratio,
  and the interval is not for a difference of proportions.
- For larger tables, support only two-sided inference. Use R's exact
  network algorithm by default and pass through workspace errors rather
  than silently substituting an approximation. Monte Carlo computation is
  available through `simulate.p.value = TRUE`, with `B = 2000` by default.
  It samples tables with the observed margins fixed and uses the plus-one
  p-value convention, giving a minimum of `1 / (B + 1)`. R's simulation
  option does not change the exact calculation for a 2 x 2 table.
- Preserve Fisher's own inference fields. There is no chi-squared statistic,
  degrees of freedom, or Pearson residual output for this method. Return
  observed and expected tables for all three methods, preserving labels.
- The Python Fisher implementation supports general RxC tables, not just
  2 x 2 tables, independently of SciPy's Fisher-test interface. It uses native
  fixed-margin enumeration for exact RxC tests, with an explicit search-node
  limit, and SciPy's Patefield sampler for Monte Carlo tests. Two-sided exact
  p-values sum probabilities of tables no more likely than the observed table under
  the fixed-margin null, rather than ordering tables by Pearson X-squared.
  R's `fisher.test()` is a cross-language reference, supplemented by direct
  hypergeometric and small RxC enumeration checks. Record probability tie
  tolerances, including R's relative `1e-7` tolerance for 2 x 2 tables.
  Seed reproducibility is required within each implementation; Monte Carlo
  samples need not match across languages.

Goodness-of-fit tests and likelihood-ratio continuity corrections remain
outside this interface.

Direct development comparisons with DescTools 0.99.60 covered both book
G-test examples, a zero-cell table, and 100 seeded random tables of sizes
2 x 2 through 5 x 5. Statistics and expected counts agreed within `1e-10`.
For the tumor table (rows `109,68,14; 72,79,37; 9,42,139`), DescTools returned
a zero p-value; the upper-tail calculation in `indep.test()` retains the
probability of approximately `5.73e-53`. Compute survival probabilities
directly, not by subtracting the lower-tail probability from one. DescTools
is not a runtime or test dependency; retain independent formula and Poisson
model checks in the test suite.

Python uses SciPy's conditional odds-ratio estimate and exact limits for
binary Fisher tests. These satisfy conditional mean and inclusive tail
checks at absolute tolerance `1e-10`; R's default `fisher.test()` root finding
can give differences of about `1.2e-4` relatively for estimates and `5e-4` for
limits in the reference cases. Keep the statistical definitions aligned,
rather than reproducing R's coarse numerical roots. Native exact RxC Fisher
enumeration is less efficient than R's network algorithm; its `max_tables`
limit raises an error instead of substituting a different method. Python
simulation options are `n_resamples` and `random_state`. General-table tests
use the same relative probability tie tolerance `1e-7` as the binary tests.

### Exact common odds-ratio tests

`or.common.test(x, theta = 1, alternative = "two.sided")` has the Python
counterpart `or_common_test()`, using the distribution engine shared with
`or_common_ci()`. Input is a list of 2 x 2 tables or a 2 x 2 x K array.
The null specifies the same positive finite OR in every stratum. Conditional
on the margins, combine noncentral hypergeometric distributions at that OR.
Use inclusive tails of the total upper-left count; the two-sided p-value is
twice the smaller tail, capped at one. This follows Jung's (2014) stratified
Fisher tail formulation, with the non-unit-null extension used in R.

Return `TestResult` with the raw integer statistic, p-value, alternative,
conditional common-OR MLE, `null_value=theta`, and stratum count. R's
`parameter` duplicates the null OR; do not translate it into degrees of
freedom. Python uses `df=None` and a separate `n_strata`. Confidence limits
remain in the CI function. The MLE depends on the data, independently of
the tested null. Empty-margin/deterministic strata follow the existing CI
engine, including R's estimate-zero convention for an unidentified OR.
Python retains exact integer count totals and sums tails in log space.

Independent checks include full allocation enumeration, hypergeometric
tails, exact conditional null rejection probabilities, CI/test inversion,
reversal invariance, the primary paper's Table 2, and tiny-tail formulas.
`data-raw/export-python-common-test-reference.R` exports 456 further R
comparisons, including nonzero support minima, boundary estimates, and
non-unit nulls. P-value tolerance is `2e-10` relative without an absolute
allowance; conditional MLE tolerance is `2e-9` relative or `1e-10` absolute.
The accompanying R bibliography correction uses the verified reference:
Jung, S.-H. (2014), Biometrical Journal 56(1), 129--140,
DOI `10.1002/bimj.201300048`.

### Nominal association measures

The R interface is `nominal.assoc(x, measure = c("phi", "contingency",
"cramer-v"))`; the Python counterpart is `nominal_assoc()` with the
same measure names. All three are selected by default; preserve requested
order and remove duplicates. Both are implemented. Book changes
remain deferred until after the final exam.

Input requirements and expected counts match `indep.test()`, including
rejection of empty margins. Let `X2` be the uncorrected Pearson statistic,
`n` the total count, and `d = min(nrow, ncol)`. Return:

- `phi = sqrt(X2 / n)`, non-negative even for 2 x 2 tables, with maximum
  `sqrt(d - 1)`; do not substitute a signed binary correlation.
- `contingency = sqrt(X2 / (X2 + n))`, without Sakoda adjustment, with
  maximum `sqrt((d - 1) / d)`.
- `cramer-v = sqrt(X2 / (n * (d - 1)))`, without bias correction, between
  zero and one. Phi and Cramer's V agree when `d = 2`.

These definitions match the default point estimates of DescTools `Phi()`,
`ContCoef()`, and `CramerV()`. Do not apply Yates correction or use G-squared.
These descriptive measures ignore category ordering and are invariant to
transposition, category permutations, and uniform replication of counts.
Suppress sparse-count approximation warnings because no p-value is returned.
Confidence intervals, bias corrections, and adjusted contingency coefficients
are outside this first implementation.

R returns an `ibist_nominal_assoc` object with `estimate` (a named vector),
`measure`, `statistic` (uncorrected Pearson X-squared), `n`, `observed`,
`expected`, and `data.name`. Its print method rounds only for display. The
Python result should expose corresponding native fields. Validate against
independent 2 x 2 cross-product formulas, exact independence and perfect
association cases, and direct DescTools comparisons, including the book's
3 x 3 tumor table, rectangular tables, sparse tables, and Phi above one.
DescTools remains a development comparison only, not a package dependency.

Direct comparisons against DescTools 0.99.60 matched all three estimates
exactly on 205 tables: the book's tumor table, positive and negative binary
association, perfect 3 x 3 association, a sparse zero-cell table, and 200
seeded random tables with dimensions from 2 through 6. The tumor estimates
are `phi = 0.63761827201204579`, `contingency = 0.53762824441746093`, and
`cramer-v = 0.45086420394816618`; these are retained as regression references
alongside independent formula checks.

### Unified association tests

`assoc.test()` and `assoc_test()` accept one or more measures and one method
for the whole call. Defaults are Pearson score correlation, multinomial
inference, two-sided alternatives, and confidence level 0.95. Existing
functions and book examples retain their interfaces without deprecation.
Agreement and trend interfaces are separately planned.

- `method="pearson"` supports nominal `phi`, `contingency`, and `cramer-v`,
  sharing one Pearson independence test. Estimates use the uncorrected
  statistic even when `correct=TRUE`/`True` requests the 2 x 2 capped Yates
  correction for testing. Only two-sided, directionless independence tests
  are supported; scores and empty margins are rejected. No coefficient CI
  is supplied. Small expected counts warn once per call.
- `method="multinomial"` supports all seven ordinal measures and preserves
  the existing ordinal estimates, estimation standard errors, intervals,
  null variances, tails, and score rules. The five concordance-based
  measures share a test. Pearson and Spearman have separate null variances.
- `method="cmh"` supports Pearson and Spearman, with signed Z=r*sqrt(n-1)
  and its square in `chisq`. No CI or estimation standard error is supplied;
  do not substitute the multinomial null variance. Scores must be supplied
  for both axes and affect Pearson only, including mixed requests. Reject
  scores for Spearman-only CMH calls. Ordinal methods drop empty margins.

Reject incompatible combinations, including nominal/ordinal mixtures;
never select different methods silently. Confidence levels are validated
for all methods but used only for multinomial intervals. One measure,
including a length-one vector/sequence, returns an htest/TestResult. Add
`measure` and `test.method`/`test_method` identifiers; standard errors and
CMH squares are method-specific optional fields. Null coefficients are zero.

Multiple measures return an ordered list/sequence preserving repetitions:
R's `ibist_assoc_tests` and Python's `AssociationTestResults`. Provide native
individual tests and a full-precision combined table, with shared-test groups
identified in metadata and printing. State that p-values are unadjusted and
tests use the same observations. Python's new optional TestResult fields
are appended, preserving existing positional arguments and defaults.
Independent tests check Pearson/Yates formulas, expanded score correlations,
CMH normal tails, the primary SAS example, invalid combinations, preserved
legacy results, and tiny p-values. The R exporter adds 274 comparison calls;
tolerances follow the ordinal specification below. Book edits remain deferred.

### Two-rater agreement tests

`agreement.test()` and `agreement_test()` initially provide Cohen's kappa
and weighted kappa. Their arguments mirror the existing R function:
`x`, `weights="unweighted"`, confidence level 0.95, two-sided alternative,
and `test.variance="fce"` / `test_variance="fce"`. One weighting scheme is
selected per call: unweighted, linear, quadratic, or a custom matrix.
Named weights use category positions. Retain empty rating categories and
require matching category order on both axes; do not sort or align labels.
Counts and their total must not exceed 2^53-1 in the new entry points.

R delegates to `cohen.kappa.test()` and retains its numerical results and
positional arguments. Keep that legacy function fully working without
deprecation warnings until after the final exam. Deprecation is a later,
separate change; do not activate it automatically by calendar date. Book
sources remain unchanged.
At deprecation, route both entry points through an internal implementation
so that the new function never emits the legacy deprecation warning.

Use the observed FCE estimation variance for CIs. The default FCE test
uses its distinct null variance; delta testing uses the observed estimation
variance instead. Preserve `null.std.error` / `null_std_error` as the
selected test SE, explicitly documenting that delta supplies an estimation
SE rather than a null SE. Return htest/TestResult with the existing kappa
inference fields, standard errors, selected variance, weight matrix,
weight scheme, and sample size. Append optional Python fields to preserve
existing positional TestResult arguments; weights are immutable internal
tuples exposed as fresh NumPy arrays.

Match R's custom-matrix symmetry tolerance and reject chance agreement one
or invalid/nonpositive test variances. Two-sided Wald intervals are
untruncated; one-sided finite limits are clipped to [-1, 1]. Exact and
multi-rater inference are outside this implementation. Independent checks
use numerical multinomial derivatives, closed-form binary variances, the
primary SAS skin example, invariance, replication, and tiny normal tails.
`data-raw/export-python-agreement-reference.R` adds 60 cases including
the book's pericardial-effusion example and expected errors. Use `2e-9`
relative tolerance with no absolute allowance for p-values; ordinary
estimates/statistics allow `1e-11` absolute error. Near-zero standard errors
allow `1e-8`, intervals `3e-8`, and serialized weight matrices `1e-15`.

### Trend and power

`trend.test()` / `trend_test()` provide asymptotic and exact conditional
Cochran--Armitage inference. Keep `ca.trend.test()` supported. Rows are
ordered groups, column one contains events, scores default to original row
positions, and empty rows drop after score assignment. Return signed Z,
its square, the raw score statistic, p-value, alternative, and description.
Exact tails are inclusive and two-sided p-values double the smaller tail.
Python follows R's 15-significant-digit score merging and tie tolerance,
with a one-million-state resource limit absent from R's legacy backend.

`power_p1s_test()` and `power_p2s_test()` mirror the existing R proportion
power functions. Exactly one of n, relevant proportions, power, or
significance level is omitted. `PowerResult` follows `power.htest` fields,
translating dots to underscores. Preserve correction/strict defaults,
alternative names, root tolerances, and first-group allocation semantics.
Supplied one-sample n is rounded up; normal solved n and two-sample n
remain continuous. Exact one-sample methods are quantile, deterministic
midp, randomized midp-rand boundaries, and inclusive cp rejection.

The R exact sample-size search was corrected before translation: scan from
one in order, honor max_n, and check actual size plus target power. The
previous local window could return 5 when 1 sufficed or max_n was 4.
Exact power is not monotone in n. Preserve close/minimal rules and report
achieved significance level, including randomized rejection contributions.
For normal one-sample designs, max_n remains an extensible initial bracket,
as in the legacy solver. Exact p0/level inversions are discontinuous and
can locate a rejection boundary without attaining the requested power;
document R's target-power convention and explicit recomputation.
Python two-sample probability inversion stays within admissible domains.

Focused validation uses direct trend enumeration, closed-form normal
designs, independently summed binomial rejection regions, and 30 R cases.
Allow 2e-4 absolute tolerance for solved parameters because R's two-sample
root tolerance is eps^0.25. Broader validation remains reserved for a push.

`power.r1s.test()` / `power_r1s_test()` and `power.r2s.test()` /
`power_r2s_test()` add normal-approximation Poisson rate planning. Omit
exactly one of n, relevant rates, power, or significance level. Exposure
is fixed per-subject follow-up; total exposure is n times follow-up.
Two-sample group.rate/group_rate is n2/n1, independently of follow-up;
exposure accepts a shared scalar or separate group values. Both supplied
and solved n remain continuous. Results retain rates, exposure, allocation,
level, power, alternative, method, and note in power.htest/PowerResult.

Use rate.test() score boundaries and the alternative Poisson variance.
Two-sample null variance uses the exposure-weighted expected pooled rate;
this replaces the test's random pooled rate for planning. Help pages give
the full derivation. Defaults match proportion power: correction and strict
are false, level is .05. Strict includes both two-sided tails. One-sample
alternatives are two.sided/less/greater; two-sample alternatives are
two.sided/one.sided, using the absolute difference. These are approximate
normal calculations, not exact Poisson or conditional-test power.

Positive inverses use log coordinates within 100 doublings/halvings from
1 for n or the other rate; level uses logit coordinates of magnitude at
most 36. Default root tolerance is sqrt(eps). Prefer r1 above r0 or r0
below r1 for one sample (reverse for less), and r1 below r2 or r2 above r1
for two samples; try the other branch if needed. Reject equal-rate sample
size designs and unattainable targets. Focused checks use closed-form
sample sizes, inverse recovery, unit/allocation invariance, and 19 R
cases, with 5e-8 relative or 1e-10 absolute comparison tolerance.
Published Minitab examples additionally verify four powers at 5e-7
absolute tolerance and two rounded sample sizes in both languages. Match
the uncorrected normal method and request strict two-sided power.
See [published rate-power validation](RATE_POWER_VALIDATION.md) for inputs,
sources, numerical comparisons, and the limits of this external check.
One-sample solved parameters allow 2e-8 absolute error; exact null size
derived from a fitted p0 allows 5e-8 because it inherits root-finding error.
The CLT helper name remains open for discussion; do not translate it yet.

### GLM goodness of fit

The agreed names are glm.gof() and glm_gof(), with the fitted model first.
Pearson and deviance support binomial and Poisson GLMs; Hosmer-Lemeshow,
Osius-Rojek, and Stukel require binomial logit models. One method returns
htest/TestResult, multiple methods return ordered ibist_glm_gof_tests /
GLMGOFResults collections. Preserve order and repetitions; p-values are
unadjusted. Default method is pearson; g=10 controls only Hosmer-Lemeshow.

Use integer counts/trial/frequency weights, grouped by identical design
rows AND offsets. Exclude zero-weight rows, including their boundary
predictions. Reject unsupported families, estimated dispersion, missing
retained inputs, and active boundary means. No automatic pooling or method
switching occurs. R's logistic.gof() retains its original implementation
and output; pooling rules are not transferred into the new interface.

All five methods are asymptotic. Use actual groups and estimable rank;
Hosmer-Lemeshow uses SAS target-frequency groups and probability ties at
14 significant digits by default (grouping="target"). R also supports
grouping="quantile"
using unique type-7 sample quantiles, right-closed groups, and frequency
weights interpreted as repeated subjects. A quantile within a tied mass
puts the whole mass in the group ending at that probability. Empty groups
are dropped, and df is actual group count minus two. Retain the selected
rule in results and diagnostics. Quantile grouping applies only to HL;
Python translation remains deferred until R stabilizes.
Osius-Rojek uses SAS J-centering. Stukel uses an
efficient score test, not an augmented likelihood-ratio fit, and accounts
for a rank-one increment when necessary. Column-scaled SVD uses 1e-10
relative tolerance; standardizing variances <=1e-12 are numerically zero.

Support asymcheck(glm.gof/glm_gof, fit, ...) without computing test
statistics and asymcheck(result) after inference. Optional asymcheck=TRUE
attaches the same diagnostics. Retain diagnostic inputs independently of
the visible flag and never retain the full fit. Pre-checks retain positive
df and usable standardization under details$computability; inference rejects
uncomputable tests. These quantities are not adequacy criteria. Report
convergence, rank, fitted boundary distance, original groups, and actual HL
risk groups separately from approximation assessment.

The 2026-11-09 book audit implemented the following corrected R diagnostics;
Python translation and shared-reference updates remain deferred until R
stabilizes. Elizabeth's editorial review of the proposed cross-chapter
binomial guideline is pending in the book's Calibration section. Binomial
Pearson/deviance apply Ch. 11's n_g*p.hat_g*(1-p.hat_g) >=5 within each
Ch. 16 covariate-pattern group. HL applies this to risk groups and requires
actual H >=6 after ties, rather than merely positive df. The binomial
variance screen explicitly operationalizes Ch. 16's qualitative ``large
groups'' condition; it is not a numeric cutoff stated in that chapter.
Report group sizes and G/n without inventing a ratio threshold.

Poisson Pearson/deviance use all fitted expected counts >=5 as Ch. 23's
conservative screen. Report the number below 2 and fraction below 5 without
inventing a cutoff for ``many.'' Residual df >50 concerns sensitivity to
minor departures, not failed asymptotic adequacy. Osius-Rojek and Stukel
remain asymptotic and applicable but have condition.met=NA: the reviewed
references establish no finite-sample adequacy cutoff. Osius-Rojek permits
individual groups; its variance-growth and other limiting conditions cannot
be certified from positive fitted variance. Hosmer et al. (1997) describes
conditional-moment and scaled chi-square refinements that are not part of
the current method. Computability must not turn unassessed adequacy into
a pass. A collection is FALSE if any method fails, NA if none fails but any
is unassessed, and TRUE only if all pass. Test statistics, p-values,
grouping, validation, and the old logistic.gof() algorithms are unchanged.

Python consumes statsmodels GLM results and requires statsmodels>=0.14.
Its new TestResult fields are appended for positional compatibility.
Published SAS examples validate all five logistic methods, grouped
binomial responses, and Poisson offsets; eight R cases and independent
formula/score checks cover interfaces and numerical behavior. See
[GLM GOF validation](GLM_GOF_VALIDATION.md) for sources, variants, data
preparation, numerical tolerances, and limits of the evidence.

### Shared approximation diagnostics

`asymcheck()` is the agreed R and Python name. A function plus its usual
inputs gives a pre-check without inference; a supported result gives a
post-check using retained inputs and settings. Initial function adapters
cover indep.test/indep_test and rate.1s.ci/rate_1s_ci. R additionally
supports stats::mantelhaen.test pre-checks on 2 x 2 x K count arrays.
Raw counts plus method="mantel-haenszel" work in both languages. Do not
evaluate data.name to recover missing inputs from arbitrary results.

Both inference families accept asymcheck=FALSE by default; TRUE attaches
the same result under diagnostics. Keep default R Pearson results identical
to chisq.test(); existing observed counts and method descriptions suffice
for post-checks. R rate CIs retain an inputs data frame and correct flag;
Python already retains rate counts/exposure, correction, and orientation.
Append optional diagnostics fields to Python result classes to preserve
positional arguments. Neither diagnostics nor their display change inference.

Return ibist_asymcheck/AsymCheckResult with procedure, method, applicable,
checks, details, notes, and condition.met/condition_met. Criterion tables
record method, criterion, value, threshold, comparison, and whether met.
Unassessed values are NA/None; summaries cover only assessed criteria.
Expected-count screening follows book Ch. 12: all four expected counts >=5
for 2x2 tables; Ch. 14 uses minimum >=1 and fraction >=5 of at least .8 for
larger tables. LR reports the same
sparsity indicators with an explicit caveat about Pearson-derived rules.
Fisher inference is marked not applicable. Mantel-Haenszel preserves the
legacy variance-sum >=5 heuristic and zero contributions for empty strata.

One-rate intervals follow book Ch. 21: ordinary Wald observed x >=50;
log-Wald and uncorrected score x >=10. WH is recommended even for small
counts; no binary cutoff is invented for WH or corrected score. Their
criterion remains unassessed. Plug-in skewness 1/sqrt(x) is descriptive,
not a pass/fail approximation criterion. Exact intervals and
zero-event log/uncorrected Wald exact fallbacks are not applicable. Preserve
method order, repeated methods, correction, confidence level, and CI side.
Sampling independence and common-effect homogeneity are not assessed.
Keep mantelhaen.asymcheck working as a wrapper with its original output;
deprecation warnings and book edits remain deferred until after the exam.
Focused checks use direct formulas, pre/attached/post equivalence, pre-checks
without inference, and 14 R reference cases (1e-12 relative tolerance).

The book is the source for approximation criteria. Input validity and
computability are separate and must not produce a successful approximation
assessment. Report the book quantity, threshold, comparison, and result;
passing a guideline is not a coverage guarantee. Pseudo-counts and
continuity corrections do not increase the observed quantities used in the
guidelines. These criteria were settled in R first, then translated to
Python with the same quantities, thresholds, defaults, and settings.

| Procedure | Book criterion | Source |
| --- | --- | --- |
| One-proportion test | n*pi0*(1-pi0) >=5 | Ch. 11 score-test summary |
| Two-proportion equality test | Four pooled expected counts >=5 | Ch. 12 |
| Proportion CIs | n_i*p.hat_i*(1-p.hat_i) >=5 per sample | Ch. 11--12 |
| RR/OR approximate CIs | Same original-sample criterion | Ch. 15 RR/OR |
| MH test | Sum of null conditional variances >=5 | Ch. 15 MH summary |
| One-rate test | Null expected events r*T >=10 | Ch. 21 test summary |
| Two-rate equality test | n*pi0*(1-pi0) >=5 | Ch. 22 test summary |
| Log-Wald rate-ratio CI | Both original observed counts >10 | Ch. 22 |
| Conditional score rate-ratio CI | Binomial CI criterion | Ch. 22 and 11 |
| Exact/mid-p procedures | No asymptotic approximation check | Exact sections |

For the two-rate equality test, n=x1+x2 and pi0=T1/(T1+T2). A test result's
check applies to its test, not its accompanying CI. Function pre-checks
support stats::prop.test and stats::binom.test (one/two-sample equality
settings only); base R results lack retained counts and need pre-checks.
Package prop.1s.ci, prop.2s.ci, rr.2x2.ci, or.2x2.ci, rate.2s.ci, and
rate.test and their Python counterparts append asymcheck=FALSE/False and
retain original inputs and settings
for optional attached diagnostics and post-checks. Preserve method order,
repetitions, correction/mid-p settings, and CI side. OR pre-checks permit
named tol through ...; no interval inversion or enumeration is run.
Positive original sample sizes are required to check binomial CI criteria.
Book examples used for focused checks include OSA 16/40, VF 38/100 versus
267/984, AF observed 411 against .0119*25800 null expected events, and
the low-HDL and MI two-rate examples. Small positive tables must fail
when their book approximation quantity is below its threshold.
The exporter data-raw/export-python-book-asymcheck-reference.R supplies
15 representative CI/test cases with R revision and source hashes. It
translates field names and input/sample indices only; existing 14 shared
independence/MH/one-rate cases are regenerated with the book's criteria.
Python root tolerance for RR now matches R's strict (0,1) interval. Focused
comparisons use 1e-12 relative and 1e-14 absolute tolerances; direct book
quantities and pre/attached/post checks supplement cross-language agreement.

### Ordinal association

`ordinal.assoc()` has the Python counterpart `ordinal_assoc()`, with all
seven measures selected by default: `gamma`, `tau-b`, `tau-c`, `somers-c|r`,
`somers-r|c`, `pearson`, and `spearman`. Preserve requested order and
repetitions. Return `OrdinalAssociationResult`, corresponding to R's
`ibist_ordinal_assoc`, with named estimates, standard errors, confidence
limits, statistics, p-values, and concordance metadata. Python uses Series
for named vectors and a lower/upper DataFrame for intervals. Dotted field
names use underscores; `measures` corresponds to R's `measure`.

Keep physical row and column order; do not sort category labels. Drop empty
margins, then require at least two nonempty rows and columns. Validate
explicit Pearson scores against the original dimensions before subsetting.
Default Pearson scores use retained numeric labels or category positions.
Array input has no labels; DataFrame numeric indexes are labels, including
gaps after dropping empty categories. Spearman uses marginal midranks,
including their estimated-margin contribution to sampling variances.

The five concordance measures share `S/sqrt(var.S)`. Pearson and Spearman
use their own multinomial null variances; do not substitute continuous-normal
correlation inference or the separate CMH linear-by-linear test. Match R's
`alternative` and `conf.level` defaults. Intervals use estimation variances:
two-sided Wald intervals are untruncated; one-sided finite bounds are clipped
to [-1, 1], with the other endpoint at 1 or -1. Zero null variance produces
missing tests. These combined estimation/test functions retain `alternative`
without adding the CI-only `side` argument.

Python sums integer pair counts exactly and uses centered influence
variances to reduce cancellation. Category score rescaling avoids overflow
and underflow without changing correlations. Independent tests enumerate
pairs, differentiate estimands under multinomial sampling, reproduce all
seven measures in the SAS pain example, and check reversal, transposition,
replication, empty margins, score choices, and numerical boundaries.
`data-raw/export-python-ordinal-reference.R` adds 222 R comparison cases.
Use `2e-9` relative tolerance with no absolute allowance for p-values;
other fields allow `1e-11` absolute error. R's cancellation near perfect
association can leave standard errors near `1e-9`; allow `1e-8` absolute
error for those near-zero standard errors and `3e-8` for affected intervals.
Statistical references remain those in R's help; Brown and Benedetti (1977)
and the primary SAS example were checked while implementing the translation.

### Odds-ratio homogeneity tests

`or.homogeneity.test(x, method = c("zelen", "breslow-day", "woolf"),
correct = FALSE)` selects one method per call and preserves Zelen's exact
test as the default. The Python counterpart is `or_homogeneity_test()`,
with the same method names and correction policy. Both are implemented;
book changes remain deferred until after the final exam.

All methods accept a list of 2 x 2 count tables or a 2 x 2 x K count array.
The two asymptotic methods require at least two strata and positive row
and column margins in every stratum; reject empty margins rather than
silently dropping strata. The existing exact method's behavior is unchanged.

- Breslow--Day uses the Mantel--Haenszel common odds ratio, fits each
  stratum's cell counts under that odds ratio and the observed margins,
  and sums squared upper-left residuals divided by their fitted variances.
  Each variance is the reciprocal of the sum of reciprocals of the four
  fitted counts. Require a finite positive MH estimate; otherwise suggest
  the exact method. Solve the fitted-count quadratic stably, including
  the linear case when the MH estimate is exactly one.
- `correct = TRUE` requests Tarone's adjustment for Breslow--Day, subtracting
  the squared sum of upper-left residuals divided by the sum of their
  variances. This is not a continuity correction. Reject this flag for
  Zelen and Woolf rather than silently ignoring it.
- Woolf uses the inverse-variance weighted dispersion of stratum log odds
  ratios. Each variance is the sum of reciprocals of the four cell counts.
  To match DescTools, if any cell is zero, add 0.5 to every cell in every
  stratum before computing the log odds ratios and variances. This adjustment
  is automatic and independent of the `correct` flag.
- Both asymptotic methods return an `htest` with a chi-squared `statistic`,
  `parameter` named `df` equal to K - 1, an upper-tail `p.value`, a two-sided
  alternative, test description, and input description. Use the survival
  function directly to retain very small p-values. Warn if any fitted
  Breslow--Day cell or adjusted Woolf cell is below 5. Retain Zelen's existing
  statistic, parameter, and exact p-value structure.

Validate against independent common-OR-one and two-stratum Woolf formulas,
the book's postpartum-haemorrhage example, direct DescTools comparisons,
and transposition, simultaneous row reversal, and stratum permutation.
DescTools is a development comparison only, not a package dependency.

Direct comparisons with DescTools 0.99.60 covered 303 sets of strata,
including the book example, a common odds ratio exactly equal to one,
a zero-cell example, and 300 seeded random sets of 2 through 5 strata.
Breslow--Day, Tarone-adjusted Breslow--Day, and Woolf all agreed: maximum
absolute statistic differences were below `3e-11` and p-value differences
below `1.1e-13`. Retain the book's reference results in the regression tests
alongside independent formula checks.

### Numerical foundations and results

Use NumPy and SciPy for numerical building blocks. Add dependencies when the
functionality needs them; do not require a statistical modeling framework
merely to implement table-based inference. Follow the same statistical
definitions across languages even when underlying libraries have different
defaults.

pandas is a regular Python dependency: a standard installation must support
the teaching examples and DataFrame output. Keep development tools in an
optional dependency group, without requiring users to select an extra merely
to run the examples.

Return structured Python results with accessible estimates, statistics,
p-values, confidence limits, confidence level, and method information as
applicable. Provide readable summaries and tabular results for requests with
multiple methods or measures. Preserve full numerical precision internally;
round only for presentation. Track the correspondence between R and Python
functions and document any intentional differences.

## Common data sources and preparation

Maintain one authoritative collection of raw inputs, provenance, cleaning
rules, and variable definitions. Initially, the existing R `data-raw/` scripts
remain the preparation workflow. They may generate both R `.rda` files and
portable distribution files for Python. Python data preparation does not need
to duplicate the R cleaning code, and installed Python users should not need
R to load the prepared datasets.

Portable exports must preserve or accompany metadata for variable labels,
units, category levels and order, reference categories, and missing values.
CSV alone is insufficient to preserve all of this information. Record the
export format and metadata schema whenever they change.

The initial implementation uses `data-raw/export-python.R` to export the
existing prepared R datasets as gzip-compressed JSON column arrays and JSON
metadata. Missing values are null; metadata retain factor order, numeric value
labels, descriptions, units, and source-file hashes. Regenerate R objects with
their existing preparation scripts before exporting changed data. The exporter
checks an exact value round-trip and records raw inputs missing from the
checkout; export alone does not claim those inputs were reprocessed. See
[the Python dataset documentation](https://github.com/ibist/ibist-py/blob/main/docs/datasets.md)
for the schema. R comparison cases are generated explicitly with
`data-raw/export-python-reference.R` and stored in the Python test suite;
they supplement independent numerical checks.

Distribute generated copies with each package, but edit the authoritative
source and regenerate both outputs when data change. Record source and
preparation versions so that each package's data can be traced to its inputs.
Compare values, row and column counts, types, categories, and missingness
across the two exports. Release timing may differ, so record which source
revision each release includes.

## Feature workflow and validation

Thorough manual validation is reserved for immediately before a push.
During development, use only the focused checks needed for the current
change; do not repeat full suites, package checks, or distribution matrices
for each addition or commit. Batch additions and run the thorough checks once
before pushing. After a failure, rerun the affected checks; repeat broader
checks only when the fix invalidates their earlier results. Keep reference
cases representative and avoid redundant Cartesian grids.

1. Agree on the statistical specification, public interfaces, defaults, and
   reference examples before implementing a feature.
2. Research and implement the feature in R first.
3. Stabilize the R interface, output, documentation, and focused validation.
   Check published results, direct formulas, independently constructed
   examples, or established software with matching assumptions. Retain the
   published SAS and Jung examples already used by R where applicable.
4. Translate the settled R design to Python using the same specification and
   reference cases, supplemented by independent numerical validation.
5. Compare R and Python results using shared cases and explicit numerical
   tolerances. Include sparse tables, zeros, boundary estimates, one-sided
   alternatives, and multiple-measure requests where relevant. Agreement
   between the two implementations supplements independent validation.
6. Update documentation, examples, the interface correspondence, and affected
   book material. Record any feature available in only one implementation.

Generate portable reference cases from documented inputs and independently
validated expectations. Do not silently replace expected values with the
latest output of either implementation when a comparison fails.

## Computational code

Keep one authoritative, language-neutral implementation of shared native
algorithms. Language bindings handle validation, object conversion, errors,
and missing-value representations; the core uses standard C++ types and
has no Rcpp or Python dependencies.

For Wang intervals, the authoritative core is `src/wang_core.h` in `ibist-R`.
`src/wang_paired_ci.cpp` is the Rcpp binding. The exporter
`data-raw/export-python-wang-core.R` copies the header byte-for-byte into
`ibist-py/src/ibist/_native/wang_core.h` and records the source path, MD5 hash,
and license in `wang_core.json`. Edit the R source, then regenerate the Python
copy; never maintain separate handwritten algorithms. Both repositories
require LF line endings for this header to preserve byte identity on Windows.

The CPython binding is `_native/wang_module.cpp`; setuptools compiles it
into `ibist._wang_native`. Source distributions include the header and
binding; binary wheels include the compiled extension. Neither package
requires the other at installation or runtime. Source installation requires
a C++ compiler; wheels are specific to their Python version and platform.
The current publishing workflows build a Linux Python 3.12 wheel and a source
distribution; broader binary-wheel publishing is separate release work.

Run `Rscript --vanilla data-raw/export-python-wang-core.R --check` from the R
root, or `python tools/check_shared_core.py --r-source ../ibist-R` from the
Python root. Python CI checks out the authoritative R repository, compares
both copies and the manifest hash, and tests builds on Linux, macOS, and
Windows. Land the R core changes before the dependent Python CI changes.
The manifest hash identifies the exact core version shipped in each package.

Validate the shared core against independently calculated cases and finite
sample-space coverage, then validate each binding's interface. Cross-language
agreement alone is not an independent algorithm check. The refactored Wang
core also agrees bit-for-bit with the original Rcpp implementation for all
collapsed tables at n=0 through 6 and three confidence levels. Python's
Wang nuisance grids and precision settings preserve R's construction; finite
grids approximate the supremum and do not certify exact coverage.

Paired score limits use accurate constrained likelihood fits in Python;
R's existing optimizer is coarser, so comparison tolerance is `3e-4` absolute
for this method. Other remaining-CI references use `2e-7` relative or `2e-8`
absolute tolerance; Wang endpoints match exactly at the recorded precision.
The exporter `data-raw/export-python-remaining-ci-reference.R` supplies 545
cases covering all newly ported CI families. Independent checks include
conditional hypergeometric tail equations and conditional MLEs, multinomial
likelihood fits, small-space exact coverage, boundary closed forms, and a
Chan--Zhang disconnected acceptance region.

Porting exposed R edge cases: zero paired sample sizes now raise a clear
error; completely discordant samples use the analytic score variance at
boundary fits; low-confidence uncorrected Wilson bounds snap to exact
probability boundaries before odds transformations. Existing two-sided book
calls remain valid. Book sources are unchanged.

## Initial implementation sequence

1. Establish packaging, dataset loading, result objects, and validation;
   implement basic odds-ratio, rate, and paired-data intervals.
2. Add ordinal association, kappa, trend tests, and exact stratified inference.
3. Add remaining power calculations, logistic goodness of fit, and
   computationally demanding paired methods.

The first complete implementation should combine datasets and `or_2x2_ci()`,
including documentation and cross-language validation, to establish the
conventions before translating the remaining methods.

## Maintaining these decisions

Update this document whenever the authors agree to change a principle or
resolve an open design choice. Update the last-updated date and explain the
reason for a substantive change. Keep one authoritative copy and link to it
from the Python repository. Keep temporary
progress and repository state in the project status document.

On 2026-11-06, the authors approved creating `ibist/ibist-py`; this resolves
the repository name. They subsequently chose `ibist` for both distribution
and import names to keep installation and usage consistent. A PyPI check on
that date found no project named `ibist`; the name is not reserved.

Open implementation decisions include Python package-name availability,
release coordination, and the eventual location of shared source material.
