# R and Python development principles

Agreed: 2026-11-06. Last updated: 2026-11-07.

This document records the agreed approach to a Python sister package for
ibist. It is the authoritative development policy for both implementations.
The Python repository is [ibist-py](https://github.com/ibist/ibist-py).
Its initial package structure uses `ibist` for both the distribution and import
names. The first implementation provides dataset loading and `or_2x2_ci()`.
These names describe the development package and do not imply a reserved
PyPI name.

## Shared statistical specification

Maintain one specification of each method and two native interfaces. The
specification defines the estimand, sampling assumptions, null hypothesis,
statistic, variance formulas, confidence interval construction, and numerical
algorithm. It also defines defaults, alternatives, table orientation,
reference categories, category order, continuity corrections, exact versus
mid-p calculations, boundary behavior, and handling of invalid inputs.

Initially, R is the reference implementation because it contains the existing
methods and validation examples. Future features need not start in R. Choose
the implementation language according to the task, and implement the other
language against the agreed specification. An implementation's behavior alone
does not establish statistical correctness.

## Native interfaces

The Python package should work without an R installation. Use snake_case for
Python functions and arguments, for example `or_2x2_ci()`,
`cohen_kappa_test()`, `conf_level`, and `ci_variance`. These examples illustrate
naming conventions; individual signatures still require design review.
Preserve the current R interface until its separately planned naming migration.

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

1. Agree on the statistical specification, public interfaces, defaults, and
   reference examples before implementing a feature.
2. Implement first in the language most appropriate to the task, usually R
   during the initial development period.
3. Implement the corresponding feature in the other language using the same
   specification and reference cases.
4. Validate each implementation against published results, direct formulas,
   independently constructed small examples, or established software with
   matching assumptions. Retain the published SAS and Jung examples already
   used by the R package where applicable.
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
