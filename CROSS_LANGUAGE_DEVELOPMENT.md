# R and Python development principles

Agreed: 2026-11-06. Last updated: 2026-11-06.

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
families are being implemented in R first; Python parity remains pending.

### Contingency-table independence tests

On 2026-11-06, the authors chose `indep.test()` in R and `indep_test()` in
Python to provide a common interface to independence tests. The first R
implementation supports `method = "pearson"` (default) and `method = "lr"`;
Python parity remains pending. This prepares replacements for the book's
use of DescTools for the likelihood-ratio test and expected cell counts;
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

Fisher's exact test, simulated p-values, goodness-of-fit tests, and
likelihood-ratio continuity corrections are outside this first interface.

Direct development comparisons with DescTools 0.99.60 covered both book
G-test examples, a zero-cell table, and 100 seeded random tables of sizes
2 x 2 through 5 x 5. Statistics and expected counts agreed within `1e-10`.
For the tumor table (rows `109,68,14; 72,79,37; 9,42,139`), DescTools returned
a zero p-value; the upper-tail calculation in `indep.test()` retains the
probability of approximately `5.73e-53`. Compute survival probabilities
directly, not by subtracting the lower-tail probability from one. DescTools
is not a runtime or test dependency; retain independent formula and Poisson
model checks in the test suite.

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

Most existing statistical code is in R and requires translation. The Wang
paired-interval implementation uses Rcpp. Sharing its C++ computational core
would require separating that core from R-specific types and error handling,
then providing language-specific bindings. Make this change when its benefit
justifies the additional build and maintenance work.

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
