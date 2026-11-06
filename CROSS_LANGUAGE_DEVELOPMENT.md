# R and Python development principles

Agreed: 2026-11-06. Last updated: 2026-11-06.

This document records the agreed approach to a Python sister package for
ibist. It is the authoritative development policy for both implementations.
The Python repository is [ibist-py](https://github.com/ibist/ibist-py).
Its initial package structure uses `ibist` for both the distribution and import
names; methods and datasets remain to be implemented. These names
describe the development package and do not imply a reserved PyPI name.

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

Use NumPy and SciPy for numerical building blocks. Add dependencies when the
functionality needs them; do not require a statistical modeling framework
merely to implement table-based inference. Follow the same statistical
definitions across languages even when underlying libraries have different
defaults.

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
CSV alone is insufficient to preserve all of this information. Select the
export format and metadata schema when implementing the first dataset.

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
portable data and metadata formats, reference-case storage, release
coordination, and the eventual location of shared source material.
