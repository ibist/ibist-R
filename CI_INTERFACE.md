# Common confidence interval interface

R defines the interface first; Python translates dotted names to snake_case.
This contract applies to CI functions. Hypothesis tests keep `alternative`.

## Arguments and defaults

Every public R CI function and every implemented Python counterpart accepts
`side` with exactly these values:

| Value | Interval |
| --- | --- |
| `"two.sided"` | Two-sided interval (default) |
| `"lower"` | Lower confidence bound; upper endpoint is the domain maximum |
| `"upper"` | Upper confidence bound; lower endpoint is the domain minimum |

An omitted `side`, R `NULL`, or Python `None` selects two-sided limits.
Using a null default distinguishes omission from an explicitly supplied
`"two.sided"`, allowing reliable conflict checks with the compatibility alias.
Both packages retain CI `alternative`: `"greater"` maps to `"lower"`,
`"less"` to `"upper"`, and `"two.sided"` maps to itself. Conflicting values
raise an error, including explicit `side="two.sided"` with a one-sided alias.

Confidence level defaults to 0.95 (`conf.level` / `conf_level`). One-sided
bounds use the full error probability `1 - confidence level`; central
intervals allocate half to each tail. This also applies below 50% confidence.
They are computed directly, without rescaling a two-sided confidence level.

Method requests preserve order and repetitions. A single method is selected
by default; an explicit vector/list requests multiple methods. Existing R
positional arguments retain their positions; use `side` by name.

## Function correspondence

| R | Python | Default method |
| --- | --- | --- |
| `prop.1s.ci` | `prop_1s_ci` | `wilson` |
| `prop.2s.ci` | `prop_2s_ci` | `newcombe` |
| `rr.2x2.ci` | `rr_2x2_ci` | `koopman` |
| `or.2x2.ci` | `or_2x2_ci` | `wald` |
| `rate.1s.ci` | `rate_1s_ci` | `exact` |
| `rate.2s.ci` | `rate_2s_ci` | `log` |
| `or.paired.ci` | `or_paired_ci` | `wald` |
| `prop.paired.ci` | `prop_paired_ci` | `score` |
| `or.common.ci` | `or_common_ci` | `exact` |

All implemented counterparts use the same method names and default
correction flags. Chan--Zhang exact intervals are implemented in both
languages, including the shared nuisance tolerance and enumeration limit.
Correction and mid-p options are separate from `side`.

One-sided BP intervals invert the conditional inclusive tail, matching
Cornfield; mid-p subtracts half the observed-table probability. One-sided
Blaker bounds use the corresponding exact binomial tail. Wang paired
bounds use the native one-sided construction. Parameter domains are `[0,1]`
for probabilities, `[-1,1]` for differences, and `[0,Inf]` for ratios/rates.
The existing R WH formula can give a negative limit: its review is separate
and recorded in [the review note](misc/rate-interval-review.md).

## Results and validation

`ibist_ci$side` and `CIResult.side` record the selected interval form.
Optional `null.value` / `null_value` reference values also match R; rate
CIs leave this field unset. Defaults are checked for all nine entry points.
Confidence matrices retain lower/upper columns and one row per requested
method (per input and method for vector inputs). Estimates remain at full
precision. Python DataFrame exports include `side` metadata. Compatibility
`alternative` fields may remain; consumers should use `side` for CI form.

Python uses scalar parameter estimates and separate sample metadata. R's
existing two-rate result additionally includes both group rates in its
estimate vector; Python exposes them through `sample_rates`. Vector input
indices are one-based in R and zero-based in Python. These existing native
result representations are documented rather than silently changed.

R tests check aliases, conflicts, full-tail Poisson equations, conditional
Poisson bounds, paired signed score inversion, and all nine CI entry points.
`data-raw/export-python-ci-interface-reference.R` generates 300 cases for
all shared methods at four confidence levels and three interval forms.
A further 545 cases cover paired, stratified, and Chan--Zhang methods.
Python checks argument defaults, result metadata, method order, and numerical
limits against those cases, alongside independent formula tests. R source
hashes identify the implementation used for each export. No R runtime is
needed in Python.

```r
rate.1s.ci(5, exposure = 10, side = "lower")
```

```python
rate_1s_ci(5, exposure=10, side="lower")
```
