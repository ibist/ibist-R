# GLM goodness-of-fit specification and validation

`glm.gof(fit, method="pearson", g=10, asymcheck=FALSE)` and Python's
`glm_gof(fit, method="pearson", g=10, asymcheck=False)` take the fitted
model first. One method returns htest/TestResult. Multiple methods return
an ordered collection, preserving repetitions, with unadjusted p-values.

## Statistical scope

Pearson and deviance support binomial and Poisson GLMs with fixed
dispersion one. Aggregate identical design rows and offsets, using
integer trial/frequency weights. Retain expected totals and use J-rank df.
Zero-weight observations are excluded. No automatic pooling is performed.
The old `logistic.gof()` keeps its behavior and output without new warnings.

Hosmer-Lemeshow, Osius-Rojek, and Stukel require binomial logit models.
Hosmer-Lemeshow uses SAS target-frequency grouping with probability ties
at 14 significant digits and actual groups minus two df. Osius-Rojek
uses SAS's profile-count centering, J rather than J-rank, and reports
squared standardized Z with one df. Stukel is an efficient score test
with nuisance projection, not a likelihood-ratio refit. Its df reflects
the estimable increment (one or two). Formula details are in function help.

## Published numerical benchmarks

SAS LOGISTIC Example 79.19 reports these results for its training sample:

| Method | Published statistic | df | Published p-value |
| --- | --- | --- | --- |
| Pearson | 538.4877 | 544 | .5586 |
| Deviance | 606.5035 | 544 | .0324 |
| Hosmer-Lemeshow | 13.2967 | 8 | .1020 |
| Osius-Rojek | .8412 | 1 | .3590 |
| Stukel | 2.2581 | 2 | .3233 |

Both packages reproduce these to the reported precision. Source:
[SAS LOGISTIC manual][s], Example 79.19, Outputs 79.19.1 and 79.19.3.
The [SAS example program](
https://support.sas.com/documentation/onlinedoc/stat/ex_code/151/logiex19.html
) supplies the data definitions and training split.

Two additional published examples check grouped responses and offsets:
LOGISTIC Example 60.7 gives Pearson 6.6020 and deviance 7.7756, both
with four df. GENMOD's introductory insurance Poisson example gives
Pearson 2.8416 and deviance 2.8207, both with two df and a log-exposure
offset. Both packages reproduce these statistics. Tolerance is half a
unit in the last published decimal, with 1e-6 fitting allowance.

## Reproducible data and independent checks

`data-raw/export-python-glm-gof-reference.R` uses the public Wooldridge
Mroz data, reproduces SAS's RANUNI(3939) split with exact modular
arithmetic, and retains 552 training rows. Income is rounded to seven
significant digits as in the SAS listing. The generated CSV is shared by
both packages. Its response is nonparticipation (1-inlf), matching SAS's
default event. The CSV hash is recorded with eight R reference cases.
The exporter needs jsonlite and wooldridge; neither is a runtime dependency
for these GOF methods. Python fitting uses statsmodels (>=0.14).

Independent checks compare Poisson statistics with native GLM residuals
and deviance, Stukel with statsmodels' score test, and Osius-Rojek with
the direct information-matrix formula. Further checks cover offsets,
probit links, aliasing, repeated requests, frequency weights, zero-weight
boundary predictions, singular tests, and pre/attached/post diagnostics.
R/Python comparisons use 1e-12 fitting tolerances, 2e-7 relative and 1e-9
absolute tolerance for statistics/groups, and 1e-10 absolute for p-values.

These are all asymptotic tests. Published numerical agreement does not
establish finite-sample calibration. In particular, SAS warns against
interpreting its unreplicated Bernoulli Pearson/deviance p-values as valid
GOF evidence; our expected-count diagnostics flag these cases. Other R
packages can use different risk grouping or Osius-Rojek variants.
Thorough package checks remain deferred until before a push.

[s]: https://support.sas.com/documentation/onlinedoc/stat/151/logistic.pdf
