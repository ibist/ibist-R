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

## Book-based approximation diagnostics

The 2026-11-09 audit corrects R diagnostics first. Python translation and
shared-reference updates are deferred until R stabilizes; Python currently
reports GLM approximation adequacy as unassessed. Elizabeth's editorial
review of the proposed cross-chapter binomial guideline is pending in the
book's Calibration section. Numerical tests and grouping algorithms are
unchanged. All methods remain asymptotic; published agreement does not
establish finite-sample calibration.

- Binomial Pearson/deviance use each covariate-pattern group's fitted
  `n_g*p_hat_g*(1-p_hat_g) >= 5`. HL uses the same screen within risk groups
  and actual `H >= 6`, after ties. Applying Chapter 11's binomial guideline
  operationalizes Chapter 16's qualitative large-group requirement; the
  numerical variance cutoff is not explicitly stated in Chapter 16.
- Group sizes, total binomial trials, pattern count, and G/n are retained.
  No ratio cutoff is invented for ``G much smaller than n.'' Unreplicated
  Bernoulli groups necessarily fail the binomial variance screen.
- Poisson Pearson/deviance use all fitted expected counts >=5 as a
  conservative Chapter 23 screen. Number below 2 and fraction below 5 are
  also reported, with no invented percentage for ``many.''
- Osius–Rojek/Stukel adequacy is unassessed (`condition.met = NA`). The book
  gives no criterion, and the reviewed papers establish no finite-sample
  cutoff. Degrees of freedom, standardizing variance, and estimable added
  terms are computability information, never
  approximation passes. Existing inference errors remain in place.
- A multiple-method summary cannot hide an unassessed method behind a
  passing one: any failure gives FALSE; otherwise any unassessed method
  gives NA; all passes give TRUE.

The book's NRS categorical-age example was checked using 90,122 observations,
32 covariate patterns, and model rank 8 (24 GOF df). Its minimum fitted
binomial variance is approximately 31.85, satisfying the proposed screen.
Focused regression cases include expected successes/failures 6/18 but
variance 4.5 (fail), the binomial threshold 5 (pass), five actual HL groups
despite ten requested (fail), and six adequate HL groups (pass). Poisson
small-count summaries use fitted means directly. Published SAS statistics
remain unchanged, and pre/attached/post diagnostics agree without running
test statistics during pre-checks.

Other R packages can use different risk grouping or Osius-Rojek variants.
Thorough package checks remain deferred until before a push.

## Reference audit, 2026-11-09

The supplied full PDFs of Osius and Rojek (1992) and Hosmer et al. (1997)
were reviewed, including the original limit conditions and variance
formulas. Stukel's original full paper remains unreviewed.

- Osius–Rojek, pp. 1146–1148: total sample size and number of independent
  multinomial groups increase; category and parameter counts stay fixed.
  Section 2.3 requires smooth probabilities, nonsingular information,
  consistency and an asymptotic expansion of the estimator, probability
  and derivative bounds, variance growth, and a group-size condition.
  These are conditions on sequences of designs, not finite-sample cutoffs.
  The probability bounds concern true probabilities; fitted extremes can
  flag concern but cannot verify the condition.
- The variance condition is bounded J/sigma^2 as J increases. Merely having
  positive fitted variance does not verify this. The size condition is
  increasing arithmetic mean group size OR QM/sqrt(n) tending to zero.
  Groups of size one satisfy the latter. No expected-count-five screen
  should be imposed on this test.
- Osius–Rojek, pp. 1150–1151: for binomial Pearson statistics,
  sigma^2 = 2*sum(1-1/N_j) + RSS from the specified weighted regression.
  This matches the implemented variance formula. With individual binary
  groups the first term vanishes; RSS can also vanish, including the
  constant model. A two-sided normal test is suggested for individual
  grouping; the implemented squared statistic with chi-square(1) tails
  gives that same p-value.
- Hosmer et al., p. 978, independently summarizes the binary-data
  (X^2-n)/sqrt(RSS) standardization. It then gives conditional moment
  corrections; p. 979 describes scaled chi-square calibration for better
  small-sample performance. These refinements are not implemented by our
  current Osius–Rojek method and must not be described as equivalent to it.
- Hosmer et al., p. 970, specifies Stukel's two added score-test terms.
  Its null simulations use n=100 and n=500, with 500 replications per
  design (pp. 971–973). These designs do not establish a universal
  sample-size threshold. The paper cautions about limited power and
  performance with more complicated models (pp. 976–977).

For now, approximation adequacy remains unassessed for these methods.
Group-size summaries, fitted probability extremes, and fitted J/sigma^2
could describe the Osius–Rojek conditions, but no numerical pass/fail
threshold is justified by these papers. Positive variance and estimable
score terms remain computability checks. No inference or API was changed
as part of this reference review.

[or]: https://doi.org/10.1080/01621459.1992.10476271
[st]: https://doi.org/10.1080/01621459.1988.10478613
[comparison]: https://pubmed.ncbi.nlm.nih.gov/9160492/

[s]: https://support.sas.com/documentation/onlinedoc/stat/151/logistic.pdf
