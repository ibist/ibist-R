# Published validation of Poisson rate power

Checked on 2026-11-08 in both R and Python. Expected values below are
transcribed from Minitab's published examples, independently of our code.

## One sample

Source: Minitab, LLC, *Example of Power and Sample Size for 1-Sample
Poisson Rate*, [Minitab Support][minitab]. Navigate through Power and
Sample Size, Hypothesis Tests, 1-Sample Poisson Rate, Example.

Use null rate 15, alternative rate 13, exposure per subject 1, significance
level .05, `alternative="less"`, and no continuity correction.

| n | Published power | ibist power (both languages) |
| --- | --- | --- |
| 25 | 0.842947 | 0.842947241487 |
| 30 | 0.898200 | 0.898200031166 |

## Two samples

Source: Minitab, LLC, *Example of Power and Sample Size for 2-Sample
Poisson Rate*, [Minitab Support][minitab]. Follow the same path for
2-Sample Poisson Rate.

Use baseline rate 35, equal sample sizes, exposure per subject 1 in each
group, significance level .05, a two-sided alternative, no correction,
and `strict=TRUE` / `strict=True` to include both rejection tails.

| Comparison rate | n per group | Published power | ibist power |
| --- | --- | --- | --- |
| 32 | 79 | 0.902793 | 0.902793409906 |
| 38 | 86 | 0.902550 | 0.902550388521 |

For target power .9, ibist returns continuous n of 78.221900146 and
85.226846158. Rounding up reproduces Minitab's 79 and 86 per group.

## Interpretation and scope

All four powers match within half a unit of the last published decimal
(`5e-7` absolute tolerance). Focused regression tests in both packages
retain these published values and check the rounded sample sizes.

Minitab's methods-and-formulas pages for these procedures use the same
normal score boundaries and alternative variances. Our default
`strict=FALSE` retains the dominant-tail convention of the proportion
power functions; explicitly request both tails for the two-sided comparison.
No numerical implementation changes were needed.

These published examples validate the uncorrected normal methods with
equal allocation and unit exposure. Unequal allocation, follow-up units,
continuity correction, and other inverses have separate formula and
invariance checks; they are not covered by these four published results.
This does not validate exact Poisson power or guarantee finite-sample size.

[minitab]: https://support.minitab.com
