# R Package **ibist**

The R package **ibist** is a companion to the book
"Introduction to Biostatistical Analysis of Proportions and Rates"
by Elizabeth D. Schifano and Jun Yan

## Features

The package **ibist** provides functions, datasets, and demos
used in the book.

`indep.test(x, method = "pearson")` wraps `stats::chisq.test()` for
contingency-table independence, with continuity correction off by default.
Use `method = "lr"` for the likelihood-ratio G-squared test. Both return
observed and expected cell counts alongside the test results. Set
`correct = TRUE` for Yates' correction in Pearson's test of a 2 x 2 table.

`nominal.assoc(x)` returns Phi, Pearson's contingency coefficient, and
Cramer's V as descriptive point estimates. Use `measure = "cramer-v"` to
select one measure. Phi is non-negative for all table sizes and can exceed
one for larger tables; Cramer's V is normalized between zero and one.

## Development

The latest version of the package is under development at
[GitHub](https://github.com/ibist/ibist-R). If it is able to pass
the automated package checks, one may install it by

``` r
if (! require(remotes)) install.packages("remotes")
remotes::install_github("ibist/ibist-R", upgrade = "never")
```

The [Python sister package](https://github.com/ibist/ibist-py) follows the
[R and Python development principles](CROSS_LANGUAGE_DEVELOPMENT.md),
including shared statistical specifications, data preparation, and validation.
After regenerating R datasets, run `Rscript data-raw/export-python.R` from
this repository to refresh the portable Python data and metadata. The exporter
reads the prepared R objects; it requires jsonlite during development.

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)

