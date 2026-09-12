# Dataset Consistency Review

This note records the initial review of datasets in the `ibist` package against
the linked book manuscript.

## Summary

The package currently includes six teaching datasets:

- `depress`: Depression, insomnia, and stress dataset.
- `fluorosis`: Dental fluorosis grouped count data.
- `midbp`: Diastolic blood pressure and myocardial infarction rate data.
- `mihdl`: HDL cholesterol and myocardial infarction rate data.
- `nrs`: Non-restorative sleep cohort data.
- `rds`: Neonatal respiratory distress syndrome grouped count data.

## Findings

### `midbp`

This is the highest-priority mismatch.

The book Rnw files use `midbp.csv` with 16 rows and DBP categories 1--8. The
package dataset currently has 14 rows and DBP categories 1--7. The manuscript
also describes eight DBP categories:

- 1 = <70 mmHg
- 2 = 70--74 mmHg
- 3 = 75--79 mmHg
- 4 = 80--84 mmHg
- 5 = 85--89 mmHg
- 6 = 90--94 mmHg
- 7 = 95--99 mmHg
- 8 = >=100 mmHg

The package object appears to have combined the case counts from DBP categories
7 and 8 but retained only category 7 person-years. This changes rates and fitted
Poisson regression results, so the package should be corrected to match the
book's 16-row data.

### `fluorosis`

The package data match the book's `amoxfluor.csv` after standardizing column
names to lowercase. The package documentation says the counts sum to 579, but
the actual data sum to 308. The book statement that the `om == 0` subset contains
178 subjects is consistent with the package data.

Recommended action: correct the package documentation count from 579 to 308.

### `mihdl`

The package data match the book's `mihdl.csv`. The book's high HDL-C example
uses the values from HDL category 11:

- men: 9 cases over 1817.6 person-years
- women: 12 cases over 7496.3 person-years

The package documentation currently identifies HDL categories as integers 1--12.
If category cutpoints are available, documenting them would improve the teaching
connection.

### `depress`

The package data are consistent with the book's `depress.csv` plus one additional
derived variable, `depress1or2yr`. The manuscript has an editorial note that
`IEStert` may be redundant with `IESgrp`; the package documentation already
explains the relationship between the variables and the optional paper analytic
sample restriction.

### `nrs`

The package provides the six core variables used in the book. The old raw CSV has
additional derived variables such as age group and exercise/physical-activity
classification. The book derives these within Rnw chunks, which is reasonable
for teaching.

### `rds`

The package documentation still says the exact source is to be specified. This
should be resolved before CRAN submission if the original source can be
identified.

## Recommended Order

1. Fix `midbp` package data, documentation, and tests to match the book's 16-row
   data with DBP categories 1--8.
2. Fix the `fluorosis` documentation count.
3. Replace manuscript Rnw `read.csv("midbp.csv")` and `read.csv("amoxfluor.csv")`
   calls with `ibist::midbp` and `ibist::fluorosis` when updating the book.
4. Identify and document the original source for `rds`.
