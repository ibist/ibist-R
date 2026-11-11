# Draft replacements for book/Rnw/pepi-vf.Rnw and pepi-coffee.Rnw.
# Run with Rscript using an installed development version of ibist.
# Tables have exposure in rows and outcome in columns, positive first.
# The displays cover the intervals and tests interpreted in the examples.
# Attributable fractions, population attributable risk, and E-values are
# omitted from this draft. Retaining them would need separate exposition.

ci_rows <- function(fit) {
  cbind(estimate = unname(fit$estimate), fit$conf.int)
}

association_tests <- function(tab) {
  pearson <- ibist::indep.test(tab, method = "pearson", correct = FALSE)
  fisher <- ibist::indep.test(tab, method = "fisher")
  c(pearson.chisq = unname(pearson$statistic),
    pearson.p.value = pearson$p.value, fisher.p.value = fisher$p.value)
}

# Cohort example: risk of fracture hospitalization by baseline fracture.
vf <- matrix(c(38, 62, 267, 717), nrow = 2, byrow = TRUE)
vf.counts <- rbind(vf, colSums(vf))
vf.risk <- ibist::prop.1s.ci(vf.counts[, 1], rowSums(vf.counts),
                            method = "clopper-pearson")
vf.risk.per100 <- 100 * ci_rows(vf.risk)
rownames(vf.risk.per100) <- c("Exposed", "Unexposed", "Total")

vf.rr <- ibist::rr.2x2.ci(vf, method = "katz")
vf.or <- ibist::or.2x2.ci(vf, method = "wald")
vf.rd <- ibist::prop.2s.ci(vf, method = "wald", correct = FALSE)
vf.effects <- rbind(
  "Risk ratio" = ci_rows(vf.rr)[1, ],
  "Odds ratio" = ci_rows(vf.or)[1, ],
  "Risk difference per 100" = 100 * ci_rows(vf.rd)[1, ]
)
vf.tests <- association_tests(vf)

print(vf.risk.per100)
print(vf.effects)
print(vf.tests)

# Case-control example: coffee exposure and oral cleft.
coffee <- matrix(c(265, 483, 112, 280), nrow = 2, byrow = TRUE)
coffee.counts <- rbind(coffee, colSums(coffee))
# These row odds describe the sampled outcomes, given the exposure group.
coffee.sample.prop <- ibist::prop.1s.ci(
  coffee.counts[, 1], rowSums(coffee.counts), method = "clopper-pearson"
)
coffee.sample.odds <- ci_rows(coffee.sample.prop)
coffee.sample.odds <- coffee.sample.odds / (1 - coffee.sample.odds)
rownames(coffee.sample.odds) <- c("Exposed", "Unexposed", "Total")

coffee.or <- ibist::or.2x2.ci(coffee, method = "wald")
coffee.effect <- ci_rows(coffee.or)
rownames(coffee.effect) <- "Exposure odds ratio"
coffee.tests <- association_tests(coffee)

print(coffee.sample.odds)
print(coffee.effect)
print(coffee.tests)
