# Lessons from the approximation-check work

Recorded 2026-11-10. Read before changing statistical diagnostics or
translating them between R and Python.

## Mistakes made

I implemented approximation checks before establishing what `asymcheck()`
was intended to assess. I treated positive degrees of freedom, positive
standardizing variance, and estimable score terms as evidence of adequate
asymptotic approximation. Those establish computability, not calibration.

I introduced criteria without first grounding them in the teaching book,
translated unsettled behavior to Python, and initially relied on abstracts
and secondary descriptions without clearly stating the limits of that
review. This created avoidable code, tests, documentation, and history
cleanup. I should have reported uncertainty before implementation.

## Establish the statistical basis first

- Start with the book's actual section, notation, assumptions, and examples.
  Present the proposed check and source before implementing an unsettled
  interpretation.
- Read the relevant original reference. Record page numbers, the method
  variant, and what was actually verified. An inaccessible paper remains
  unverified; an abstract is not a review of its full conditions.
- Verify that each citation supports the specific measure, formula, or
  approximation claim. A familiar textbook title is not sufficient evidence.
- Separate applicability, computability, approximation adequacy, and power.
  Numerical agreement with another implementation verifies neither coverage
  nor test size.
- Do not convert limiting assumptions into invented finite-sample cutoffs.
  Report unassessed adequacy explicitly when no defensible screen is settled.

For example, Osius–Rojek permits individual binary groups. Positive fitted
variance does not verify its variance-growth condition, and expected counts
of at least five are not its adequacy criterion. Hosmer et al.'s conditional
moment corrections must not be confused with the original standardization.
See [the reference audit](GLM_GOF_VALIDATION.md#reference-audit-2026-11-09).

## Control implementation and validation

- Stabilize statistical behavior, interfaces, output, documentation, and
  focused validation in R before translating to Python.
- Use representative tests that distinguish plausible errors, boundaries,
  and published results. Avoid large grids that merely repeat implementation.
- Run affected checks during development; reserve thorough checks for before
  pushing. Do not repeat successful checks after unrelated changes.
- Keep edited source and documentation within 80 characters per line.
- Preserve existing book teaching material unless its revision is authorized.

## Report and preserve repository state accurately

- Inspect the actual diff and history before claiming work is complete or
  removed. Distinguish local branch history from GitHub history.
- Identify execution-environment network failures precisely; do not imply
  that the user's browser or terminal has the same restriction.
- Check branch protection before proposing a history rewrite. Rewriting
  local commits does not remove their published counterparts.
- The decision for this incident is to retain the original history and
  record corrections in ordinary commits. Do not reopen history cleanup
  without a new explicit instruction.
