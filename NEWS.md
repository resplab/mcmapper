# mcmapper 0.0.12

## Bug fixes

* `mcmap()` no longer fails with a cryptic error ("missing value where
  TRUE/FALSE needed") when the `target` is a named vector whose names are not
  exactly `"m"` and `"c"`. Named targets are now used only when both `"m"` and
  `"c"` are present; otherwise the documented positional reading
  `(mean, c-statistic)` is used.

* `mcmap(type = "logitnorm")` and `mcmap_logitnorm()` (default method) were
  unreliable when the mean was close to 0.5. Nelder-Mead's initial simplex
  scales with `logit(m)`, which is then near 0, so the optimiser either failed
  to converge and `NULL` was returned without a warning, or stopped at its
  starting point while reporting convergence, giving parameters whose
  c-statistic could be off by more than 0.1. A fit is now accepted only if it
  reproduces the target `(mean, c-statistic)` to within 1e-3; otherwise the
  optimiser restarts from `(0, 0)`, and a warning is issued if no fit is found.
  Results that were already correct are unchanged.

## New features

* `elogitnorm()` and `eprobitnorm()` are now exported. They return the mean
  (expected value) of the logit-normal and probit-normal distributions,
  respectively, and are vectorised over `mu` and `sigma`.

## Internal

* Added a `testthat` test suite covering the distribution functions, the
  internal logit-normal mean (including its symmetry and exponential-tilt
  identities), input validation, and round-trip recovery of `(mean, c-statistic)`
  targets for all three families.
