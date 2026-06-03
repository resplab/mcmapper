# mcmapper 0.0.12

## Bug fixes

* `mcmap()` no longer fails with a cryptic error ("missing value where
  TRUE/FALSE needed") when the `target` is a named vector whose names are not
  exactly `"m"` and `"c"`. Named targets are now used only when both `"m"` and
  `"c"` are present; otherwise the documented positional reading
  `(mean, c-statistic)` is used.

## New features

* `elogitnorm()` and `eprobitnorm()` are now exported. They return the mean
  (expected value) of the logit-normal and probit-normal distributions,
  respectively, and are vectorised over `mu` and `sigma`.

## Internal

* Added a `testthat` test suite covering the distribution functions, the
  internal logit-normal mean (including its symmetry and exponential-tilt
  identities), input validation, and round-trip recovery of `(mean, c-statistic)`
  targets for all three families.
