# Tests for the main mcmap() interface: input handling, validation,
# type matching, and the structure of the returned object.

test_that("mcmap rejects targets that are not length-2 vectors", {
  expect_error(mcmap(0.5, "beta"))
  expect_error(mcmap(c(0.5), "beta"))
  expect_error(mcmap(c(0.3, 0.7, 0.9), "beta"))
})

test_that("mcmap rejects inadmissible (m, c) values", {
  expect_error(mcmap(c(0,    0.7), "beta"), "inadmissible")  # m <= 0
  expect_error(mcmap(c(1,    0.7), "beta"), "inadmissible")  # m >= 1
  expect_error(mcmap(c(0.3,  0.5), "beta"), "inadmissible")  # c <= 0.5
  expect_error(mcmap(c(0.3,  1.0), "beta"), "inadmissible")  # c >= 1
})

test_that("mcmap does partial matching of the type argument", {
  expect_equal(mcmap(c(0.3, 0.7), "b")$type, "beta")
  expect_equal(mcmap(c(0.3, 0.7), "logit")$type, "logitnorm")
  expect_equal(mcmap(c(0.3, 0.7), "probit")$type, "probitnorm")
  expect_error(mcmap(c(0.3, 0.7), "gamma"))  # not a known type
})

test_that("mcmap returns a well-formed mcmapper_output object", {
  res <- mcmap(c(0.3, 0.7), "beta")
  expect_s3_class(res, "mcmapper_output")
  expect_named(res, c("type", "value"))
  expect_equal(res$type, "beta")
  expect_length(res$value, 2)
  expect_true(all(is.finite(res$value)))
})

test_that("parameter names match the family", {
  expect_named(mcmap(c(0.3, 0.7), "beta")$value,       c("alpha", "beta"))
  expect_named(mcmap(c(0.3, 0.7), "logitnorm")$value,  c("mu", "sigma"))
  expect_named(mcmap(c(0.3, 0.7), "probitnorm")$value, c("mu", "sigma"))
})

test_that("named target is order-independent (uses 'm' and 'c' names)", {
  a <- mcmap(c(m = 0.3, c = 0.7), "beta")$value
  b <- mcmap(c(c = 0.7, m = 0.3), "beta")$value
  expect_equal(a, b, tolerance = 1e-8)
})

test_that("unnamed target is read as (m, c) in order", {
  named   <- mcmap(c(m = 0.3, c = 0.7), "beta")$value
  unnamed <- mcmap(c(0.3, 0.7), "beta")$value
  expect_equal(named, unnamed, tolerance = 1e-8)
})

test_that("extreme-but-admissible values emit an informational message", {
  # general 'Extreme' warning fires at |m-0.5|>0.499 or |c-0.75|>0.24
  expect_message(suppressWarnings(mcmap(c(0.3, 0.992), "probitnorm")), "Extreme")
  # logit-normal additionally announces the 1-D fallback at |c-0.75|>0.19
  expect_message(mcmap(c(0.3, 0.95), "logitnorm"), "1D")
})

test_that("named target with non-(m,c) names falls back to positional order", {
  # Regression: previously crashed with 'missing value where TRUE/FALSE needed'
  # because target['m'] was NA. Should now read positionally (mean, c-stat).
  ref <- mcmap(c(0.3, 0.7), "beta")$value
  expect_equal(mcmap(c(mean = 0.3, cstat = 0.7), "beta")$value, ref,
               tolerance = 1e-8)
  expect_equal(mcmap(c(0.3, 0.7), "beta")$value, ref, tolerance = 1e-8)
})
