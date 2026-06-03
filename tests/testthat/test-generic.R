# Tests for mcmap_generic(), the family-agnostic mapper that takes an arbitrary
# two-parameter CDF.

test_that("mcmap_generic with pbeta matches the dedicated beta mapper", {
  for (tc in list(c(0.2, 0.7), c(0.3, 0.75), c(0.4, 0.8))) {
    generic <- mcmap_generic(tc, stats::pbeta)
    beta    <- suppressMessages(mcmap_beta(tc))
    expect_equal(generic, unname(beta), tolerance = 1e-2,
                 info = sprintf("m=%g c=%g", tc[1], tc[2]))
  }
})

test_that("mcmap_generic output reproduces the target (m, c)", {
  tc <- c(0.3, 0.75)
  par <- mcmap_generic(tc, stats::pbeta)
  got <- mc_from_cdf(function(x) stats::pbeta(x, par[1], par[2]))
  expect_equal(unname(got["m"]), tc[1], tolerance = 5e-3)
  expect_equal(unname(got["c"]), tc[2], tolerance = 5e-3)
})

test_that("mcmap_beta and mcmap_logitnorm reject unknown methods", {
  expect_error(mcmap_beta(c(0.3, 0.7), method = "nope"))
  expect_error(mcmap_logitnorm(c(0.3, 0.7), method = "nope"))
  expect_error(mcmap_probitnorm(c(0.3, 0.7), method = "nope"))
})
