# Core correctness tests for the mappers.
#
# Strategy: for a target (m, c), run mcmap to recover the parameters, then
# independently recompute the (m, c) implied by those parameters using the
# integrate()-based oracle in helper-oracle.R. A correct mapper reproduces the
# target. This exercises the full pipeline (objective, integration, optimiser,
# parameter transforms, m > 0.5 symmetry handling) without assuming anything
# about the package's internals.

# (m, c) targets inside the region the package treats as reliable
# (|m - 0.5| < 0.49 and |c - 0.75| < 0.19, so the default logit-normal path is
# used rather than the 1-D fallback).
targets <- expand.grid(m = c(0.15, 0.3, 0.45, 0.6, 0.7),
                       c = c(0.62, 0.7, 0.78, 0.85))

for (type in c("beta", "logitnorm", "probitnorm")) {
  test_that(sprintf("round-trip recovers the target (m, c) for %s", type), {
    for (i in seq_len(nrow(targets))) {
      m <- targets$m[i]; cc <- targets$c[i]
      got <- suppressMessages(mc_roundtrip(m, cc, type))
      expect_equal(unname(got["m"]), m, tolerance = 5e-3,
                   info = sprintf("%s: target m=%g c=%g", type, m, cc))
      expect_equal(unname(got["c"]), cc, tolerance = 5e-3,
                   info = sprintf("%s: target m=%g c=%g", type, m, cc))
    }
  })
}

test_that("logit-normal recovers the README example parameters", {
  res <- suppressMessages(mcmap(c(m = 0.43, c = 0.77), "logitnorm"))
  expect_equal(unname(res$value), c(-0.3605462, 1.1836304), tolerance = 1e-3)
})

test_that("beta recovers known (alpha, beta) from their implied (m, c)", {
  truth <- c(2, 5)
  mc <- mc_from_cdf(function(x) stats::pbeta(x, truth[1], truth[2]))
  res <- suppressMessages(mcmap(unname(mc), "beta"))
  expect_equal(unname(res$value), truth, tolerance = 1e-2)
})

test_that("probit-normal recovers known (mu, sigma) from their implied (m, c)", {
  truth <- c(0.4, 0.8)
  mc <- mc_from_cdf(function(x) pprobitnorm(x, truth[1], truth[2]))
  res <- suppressMessages(mcmap(unname(mc), "probitnorm"))
  expect_equal(unname(res$value), truth, tolerance = 1e-2)
})

test_that("m and 1 - m give mirror-image parameters", {
  # logit-/probit-normal: mu flips sign, sigma unchanged.
  for (type in c("logitnorm", "probitnorm")) {
    lo <- suppressMessages(mcmap(c(m = 0.3, c = 0.72), type)$value)
    hi <- suppressMessages(mcmap(c(m = 0.7, c = 0.72), type)$value)
    expect_equal(hi[1], -lo[1], tolerance = 1e-3, info = type)  # mu
    expect_equal(hi[2],  lo[2], tolerance = 1e-3, info = type)  # sigma
  }
  # beta: alpha and beta swap.
  lo <- suppressMessages(mcmap(c(m = 0.3, c = 0.72), "beta")$value)
  hi <- suppressMessages(mcmap(c(m = 0.7, c = 0.72), "beta")$value)
  expect_equal(unname(hi[1]), unname(lo[2]), tolerance = 1e-3)
  expect_equal(unname(hi[2]), unname(lo[1]), tolerance = 1e-3)
})
