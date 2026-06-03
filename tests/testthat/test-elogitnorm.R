# Tests for the internal logit-normal mean E[expit(Z)], Z ~ N(mu, sigma).
#
# elogitnorm() is the numerical heart of the logit-normal mapper. It evaluates a
# theta-type series (elogitnorm_base), preconditioned by two exact identities
# that "tame" extreme mu before the series is summed:
#
#   symmetry: m(mu, sigma) = 1 - m(-mu, sigma)
#   tilt:     m(mu, sigma) = exp(mu + sigma^2/2) * m(-mu - sigma^2, sigma)
#
# We test (a) those identities directly -- an oracle-free check that pins the
# recursion against itself -- and (b) the result against high-order adaptive
# integration across a grid that includes the large-sigma corner.

# E[expit(mu + sigma * S)], S ~ N(0,1), by high-accuracy adaptive integration.
elogit_oracle <- function(mu, sigma) {
  stats::integrate(function(s) (1 / (1 + exp(-(mu + sigma * s)))) * stats::dnorm(s),
                   -Inf, Inf, rel.tol = 1e-12)$value
}

test_that("elogitnorm matches adaptive-integration oracle across (mu, sigma) grid", {
  grid <- expand.grid(mu = c(0, -0.7, -2, -5),
                      sigma = c(0.3, 1, 3, 6, 9))
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; sigma <- grid$sigma[i]
    expect_equal(elogitnorm(mu, sigma), elogit_oracle(mu, sigma),
                 tolerance = 1e-7,
                 info = sprintf("mu=%g sigma=%g", mu, sigma))
  }
})

test_that("symmetry identity holds: m(mu,s) = 1 - m(-mu,s)", {
  grid <- expand.grid(mu = c(0.2, 0.7, 2, 5), sigma = c(0.4, 1, 2.5, 6))
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; sigma <- grid$sigma[i]
    expect_equal(elogitnorm(mu, sigma), 1 - elogitnorm(-mu, sigma),
                 tolerance = 1e-9,
                 info = sprintf("mu=%g sigma=%g", mu, sigma))
  }
})

test_that("tilt identity holds: m(mu,s) = exp(mu+s^2/2) * m(-mu-s^2,s)", {
  grid <- expand.grid(mu = c(-0.7, -2, -5), sigma = c(0.4, 1, 2.5))
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; sigma <- grid$sigma[i]
    expect_equal(elogitnorm(mu, sigma),
                 exp(mu + sigma^2 / 2) * elogitnorm(-mu - sigma^2, sigma),
                 tolerance = 1e-7,
                 info = sprintf("mu=%g sigma=%g", mu, sigma))
  }
})

test_that("elogitnorm(0, sigma) = 0.5 and is monotone increasing in mu", {
  for (sigma in c(0.3, 1, 3)) {
    expect_equal(elogitnorm(0, sigma), 0.5, tolerance = 1e-12)
    vals <- vapply(c(-3, -1, 0, 1, 3), elogitnorm, numeric(1), sigma = sigma)
    expect_true(all(diff(vals) > 0))
  }
})

test_that("value is continuous across the steps = 0 / steps >= 1 boundary", {
  # steps = floor(-mu / sigma^2); cross the boundary at mu = -sigma^2.
  sigma <- 0.8
  eps <- 1e-6
  mu0 <- -sigma^2
  expect_equal(elogitnorm(mu0 + eps, sigma), elogitnorm(mu0 - eps, sigma),
               tolerance = 1e-5)
})

test_that("elogitnorm is vectorised over mu and sigma (with recycling)", {
  mu <- c(-2, -0.5, 0, 0.5, 2)
  sigma <- c(0.5, 1, 2, 1, 0.5)
  v <- elogitnorm(mu, sigma)
  expect_length(v, 5)
  expect_equal(v, mapply(elogitnorm, mu, sigma), tolerance = 1e-12)

  # sigma recycled against a longer mu
  expect_equal(elogitnorm(c(-1, 0, 1), 1),
               c(elogitnorm(-1, 1), elogitnorm(0, 1), elogitnorm(1, 1)),
               tolerance = 1e-12)
})

test_that("extreme inputs use the integrate() fallback without error", {
  # steps = floor(-mu/sigma^2) > 1000 routes through stats::integrate().
  mu <- -50; sigma <- 0.2
  expect_gt(floor(-mu / sigma^2), 1000)
  val <- elogitnorm(mu, sigma)
  expect_true(is.finite(val))
  expect_true(val >= 0 && val <= 1)
})
