# Tests for the logit-normal and probit-normal distribution functions
# (d/p/q/r) exported by the package.

test_that("plogitnorm is a valid CDF: in [0,1], monotone, correct limits", {
  x <- seq(0.001, 0.999, length.out = 200)
  for (par in list(c(0, 1), c(-0.5, 0.8), c(1.2, 1.5))) {
    p <- plogitnorm(x, par[1], par[2])
    expect_true(all(p >= 0 & p <= 1))
    expect_true(all(diff(p) > 0))                       # strictly increasing
    expect_lt(plogitnorm(1e-10, par[1], par[2]), 0.05)  # -> 0 at left edge
    expect_gt(plogitnorm(1 - 1e-10, par[1], par[2]), 0.95)
  }
})

test_that("pprobitnorm is a valid CDF: in [0,1], monotone, correct limits", {
  x <- seq(0.001, 0.999, length.out = 200)
  for (par in list(c(0, 1), c(-0.5, 0.8), c(1.2, 1.5))) {
    p <- pprobitnorm(x, par[1], par[2])
    expect_true(all(p >= 0 & p <= 1))
    expect_true(all(diff(p) > 0))
    expect_lt(pprobitnorm(1e-10, par[1], par[2]), 0.05)
    expect_gt(pprobitnorm(1 - 1e-10, par[1], par[2]), 0.95)
  }
})

test_that("quantile is the inverse of the CDF (logit-normal)", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (par in list(c(0, 1), c(-0.5, 0.8), c(1.2, 1.5))) {
    q <- qlogitnorm(probs, par[1], par[2])
    expect_equal(plogitnorm(q, par[1], par[2]), probs, tolerance = 1e-8)
  }
})

test_that("quantile is the inverse of the CDF (probit-normal)", {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (par in list(c(0, 1), c(-0.5, 0.8), c(1.2, 1.5))) {
    q <- qprobitnorm(probs, par[1], par[2])
    expect_equal(pprobitnorm(q, par[1], par[2]), probs, tolerance = 1e-8)
  }
})

test_that("density integrates to 1 (logit-normal and probit-normal)", {
  for (par in list(c(0, 1), c(-0.5, 0.8), c(1.2, 1.5))) {
    expect_equal(
      stats::integrate(dlogitnorm, 0, 1, mu = par[1], sigma = par[2])$value,
      1, tolerance = 1e-5
    )
    expect_equal(
      stats::integrate(dprobitnorm, 0, 1, mu = par[1], sigma = par[2])$value,
      1, tolerance = 1e-5
    )
  }
})

test_that("density equals the numerical derivative of the CDF", {
  x <- c(0.2, 0.4, 0.6, 0.8)
  h <- 1e-5
  for (par in list(c(0, 1), c(-0.5, 0.8))) {
    num <- (plogitnorm(x + h, par[1], par[2]) -
              plogitnorm(x - h, par[1], par[2])) / (2 * h)
    expect_equal(dlogitnorm(x, par[1], par[2]), num, tolerance = 1e-4)

    num <- (pprobitnorm(x + h, par[1], par[2]) -
              pprobitnorm(x - h, par[1], par[2])) / (2 * h)
    expect_equal(dprobitnorm(x, par[1], par[2]), num, tolerance = 1e-4)
  }
})

test_that("probit-normal mean has the known closed form pnorm(mu/sqrt(1+sigma^2))", {
  for (par in list(c(0, 1), c(0.4, 0.8), c(-0.7, 1.3))) {
    m_int <- stats::integrate(function(x) x * dprobitnorm(x, par[1], par[2]),
                              0, 1)$value
    expect_equal(m_int, stats::pnorm(par[1] / sqrt(1 + par[2]^2)),
                 tolerance = 1e-5)
  }
})

test_that("eprobitnorm matches numerical integration and is vectorised", {
  for (par in list(c(0, 1), c(0.4, 0.8), c(-0.7, 1.3))) {
    m_int <- stats::integrate(function(x) x * dprobitnorm(x, par[1], par[2]),
                              0, 1)$value
    expect_equal(eprobitnorm(par[1], par[2]), m_int, tolerance = 1e-5)
  }
  # symmetry: m(mu) = 1 - m(-mu)
  expect_equal(eprobitnorm(0, 1), 0.5, tolerance = 1e-12)
  expect_equal(eprobitnorm(0.6, 0.9), 1 - eprobitnorm(-0.6, 0.9), tolerance = 1e-12)
  # vectorised over mu and sigma, with recycling
  mu <- c(-1, 0, 1); sigma <- c(0.5, 1, 2)
  expect_equal(eprobitnorm(mu, sigma),
               mapply(eprobitnorm, mu, sigma), tolerance = 1e-12)
  expect_length(eprobitnorm(mu, 1), 3)
})

test_that("random generators match the theoretical CDF (seeded)", {
  set.seed(123)
  n <- 2e5

  x <- rlogitnorm(n, -0.4, 1.1)
  expect_true(all(x > 0 & x < 1))
  for (p in c(0.25, 0.5, 0.75)) {
    expect_equal(unname(stats::quantile(x, p)),
                 qlogitnorm(p, -0.4, 1.1), tolerance = 0.01)
  }

  x <- rprobitnorm(n, 0.4, 0.8)
  expect_true(all(x > 0 & x < 1))
  for (p in c(0.25, 0.5, 0.75)) {
    expect_equal(unname(stats::quantile(x, p)),
                 qprobitnorm(p, 0.4, 0.8), tolerance = 0.01)
  }
})
