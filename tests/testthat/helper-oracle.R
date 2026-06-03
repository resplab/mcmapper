# Independent "oracle" utilities shared across test files.
#
# The package maps a (mean, c-statistic) pair to distribution parameters by
# matching two moments of the risk distribution's CDF F on [0, 1]:
#
#   F1 = integral_0^1 F(x) dx           = 1 - m
#   F2 = integral_0^1 F(x)^2 dx         = 1 - (2 c m - (2c - 1) m^2)
#
# The second identity rearranges to a closed form for the c-statistic:
#
#   c = (1 - m^2 - F2) / (2 m (1 - m))
#
# We recompute (m, c) here directly from a CDF using high-order adaptive
# integration -- deliberately NOT Gauss-Hermite, which loses accuracy for the
# near-step integrands that arise at large sigma. This gives a ground truth that
# is independent of the package's own optimisation/integration code, so a
# round-trip (params -> (m, c) -> mcmap -> params) genuinely tests the mapper.

# Mean and c-statistic implied by a CDF on [0, 1].
mc_from_cdf <- function(cdf) {
  F1 <- stats::integrate(cdf, 0, 1, rel.tol = 1e-10)$value
  F2 <- stats::integrate(function(x) cdf(x)^2, 0, 1, rel.tol = 1e-10)$value
  m  <- 1 - F1
  cc <- (1 - m^2 - F2) / (2 * m * (1 - m))
  c(m = m, c = cc)
}

# CDF closures for each supported family, given a parameter vector.
cdf_for <- function(type, par) {
  switch(
    type,
    "beta"       = function(x) stats::pbeta(x, par[1], par[2]),
    "logitnorm"  = function(x) plogitnorm(x, par[1], par[2]),
    "probitnorm" = function(x) pprobitnorm(x, par[1], par[2]),
    stop("unknown type")
  )
}

# Given a target (m, c) and a type, run mcmap and report the (m, c) actually
# implied by the recovered parameters. A correct mapper returns ~ the target.
mc_roundtrip <- function(m, c, type) {
  res <- mcmap(c(m = m, c = c), type)
  mc_from_cdf(cdf_for(type, res$value))
}
