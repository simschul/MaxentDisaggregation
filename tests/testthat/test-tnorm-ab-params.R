# tests/testthat/test-tnorm-ab-params.R
set.seed(10101)

# Helper: SE-based tolerance checks (robust across RNG/platforms)
.expect_moments_close <- function(x, m, s, k = 6) {
  n  <- length(x)
  sm <- mean(x)
  ss <- sd(x)
  se_mean <- s / sqrt(n)
  se_sd   <- s / sqrt(2 * (n - 1))
  expect_true(abs(sm - m) <= k * se_mean,
              info = sprintf("mean mismatch: got %.6f vs %.6f (tol=%.3g)", sm, m, k * se_mean))
  expect_true(abs(ss - s) <= k * se_sd,
              info = sprintf("sd mismatch: got %.6f vs %.6f (tol=%.3g)", ss, s, k * se_sd))
}

test_that("Untruncated case returns (mu, sigma) = (m, s)", {
  out <- tnorm_ab_params_from_moments(m = 5, s = 2, a = -Inf, b = Inf)
  expect_equal(out$mu, 5)
  expect_equal(out$sigma, 2)
  expect_identical(out$method, "closed")
})

test_that("Lower-truncated [a, Inf) reduces to shift of tnorm_params_from_moments()", {
  a   <- 1.5
  m   <- 2.3
  s   <- 0.8
  out <- tnorm_ab_params_from_moments(m = m, s = s, a = a, b = Inf)
  out0 <- tnorm_params_from_moments(m - a, s)
  expect_equal(out$sigma, out0$sigma, tolerance = 1e-10)
  expect_equal(out$mu, out0$mu + a, tolerance = 1e-10)
  expect_true(is.finite(out$alpha))
  expect_true(is.infinite(out$beta))
})

test_that("Upper-truncated (-Inf, b] reduces to reflected lower-truncation", {
  b   <- 0.7
  m   <- 0.2
  s   <- 0.5
  out <- tnorm_ab_params_from_moments(m = m, s = s, a = -Inf, b = b)
  # Check basic structure
  expect_true(is.infinite(out$alpha))
  expect_true(is.finite(out$beta))
  # Simulate and verify moments
  n <- 60000
  x <- truncnorm::rtruncnorm(n, a = -Inf, b = b, mean = out$mu, sd = out$sigma)
  .expect_moments_close(x, m, s)
})

test_that("Two-sided truncation [a,b] produces matching moments via simulation", {
  a <- 0
  b <- 3
  m <- 1.2
  s <- 0.6
  out <- tnorm_ab_params_from_moments(m = m, s = s, a = a, b = b)
  expect_true(is.finite(out$mu))
  expect_gt(out$sigma, 0)
  expect_lt(out$alpha, out$beta)

  n <- 8e4
  x <- truncnorm::rtruncnorm(n, a = a, b = b, mean = out$mu, sd = out$sigma)
  expect_true(all(x >= a & x <= b))
  .expect_moments_close(x, m, s)
})

test_that("Consistency with one-sided baseline when [0, Inf)", {
  m <- 1.0
  s <- 0.7
  out_ab <- tnorm_ab_params_from_moments(m = m, s = s, a = 0, b = Inf)
  out_0  <- tnorm_params_from_moments(m, s)
  expect_equal(out_ab$sigma, out_0$sigma, tolerance = 1e-10)
  expect_equal(out_ab$mu, out_0$mu, tolerance = 1e-10)
})

test_that("Edge: tight interval still solves and matches moments", {
  a <- 0.0
  b <- 0.02
  # Choose feasible moments within [a,b]: mean near mid; small sd
  m <- 0.011
  s <- 0.003
  out <- tnorm_ab_params_from_moments(m, s, a, b)
  expect_gt(out$sigma, 0)
  n <- 1e5
  x <- truncnorm::rtruncnorm(n, a = a, b = b, mean = out$mu, sd = out$sigma)
  expect_true(all(x >= a & x <= b))
  .expect_moments_close(x, m, s, k = 7) # slightly looser due to extreme truncation
})

test_that("Edge: far-tail bounds are handled stably", {
  a <- -5
  b <-  2
  m <- -1.0
  s <-  0.35
  out <- tnorm_ab_params_from_moments(m, s, a, b)
  expect_true(is.finite(out$mu))
  expect_gt(out$sigma, 0)
  n <- 7e4
  x <- truncnorm::rtruncnorm(n, a = a, b = b, mean = out$mu, sd = out$sigma)
  .expect_moments_close(x, m, s)
})

test_that("Input validation: bad arguments throw errors", {
  expect_error(tnorm_ab_params_from_moments(m = NA_real_, s = 0.5, a = 0, b = 1),
               "m must be finite")
  expect_error(tnorm_ab_params_from_moments(m = 0.2, s = -1, a = 0, b = 1),
               "strictly positive")
  expect_error(tnorm_ab_params_from_moments(m = 0.2, s = 0.5, a = 1, b = 1),
               "Require a < b")
  expect_error(tnorm_ab_params_from_moments(m = 0.2, s = 0.5, a = 2, b = 1.5),
               "Require a < b")
})

test_that("Infeasible targets trigger an error", {
  # Mean outside the interval is infeasible
  expect_error(tnorm_ab_params_from_moments(m = 5, s = 0.1, a = 0, b = 1))

  # Variance way too large for a very narrow interval
  a <- 0; b <- 0.01
  expect_error(tnorm_ab_params_from_moments(m = 0.005, s = 1, a = a, b = b))
})

test_that("Large-n simulation confirms accuracy across random scenarios", {
  set.seed(2024)
  cases <- list(
    list(a = 0,    b = Inf,  m = 1.1, s = 0.8),
    list(a = -Inf, b = 1.5,  m = 0.7, s = 0.5),
    list(a = -2.0, b = 2.5,  m = 0.2, s = 0.9),
    list(a = 0.5,  b = 3.0,  m = 1.8, s = 0.6)
  )
  for (cs in cases) {
    out <- tnorm_ab_params_from_moments(m = cs$m, s = cs$s, a = cs$a, b = cs$b)
    expect_gt(out$sigma, 0)
    n <- 6e4
    x <- truncnorm::rtruncnorm(n, a = cs$a, b = cs$b, mean = out$mu, sd = out$sigma)
    .expect_moments_close(x, cs$m, cs$s)
  }
})
