set.seed(123)

# small helper: check sample mean/sd are close to targets using SE-based bounds
expect_moments_close <- function(x, m, s, k = 6) {
  n  <- length(x)
  sm <- mean(x)
  ss <- sd(x)
  se_mean <- s / sqrt(n)
  # For SD, use delta-method-ish tolerance: ~ s / sqrt(2*(n-1))
  se_sd <- s / sqrt(2 * (n - 1))
  expect_true(abs(sm - m) <= k * se_mean,
              info = sprintf("mean mismatch: got %.4f vs %.4f (tol=%.4g)", sm, m, k * se_mean))
  expect_true(abs(ss - s) <= k * se_sd,
              info = sprintf("sd mismatch: got %.4f vs %.4f (tol=%.4g)", ss, s, k * se_sd))
}

test_that("Case 1: Normal distribution branch", {
  n <- 20000
  m <- 5
  s <- 2
  x <- ragg(n, mean = m, sd = s, min = -Inf, max = Inf, log = FALSE)
  expect_length(x, n)
  expect_true(all(is.finite(x)))
  expect_moments_close(x, m, s)
})

test_that("Case 2: Lognormal via linear-scale mean/sd with log=TRUE and min<=0", {
  n <- 30000
  m <- 1
  s <- 0.3
  x <- ragg(n, mean = m, sd = s, min = 0, max = Inf, log = TRUE)
  expect_length(x, n)
  expect_true(all(x >= 0))
  # For a lognormal parameterized by mean/sd on linear scale,
  # sample moments should match (within SE).
  # True variance on linear scale:
  expect_moments_close(x, m, s)
})

test_that("Case 3: Truncated Normal branch (finite bounds) hits and respects bounds", {
  n <- 30000
  m <- 3
  s <- 1
  a <- 0    # important: tnorm_params_from_moments assumes lower truncation at 0
  b <- 6
  x <- ragg(n, mean = m, sd = s, min = a, max = b)
  expect_true(all(x >= a & x <= b))
  # Moments should be close to targets (SE-based). For truncated distributions,
  # sampling error can be a bit larger; keep default k=6 which is generous.
  expect_moments_close(x, m, s)
})

test_that("Case 4: Exponential when sd is NULL and min==0, max==Inf", {
  n <- 40000
  mu <- 2  # mean
  x <- ragg(n, mean = mu, sd = NULL, min = 0, max = Inf)
  expect_true(all(x >= 0))
  # For Exp(rate = 1/mu): mean = mu, sd = mu
  expect_moments_close(x, mu, mu)
})

test_that("Case 5: Uniform when mean/sd are NULL and bounds finite", {
  n <- 20000
  a <- -2
  b <- 5
  x <- ragg(n, min = a, max = b)
  expect_true(all(x >= a & x <= b))
  m_true <- (a + b) / 2
  s_true <- (b - a) / sqrt(12)
  expect_moments_close(x, m_true, s_true)
})

test_that("Case 6: Standard lognormal via meanlog/sdlog", {
  n <- 25000
  ml <- 0.3
  sl <- 0.7
  x <- ragg(n, meanlog = ml, sdlog = sl)
  expect_true(all(x > 0))
  # Check log-scale moments (because rlnorm is exact if we look at log(x))
  lx <- log(x)
  expect_true(abs(mean(lx) - ml) < 6 * sl / sqrt(n))
  expect_true(abs(sd(lx)   - sl) < 6 * sl / sqrt(2 * (n - 1)))
})

test_that("Error when no signature matches", {
  # mean given, sd NULL, but min != 0 ⇒ not exponential branch
  expect_error(ragg(10, mean = 2, sd = NULL, min = 1, max = Inf),
               "No matching distribution case")
  # mean/sd both NULL and one bound infinite ⇒ not uniform
  expect_error(ragg(10, mean = NULL, sd = NULL, min = -Inf, max = 1),
               "No matching distribution case")
})
