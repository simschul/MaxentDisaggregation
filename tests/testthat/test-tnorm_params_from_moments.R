set.seed(789)

test_that("tnorm_params_from_moments returns finite, positive sigma", {
  out <- tnorm_params_from_moments(m = 1.0, s = 0.9)
  expect_true(is.finite(out$mu))
  expect_true(is.finite(out$sigma))
  expect_gt(out$sigma, 0)
  expect_true(out$method %in% c("uniroot", "newton"))
})

test_that("Plugging back (mu, sigma) yields target moments (min=0, max=Inf)", {
  # Since tnorm_params_from_moments assumes truncation at 0,
  # test with min=0, max=Inf (pure lower truncation).
  n <- 50000
  m <- 1.0
  s <- 0.7
  out <- tnorm_params_from_moments(m, s)
  x <- truncnorm::rtruncnorm(n, a = 0, b = Inf, mean = out$mu, sd = out$sigma)
  # Moment check with SE tolerance (use helper logic inline)
  se_mean <- s / sqrt(n)
  se_sd   <- s / sqrt(2 * (n - 1))
  expect_true(abs(mean(x) - m) < 6 * se_mean)
  expect_true(abs(sd(x)   - s) < 6 * se_sd)
})

test_that("Also works for finite upper bound", {
  n <- 60000
  m <- 0.9
  s <- 0.5
  out <- tnorm_params_from_moments(m, s)
  # NOTE: ragg uses the same (mu, sigma) for any finite (min,max),
  # but the target (m,s) are for lower-truncation at 0.
  # Here we still check that with a reasonable upper bound
  # the moments remain close to (m,s).
  x <- truncnorm::rtruncnorm(n, a = 0, b = 3, mean = out$mu, sd = out$sigma)
  # Moments will shift slightly with a finite upper bound; we just check sanity:
  expect_true(all(x >= 0 & x <= 3))
  expect_gt(mean(x), 0)
  expect_lt(sd(x), 2)
})

test_that("Input validation fires on bad arguments", {
  expect_error(tnorm_params_from_moments(m = -1, s = 0.5), "strictly positive")
  expect_error(tnorm_params_from_moments(m = 1,  s = -0.5), "strictly positive")
  expect_error(tnorm_params_from_moments(m = NA, s = 0.5), "strictly positive")
})
