set.seed(456)

test_that("rlnorm2 matches requested linear-scale mean & sd", {
  n <- 40000
  m <- 1.2
  s <- 0.5
  x <- rlnorm2(n, mean = m, sd = s)
  # Check linear-scale moments
  # (SE bounds as in helper above, re-implemented locally)
  se_mean <- s / sqrt(n)
  se_sd   <- s / sqrt(2 * (n - 1))
  expect_true(abs(mean(x) - m) < 6 * se_mean)
  expect_true(abs(sd(x)   - s) < 6 * se_sd)
})

test_that("rlnorm2 implements the correct (m2, sd2) transformation", {
  n <- 35000
  m <- 2.0
  s <- 0.8
  m2  <- log(m^2 / sqrt(s^2 + m^2))
  sd2 <- sqrt(log(1 + (s^2 / m^2)))
  x <- rlnorm2(n, mean = m, sd = s)
  # On the log scale, samples should be N(m2, sd2^2)
  lx <- log(x)
  expect_true(abs(mean(lx) - m2) < 6 * sd2 / sqrt(n))
  expect_true(abs(sd(lx)   - sd2) < 6 * sd2 / sqrt(2 * (n - 1)))
})
