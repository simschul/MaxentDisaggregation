test_that("rdirichlet_maxent returns a matrix with correct dimensions and column names", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  n <- 10
  result <- rdirichlet_maxent(n, shares)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_equal(colnames(result), names(shares))
})

test_that("rdirichlet_uniform returns a matrix with correct dimensions", {
  length <- 3
  n <- 10
  result <- rdirichlet_uniform(n, length, names = letters[1:length])
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length)
  expect_equal(colnames(result), letters[1:length])
})

test_that("rdirichlet_generalised returns a matrix with correct dimensions and column names", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  sds <- c(a = 0.1, b = 0.1, c = 0.1)
  n <- 10
  result <- rdirichlet_generalised(n, shares, sds)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_equal(colnames(result), names(shares))
})

test_that("rdirichlet returns a matrix with correct dimensions and column names", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  gamma <- 2
  n <- 10
  result <- rdirichlet(n, shares, gamma)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_equal(colnames(result), names(shares))
})

test_that("rdirichlet_standard returns a matrix with correct dimensions and column names", {
  alpha <- c(a = 1, b = 1, c = 1)
  n <- 10
  result <- rdirichlet_standard(n, alpha)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(alpha))
  expect_equal(colnames(result), names(alpha))
})

test_that("rshares returns a matrix with correct dimensions and sums to 1", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  sds <- c(a = 0.1, b = 0.1, c = 0.1)
  n <- 10
  result <- rshares(n, shares, sds)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_true(all(abs(rowSums(result) - 1) < 1e-6))
  expect_equal(colnames(result), names(shares))
})

test_that("rbeta3 returns a matrix with correct dimensions and sums <= 1", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  sds <- c(a = 0.05, b = 0.05, c = 0.05)
  n <- 10
  result <- rbeta3(n, shares, sds)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_true(all(rowSums(result) <= 1))
  expect_equal(colnames(result), names(shares))
})

test_that("rshares returns an error when beta is only provided for some elements", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  sds <- c(a = 0.05, b = 0.05, c = NA)
  n <- 10
  expect_error(rshares(n, shares, sds))

#
#   shares <- c(a = 0.2, b = 0.3, c = NA, d = NA)
#   sds <- c(a = 0.05, b = 0.05, c = NA, d = NA)
#   rshares(n, shares, sds)

})

# helper: reproducible draws
draw <- function(n, shares, sds, ...) {
  set.seed(42)
  rdirichlet_generalised(n, shares, sds, ...)
}

test_that("sample means match target shares within tolerance", {
  m <- draw(1e4,
            shares = c(.2, .3, .5),
            sds    = c(.02, .03, .04))
  expect_equal(colMeans(m),
               c(.2, .3, .5),
               tolerance = 0.01)          # 1 % relative tolerance
})

test_that("warning is issued when relative bias exceeds threshold", {
  # pump up SDs OR tighten max_rel_bias to trigger the warning
  expect_warning(
    draw(3e4,
         shares = c(.1, .3, .6),
         sds    = c(.1, .3, .6),
         max_rel_bias = 0.05),
    regexp = "relative bias"
  )
})

test_that("no warning when bias is within threshold", {
  expect_no_warning(
    draw(5e4,
         shares = c(.4, .6),
         sds    = c(.02, .02),            # small SDs
         max_rel_bias = 0.1)
  )
})

test_that("mismatched names in shares and sds throw an error", {
  expect_error(
    rdirichlet_generalised(10,
                           shares = c(A = .5, B = .5),
                           sds    = c(.05, .05)),          # unnamed sds
    "`shares` and `sds` must have the same names"
  )
})


