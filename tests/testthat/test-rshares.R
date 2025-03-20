test_that("rdir_maxent returns a matrix with correct dimensions and column names", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  n <- 10
  result <- rdir_maxent(n, shares)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_equal(colnames(result), names(shares))
})

test_that("rdir1 returns a matrix with correct dimensions", {
  length <- 3
  n <- 10
  result <- rdir1(n, length, names = letters[1:length])
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length)
  expect_equal(colnames(result), letters[1:length])
})

test_that("rdirg returns a matrix with correct dimensions and column names", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  sds <- c(a = 0.1, b = 0.1, c = 0.1)
  n <- 10
  result <- rdirg(n, shares, sds)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_equal(colnames(result), names(shares))
})

test_that("rdir returns a matrix with correct dimensions and column names", {
  shares <- c(a = 0.2, b = 0.3, c = 0.5)
  gamma <- 2
  n <- 10
  result <- rdir(n, shares, gamma)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), length(shares))
  expect_equal(colnames(result), names(shares))
})

test_that("rdirichlet returns a matrix with correct dimensions and column names", {
  alpha <- c(a = 1, b = 1, c = 1)
  n <- 10
  result <- rdirichlet(n, alpha)
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

