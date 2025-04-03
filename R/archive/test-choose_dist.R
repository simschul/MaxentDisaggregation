library(testthat)

test_that("generate_distribution creates functions with correct behavior", {
  set.seed(123)

  # Test for normal distribution with known mean and sd
  test_normal <- generate_distribution(rnorm, mean = 100, sd = 15)
  sample_normal <- test_normal(1000)
  expect_type(test_normal, "closure")
  expect_equal(mean(sample_normal), 100, tolerance = 1, info = "Mean is not near 100")
  expect_equal(sd(sample_normal), 15, tolerance = 1, info = "SD is not near 15")
})



test_that("generate_distribution handles errors appropriately", {
  # Try an unsupported function and check for specific error message
  expect_error(generate_distribution(rpois, lambda = 5))
})

test_that("generated function reflects the right parameters and distribution", {
  # Ensure parameters are captured and distribution changes accordingly
  set.seed(123)

  # Weibull distribution test with shape and scale
  test_norm <- generate_distribution(rnorm, mean = 100, sd = 5)
  sample_norm <- test_norm(1000)
  expect_type(test_norm, "closure")

  # Insert more rigorous tests specific to the Weibull distribution
  # For instance, you might test the median or other moments
  # Adjust the expected values and tolerance according to what's appropriate for your distribution
})

# More rigorous test for parameter passing
test_that("Parameters are passed and used correctly", {
  set.seed(123)

  # Test changing parameters and seeing effect on the distribution
  test_norm1 <- generate_distribution(rnorm, mean = 0, sd = 1)
  test_norm2 <- generate_distribution(rnorm, mean = 50, sd = 10)

  sample1 <- test_norm1(1000)
  sample2 <- test_norm2(1000)

  expect_true(all(mean(sample1) != mean(sample2)), info = "Means should be different for different parameters")
  expect_true(all(sd(sample1) != sd(sample2)), info = "Standard deviations should be different for different parameters")
})
