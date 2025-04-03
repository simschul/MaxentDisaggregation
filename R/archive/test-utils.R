# Test for a basic matrix without column names
test_that("melt_matrix with basic matrix", {
  mat <- matrix(1:4, ncol = 2)
  result <- melt_matrix(mat)

  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 4)
  expect_true(all(c("x", "y", "value") %in% names(result)))
  expect_equal(result$x, c(1, 2, 1, 2)) # Update this line according to actual behavior
  expect_equal(as.integer(result$y), c(1, 1, 2, 2)) # Assuming y is factor, convert to integer
  expect_equal(result$value, c(1, 2, 3, 4))
})

# Test for removing non-finite values
test_that("melt_matrix with remove_non_finite = TRUE", {
  mat <- matrix(c(NA, 1, Inf, 2), ncol = 2)
  result <- melt_matrix(mat, remove_non_finite = TRUE)

  expect_equal(nrow(result), 2) # Update this line according to actual behavior
  expect_false(any(!is.finite(result$value)))
})
