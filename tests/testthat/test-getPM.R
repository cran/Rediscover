library(testthat)
library(Rediscover)

# Test getPM function

test_that("getPM works with basic example data", {
  # Load example data
  data("A_example")
  
  # Run getPM
  result <- getPM(A_example)
  
  # Check that result is not null
  expect_true(!is.null(result))
  
  # Check that result is a PMatrix object
  expect_true(is(result, "PMatrix"))
  
  # Check dimensions match input
  expect_equal(nrow(result), nrow(A_example))
  expect_equal(ncol(result), ncol(A_example))
  
  # Convert to matrix and check probabilities are in [0, 1]
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM works with Matrix class input", {
  data("A_Matrix")
  
  # Use smaller subset for faster testing
  A_small <- A_Matrix[1:100, 1:50]
  
  result <- getPM(A_small)
  
  # Check that result is PMatrix
  expect_true(is(result, "PMatrix"))
  
  # Check dimensions
  expect_equal(nrow(result), nrow(A_small))
  expect_equal(ncol(result), ncol(A_small))
  
  # Check probabilities are valid
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM works with standard matrix class", {
  # Create a simple binary matrix
  A <- matrix(c(0, 1, 0, 1, 1, 0, 0, 1, 1), nrow = 3, ncol = 3)
  
  result <- getPM(A)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), nrow(A))
  expect_equal(ncol(result), ncol(A))
  
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM handles NULL input", {
  expect_error(getPM(A = NULL), "not input matrix A")
})

test_that("getPM handles non-matrix input", {
  expect_error(getPM(A = data.frame(a = 1, b = 2)), 
               "input A must be a Matrix or a matrix class")
})

test_that("getPM handles empty matrix", {
  empty_matrix <- matrix(nrow = 0, ncol = 0)
  expect_error(getPM(A = empty_matrix), 
               "input A must have at least 1 row and 1 column")
})

test_that("getPM handles matrix with zero rows", {
  zero_row_matrix <- matrix(nrow = 0, ncol = 5)
  expect_error(getPM(A = zero_row_matrix), 
               "input A must have at least 1 row and 1 column")
})

test_that("getPM handles matrix with zero columns", {
  zero_col_matrix <- matrix(nrow = 5, ncol = 0)
  expect_error(getPM(A = zero_col_matrix), 
               "input A must have at least 1 row and 1 column")
})

test_that("getPM handles non-binary matrix", {
  non_binary <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 2, ncol = 3)
  expect_error(getPM(A = non_binary), "input A must be binary")
})

test_that("getPM returns PMatrix with correct structure", {
  data("A_example")
  
  result <- getPM(A_example)
  
  # Check that PMatrix has rowExps and colExps slots
  expect_true(length(result@rowExps) > 0)
  expect_true(length(result@colExps) > 0)
  
  # rowExps should have length equal to number of rows
  expect_equal(length(result@rowExps), nrow(A_example))
  
  # colExps should have length equal to number of columns
  expect_equal(length(result@colExps), ncol(A_example))
})

test_that("getPM works with all zeros matrix", {
  # Matrix with all zeros
  A_zeros <- matrix(0, nrow = 10, ncol = 10)
  
  result <- getPM(A_zeros)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), nrow(A_zeros))
  expect_equal(ncol(result), ncol(A_zeros))
  
  # Probabilities should be very low (close to 0) for all zeros
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1))
  expect_true(all(result_matrix < 0.5))  # Should be low probabilities
})

test_that("getPM works with all ones matrix", {
  # Matrix with all ones
  A_ones <- matrix(1, nrow = 10, ncol = 10)
  
  result <- getPM(A_ones)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), nrow(A_ones))
  expect_equal(ncol(result), ncol(A_ones))
  
  # Probabilities should be high (close to 1) for all ones
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1))
  expect_true(all(result_matrix > 0.5))  # Should be high probabilities
})

test_that("getPM works with sparse matrix pattern", {
  # Create a sparse binary matrix
  A_sparse <- matrix(0, nrow = 50, ncol = 50)
  A_sparse[sample(1:2500, 100)] <- 1  # Randomly set 100 entries to 1
  
  result <- getPM(A_sparse)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), nrow(A_sparse))
  expect_equal(ncol(result), ncol(A_sparse))
  
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM works with dense matrix pattern", {
  # Create a dense binary matrix (mostly 1s)
  A_dense <- matrix(1, nrow = 20, ncol = 20)
  A_dense[sample(1:400, 50)] <- 0  # Randomly set 50 entries to 0
  
  result <- getPM(A_dense)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), nrow(A_dense))
  expect_equal(ncol(result), ncol(A_dense))
  
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM works with small matrix", {
  # Very small matrix
  A_small <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
  
  result <- getPM(A_small)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM works with rectangular matrix (more rows)", {
  data("A_example")
  
  # More rows than columns
  A_rect <- A_example[1:200, 1:50]
  
  result <- getPM(A_rect)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), 200)
  expect_equal(ncol(result), 50)
  
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM works with rectangular matrix (more columns)", {
  data("A_example")
  
  # More columns than rows
  A_rect <- A_example[1:50, 1:200]
  
  result <- getPM(A_rect)
  
  expect_true(is(result, "PMatrix"))
  expect_equal(nrow(result), 50)
  expect_equal(ncol(result), 200)
  
  result_matrix <- as.matrix(result)
  expect_true(all(result_matrix >= 0 & result_matrix <= 1, na.rm = TRUE))
})

test_that("getPM can be used in downstream functions", {
  data("A_example")
  
  # Get probability matrix
  PM <- getPM(A_example)
  
  # Should be usable in getMutex without errors
  expect_silent(result <- getMutex(A = A_example, PM = PM))
  expect_true(!is.null(result))
})
