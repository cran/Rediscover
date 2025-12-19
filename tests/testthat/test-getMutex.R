library(testthat)
library(Rediscover)

# Test getMutex function

test_that("getMutex works with basic example data", {
  # Load example data
  data("A_example")
  
  # Compute probability matrix
  PMA <- getPM(A_example)
  
  # Run getMutex with default parameters
  result <- getMutex(A = A_example, PM = PMA)
  
  # Check that result is not null
  expect_true(!is.null(result))
  
  # Check that result is a matrix
  expect_true(is(result, "Matrix") || is(result, "matrix"))
  
  # Check dimensions (should be square matrix of genes)
  expect_equal(nrow(result), nrow(A_example))
  expect_equal(ncol(result), nrow(A_example))
  
  # Check that p-values are in valid range [0, 1]
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
  
  # Check that matrix is symmetric
  result_matrix <- as.matrix(result)
  expect_true(all(abs(result_matrix - t(result_matrix)) < 1e-10, na.rm = TRUE))
})

test_that("getMutex works with ShiftedBinomial method", {
  data("A_example")
  PMA <- getPM(A_example)
  
  result <- getMutex(A = A_example, PM = PMA, method = "ShiftedBinomial")
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("getMutex works with Binomial method", {
  data("A_example")
  PMA <- getPM(A_example)
  
  result <- getMutex(A = A_example, PM = PMA, method = "Binomial")
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("getMutex works with RefinedNormal method", {
  data("A_example")
  PMA <- getPM(A_example)
  
  result <- getMutex(A = A_example, PM = PMA, method = "RefinedNormal")
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})


test_that("getMutex handles NULL input matrix", {
  expect_error(getMutex(A = NULL), "not input matrix A")
})

test_that("getMutex handles non-matrix input", {
  expect_error(getMutex(A = data.frame(a = 1, b = 2)), 
               "input A must be a Matrix or a matrix class")
})

test_that("getMutex handles empty matrix", {
  empty_matrix <- matrix(nrow = 0, ncol = 0)
  expect_error(getMutex(A = empty_matrix), 
               "input A must have at least 1 row and 1 column")
})

test_that("getMutex handles non-binary matrix", {
  non_binary <- matrix(c(0, 1, 2, 3), nrow = 2)
  expect_error(getMutex(A = non_binary), "input A must be binary")
})

test_that("getMutex handles invalid method", {
  data("A_example")
  PMA <- getPM(A_example)
  
  expect_error(getMutex(A = A_example, PM = PMA, method = "InvalidMethod"),
               'method must be "Exact", "RefinedNormal", "Binomial", "ShiftedBinomial"')
})

test_that("getMutex works with Matrix class input", {
  data("A_Matrix")
  
  # Use smaller subset for faster testing
  A_small <- A_Matrix[1:50, 1:30]
  PMA_small <- getPM(A_small)
  
  result <- getMutex(A = A_small, PM = PMA_small)
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
  expect_equal(nrow(result), nrow(A_small))
})


test_that("getMutex handles mixed option", {
  data("A_example")
  PMA <- getPM(A_example)
  
  result_mixed_true <- getMutex(A = A_example, PM = PMA, mixed = TRUE)
  result_mixed_false <- getMutex(A = A_example, PM = PMA, mixed = FALSE)
  
  expect_true(!is.null(result_mixed_true))
  expect_true(!is.null(result_mixed_false))
})
