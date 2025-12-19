library(testthat)
library(Rediscover)

# Test getMutexAB function

test_that("getMutexAB works with basic example data", {
  # Load example data
  data("A_example")
  data("B_example")
  
  # Compute probability matrices
  PMA <- getPM(A_example)
  PMB <- getPM(B_example)
  
  # Run getMutexAB with default parameters
  result <- getMutexAB(A = A_example, PMA = PMA, B = B_example, PMB = PMB)
  
  # Check that result is not null
  expect_true(!is.null(result))
  
  # Check that result is a matrix
  expect_true(is(result, "Matrix") || is(result, "matrix"))
  
  # Check dimensions (should be nrow(A) x nrow(B))
  expect_equal(nrow(result), nrow(A_example))
  expect_equal(ncol(result), nrow(B_example))
  
  # Check that p-values are in valid range [0, 1]
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("getMutexAB works with ShiftedBinomial method", {
  data("A_example")
  data("B_example")
  
  PMA <- getPM(A_example)
  PMB <- getPM(B_example)
  
  result <- getMutexAB(A = A_example, PMA = PMA, B = B_example, PMB = PMB, 
                       method = "ShiftedBinomial")
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("getMutexAB works with Binomial method", {
  data("A_example")
  data("B_example")
  
  PMA <- getPM(A_example)
  PMB <- getPM(B_example)
  
  result <- getMutexAB(A = A_example, PMA = PMA, B = B_example, PMB = PMB, 
                       method = "Binomial")
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("getMutexAB works with RefinedNormal method", {
  data("A_example")
  data("B_example")
  
  PMA <- getPM(A_example)
  PMB <- getPM(B_example)
  
  result <- getMutexAB(A = A_example, PMA = PMA, B = B_example, PMB = PMB, 
                       method = "RefinedNormal")
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("getMutexAB works with lower.tail = FALSE (co-occurrence)", {
  data("A_example")
  data("B_example")
  
  PMA <- getPM(A_example)
  PMB <- getPM(B_example)
  
  result <- getMutexAB(A = A_example, PMA = PMA, B = B_example, PMB = PMB, 
                       lower.tail = FALSE)
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("getMutexAB handles NULL input matrix A", {
  data("B_example")
  PMB <- getPM(B_example)
  
  expect_error(getMutexAB(A = NULL, B = B_example, PMB = PMB), 
               "not input matrix A")
})

test_that("getMutexAB handles NULL input matrix B", {
  data("A_example")
  PMA <- getPM(A_example)
  
  expect_error(getMutexAB(A = A_example, PMA = PMA, B = NULL), 
               "not input matrix B")
})

test_that("getMutexAB handles non-matrix input A", {
  data("B_example")
  PMB <- getPM(B_example)
  
  expect_error(getMutexAB(A = data.frame(a = 1, b = 2), B = B_example, PMB = PMB), 
               "input A must be a Matrix or a matrix class")
})

test_that("getMutexAB handles non-matrix input B", {
  data("A_example")
  PMA <- getPM(A_example)
  
  expect_error(getMutexAB(A = A_example, PMA = PMA, B = data.frame(a = 1, b = 2)), 
               "input B must be a Matrix or a matrix class")
})

test_that("getMutexAB handles empty matrix A", {
  data("B_example")
  PMB <- getPM(B_example)
  empty_matrix <- matrix(nrow = 0, ncol = 0)
  
  expect_error(getMutexAB(A = empty_matrix, B = B_example, PMB = PMB), 
               "input A must have at least 1 row and 1 column")
})

test_that("getMutexAB handles empty matrix B", {
  data("A_example")
  PMA <- getPM(A_example)
  empty_matrix <- matrix(nrow = 0, ncol = 0)
  
  expect_error(getMutexAB(A = A_example, PMA = PMA, B = empty_matrix), 
               "input B must have at least 1 row and 1 column")
})

test_that("getMutexAB handles non-binary matrix A", {
  data("B_example")
  PMB <- getPM(B_example)
  non_binary <- matrix(c(0, 1, 2, 3), nrow = 2, ncol = 2)
  
  expect_error(getMutexAB(A = non_binary, B = B_example, PMB = PMB), 
               "input A must be binary")
})

test_that("getMutexAB handles non-binary matrix B", {
  data("A_example")
  PMA <- getPM(A_example)
  non_binary <- matrix(c(0, 1, 2, 3), nrow = 2, ncol = 2)
  
  expect_error(getMutexAB(A = A_example, PMA = PMA, B = non_binary), 
               "input B must be binary")
})

test_that("getMutexAB handles invalid method", {
  data("A_example")
  data("B_example")
  
  PMA <- getPM(A_example)
  PMB <- getPM(B_example)
  
  expect_error(getMutexAB(A = A_example, PMA = PMA, B = B_example, PMB = PMB, 
                          method = "InvalidMethod"),
               'method must be "Exact", "RefinedNormal", "Binomial", "ShiftedBinomial"')
})

test_that("getMutexAB works with Matrix class input", {
  data("A_Matrix")
  data("B_Matrix")
  
  # Use smaller subset for faster testing
  A_small <- A_Matrix[1:50, 1:30]
  B_small <- B_Matrix[1:50, 1:30]
  
  PMA_small <- getPM(A_small)
  PMB_small <- getPM(B_small)
  
  result <- getMutexAB(A = A_small, PMA = PMA_small, B = B_small, PMB = PMB_small)
  
  expect_true(!is.null(result))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
  expect_equal(nrow(result), nrow(A_small))
  expect_equal(ncol(result), nrow(B_small))
})

test_that("getMutexAB handles mixed option", {
  data("A_example")
  data("B_example")
  
  PMA <- getPM(A_example)
  PMB <- getPM(B_example)
  
  result_mixed_true <- getMutexAB(A = A_example, PMA = PMA, B = B_example, 
                                  PMB = PMB, mixed = TRUE)
  result_mixed_false <- getMutexAB(A = A_example, PMA = PMA, B = B_example, 
                                   PMB = PMB, mixed = FALSE)
  
  expect_true(!is.null(result_mixed_true))
  expect_true(!is.null(result_mixed_false))
})

test_that("getMutexAB result dimensions match input dimensions", {
  data("A_example")
  data("B_example")
  
  # Use smaller subsets with different dimensions
  A_subset <- A_example[1:100, ]
  B_subset <- B_example[1:150, ]
  
  PMA <- getPM(A_subset)
  PMB <- getPM(B_subset)
  
  result <- getMutexAB(A = A_subset, PMA = PMA, B = B_subset, PMB = PMB)
  
  # Result should be 100 rows (genes in A) x 150 cols (genes in B)
  expect_equal(nrow(result), 100)
  expect_equal(ncol(result), 150)
})
