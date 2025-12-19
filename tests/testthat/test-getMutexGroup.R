library(testthat)
library(Rediscover)

# Test getMutexGroup function

test_that("getMutexGroup works with Impurity type", {
  # Load and prepare example data as in the function documentation
  data("A_example")
  A2 <- A_example[, 1:30]
  A2[1, 1:10] <- 1
  A2[2, 1:10] <- 0
  A2[3, 1:10] <- 0
  A2[1, 11:20] <- 0
  A2[2, 11:20] <- 1
  A2[3, 11:20] <- 0
  A2[1, 21:30] <- 0
  A2[2, 21:30] <- 0
  A2[3, 21:30] <- 1
  
  PM2 <- getPM(A2)
  A <- A2[1:3, ]
  PM <- PM2[1:3, ]
  
  result <- getMutexGroup(A, PM, type = "Impurity")
  
  # Check that result is not null
  expect_true(!is.null(result))
  
  # Check that result is a single p-value (numeric of length 1)
  expect_true(is.numeric(result))
  expect_length(result, 1)
  
  # Check that p-value is in valid range [0, 1]
  expect_true(result >= 0 & result <= 1)
})

test_that("getMutexGroup works with Coverage type", {
  data("A_example")
  A2 <- A_example[, 1:30]
  A2[1, 1:10] <- 1
  A2[2, 1:10] <- 0
  A2[3, 1:10] <- 0
  A2[1, 11:20] <- 0
  A2[2, 11:20] <- 1
  A2[3, 11:20] <- 0
  A2[1, 21:30] <- 0
  A2[2, 21:30] <- 0
  A2[3, 21:30] <- 1
  
  PM2 <- getPM(A2)
  A <- A2[1:3, ]
  PM <- PM2[1:3, ]
  
  result <- getMutexGroup(A, PM, type = "Coverage")
  
  expect_true(!is.null(result))
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result >= 0 & result <= 1)
})

test_that("getMutexGroup works with Exclusivity type", {
  data("A_example")
  A2 <- A_example[, 1:30]
  A2[1, 1:10] <- 1
  A2[2, 1:10] <- 0
  A2[3, 1:10] <- 0
  A2[1, 11:20] <- 0
  A2[2, 11:20] <- 1
  A2[3, 11:20] <- 0
  A2[1, 21:30] <- 0
  A2[2, 21:30] <- 0
  A2[3, 21:30] <- 1
  
  PM2 <- getPM(A2)
  A <- A2[1:3, ]
  PM <- PM2[1:3, ]
  
  result <- getMutexGroup(A, PM, type = "Exclusivity")
  
  expect_true(!is.null(result))
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result >= 0 & result <= 1)
})

test_that("getMutexGroup works with default type (Impurity)", {
  data("A_example")
  A2 <- A_example[, 1:30]
  PM2 <- getPM(A2)
  A <- A2[1:5, ]
  PM <- PM2[1:5, ]
  
  result <- getMutexGroup(A, PM)
  
  expect_true(!is.null(result))
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result >= 0 & result <= 1)
})

test_that("getMutexGroup works with lower.tail = FALSE", {
  data("A_example")
  A2 <- A_example[, 1:30]
  PM2 <- getPM(A2)
  A <- A2[1:3, ]
  PM <- PM2[1:3, ]
  
  result <- getMutexGroup(A, PM, type = "Impurity", lower.tail = FALSE)
  
  expect_true(!is.null(result))
  expect_true(is.numeric(result))
  expect_true(result >= 0 & result <= 1)
})

test_that("getMutexGroup handles NULL input matrix A", {
  expect_error(getMutexGroup(A = NULL), "not input matrix A")
})

test_that("getMutexGroup handles NULL input matrix PM", {
  data("A_example")
  expect_error(getMutexGroup(A = A_example, PM = NULL), 
               "not probability matrix PM")
})

test_that("getMutexGroup handles non-matrix input A", {
  data("A_example")
  PM <- getPM(A_example)
  
  expect_error(getMutexGroup(A = data.frame(a = 1, b = 2), PM = PM), 
               "input A must be a Matrix or a matrix class")
})

test_that("getMutexGroup handles non-matrix input PM", {
  data("A_example")
  
  expect_error(getMutexGroup(A = A_example, PM = data.frame(a = 1, b = 2)), 
               "input PM must be Matrix, matrix or PMatrix class")
})

test_that("getMutexGroup handles empty matrix A", {
  empty_matrix <- matrix(nrow = 0, ncol = 0)
  PM <- matrix(0.5, nrow = 1, ncol = 1)
  
  expect_error(getMutexGroup(A = empty_matrix, PM = PM), 
               "input A must have at least 1 row and 1 column")
})

test_that("getMutexGroup handles empty matrix PM", {
  data("A_example")
  empty_matrix <- matrix(nrow = 0, ncol = 0)
  
  expect_error(getMutexGroup(A = A_example, PM = empty_matrix), 
               "input PM must have at least 1 row and 1 column")
})

test_that("getMutexGroup handles non-binary matrix A", {
  non_binary <- matrix(c(0, 1, 2, 3), nrow = 2, ncol = 2)
  PM <- matrix(0.5, nrow = 2, ncol = 2)
  
  expect_error(getMutexGroup(A = non_binary, PM = PM), 
               "input A must be binary")
})

test_that("getMutexGroup handles non-binary matrix PM", {
  A <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)
  non_binary_PM <- matrix(c(0, 1, 2, 3), nrow = 2, ncol = 2)
  
  expect_error(getMutexGroup(A = A, PM = non_binary_PM), 
               "input PM must be binary")
})

test_that("getMutexGroup works with Matrix class input", {
  data("A_Matrix")
  
  # Use smaller subset for faster testing
  A_small <- A_Matrix[1:10, 1:50]
  PM_small <- getPM(A_small)
  
  result <- getMutexGroup(A = A_small, PM = PM_small, type = "Impurity")
  
  expect_true(!is.null(result))
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result >= 0 & result <= 1)
})

test_that("getMutexGroup returns consistent results for all types", {
  data("A_example")
  A2 <- A_example[, 1:30]
  PM2 <- getPM(A2)
  A <- A2[1:3, ]
  PM <- PM2[1:3, ]
  
  result_impurity <- getMutexGroup(A, PM, type = "Impurity")
  result_coverage <- getMutexGroup(A, PM, type = "Coverage")
  result_exclusivity <- getMutexGroup(A, PM, type = "Exclusivity")
  
  # All should return valid p-values
  expect_true(all(c(result_impurity, result_coverage, result_exclusivity) >= 0))
  expect_true(all(c(result_impurity, result_coverage, result_exclusivity) <= 1))
})

test_that("getMutexGroup works with small group of genes", {
  data("A_example")
  
  # Test with just 2 genes
  A_small <- A_example[1:2, 1:20]
  PM_small <- getPM(A_small)
  
  result <- getMutexGroup(A_small, PM_small, type = "Impurity")
  
  expect_true(!is.null(result))
  expect_true(is.numeric(result))
  expect_true(result >= 0 & result <= 1)
})

test_that("getMutexGroup works with larger group of genes", {
  data("A_example")
  
  # Test with larger group
  A_large <- A_example[1:20, 1:50]
  PM_large <- getPM(A_large)
  
  result <- getMutexGroup(A_large, PM_large, type = "Coverage")
  
  expect_true(!is.null(result))
  expect_true(is.numeric(result))
  expect_true(result >= 0 & result <= 1)
})
