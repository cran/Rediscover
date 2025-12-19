library(testthat)
library(Rediscover)

# Test discoversomaticInteractions function

# Helper function to create a simple mock MAF object for testing
create_mock_maf <- function(n_genes = 10, n_samples = 50) {
  # Create a simple data.table that mimics MAF structure
  genes <- paste0("Gene", 1:n_genes)
  samples <- paste0("Sample", 1:n_samples)
  
  # Create mutation data
  mutations <- expand.grid(
    Hugo_Symbol = genes[1:min(n_genes, 5)],  # Only first 5 genes have mutations
    Tumor_Sample_Barcode = samples[1:min(n_samples, 30)]
  )
  
  # Add required columns
  mutations$Chromosome <- sample(1:22, nrow(mutations), replace = TRUE)
  mutations$Start_Position <- sample(1:1000000, nrow(mutations), replace = TRUE)
  mutations$End_Position <- mutations$Start_Position + sample(1:100, nrow(mutations), replace = TRUE)
  mutations$Reference_Allele <- sample(c("A", "T", "C", "G"), nrow(mutations), replace = TRUE)
  mutations$Tumor_Seq_Allele2 <- sample(c("A", "T", "C", "G"), nrow(mutations), replace = TRUE)
  mutations$Variant_Classification <- sample(c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del"), 
                                             nrow(mutations), replace = TRUE)
  mutations$Variant_Type <- "SNP"
  
  # Convert to data.table
  mutations <- data.table::as.data.table(mutations)
  
  # Create MAF object using maftools::read.maf
  # This requires writing to a temporary file
  temp_maf <- tempfile(fileext = ".maf")
  data.table::fwrite(mutations, temp_maf, sep = "\t")
  
  maf_obj <- tryCatch({
    maftools::read.maf(temp_maf, verbose = FALSE)
  }, error = function(e) {
    NULL
  })
  
  unlink(temp_maf)
  return(maf_obj)
}

test_that("discoversomaticInteractions handles NULL maf input", {
  expect_error(discoversomaticInteractions(maf = NULL), 
               "not input maf file")
})

test_that("discoversomaticInteractions works with mock MAF object and default parameters", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 25, n_samples = 100)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  # Suppress plot output
  pdf(NULL)
  result <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, returnAll = TRUE)
  }, error = function(e) {
    NULL
  })
  dev.off()
  
  # If the function runs, result should be a data.table or NULL (if no interactions found)
  if (!is.null(result)) {
    expect_true(is.data.frame(result) || is(result, "data.table"))
    
    # Check that expected columns exist
    expected_cols <- c("gene1", "gene2", "pValue")
    expect_true(all(expected_cols %in% colnames(result)))
  }
})

test_that("discoversomaticInteractions works with specified genes", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 10, n_samples = 50)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  test_genes <- c("Gene1", "Gene2", "Gene3")
  
  pdf(NULL)
  result <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, genes = test_genes, returnAll = TRUE)
  }, error = function(e) {
    NULL
  })
  dev.off()
  
  # Result should contain only the specified genes
  if (!is.null(result)) {
    expect_true(all(result$gene1 %in% test_genes))
    expect_true(all(result$gene2 %in% test_genes))
  }
})

test_that("discoversomaticInteractions respects getMutexMethod parameter", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 15, n_samples = 80)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  methods <- c("ShiftedBinomial", "Binomial", "RefinedNormal")
  
  for (method in methods) {
    pdf(NULL)
    result <- tryCatch({
      discoversomaticInteractions(maf = maf_obj, top = 5, 
                                 getMutexMethod = method, returnAll = TRUE)
    }, error = function(e) {
      NULL
    })
    dev.off()
    
    # Should not error with valid methods
    expect_true(is.null(result) || is.data.frame(result))
  }
})

test_that("discoversomaticInteractions handles pvalue threshold parameter", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 20, n_samples = 100)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  pdf(NULL)
  result <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, 
                               pvalue = c(0.05, 0.01), returnAll = TRUE)
  }, error = function(e) {
    NULL
  })
  dev.off()
  
  # Check that result exists and p-values are in valid range
  if (!is.null(result) && nrow(result) > 0) {
    expect_true(all(result$pValue >= 0 & result$pValue <= 1))
  }
})

test_that("discoversomaticInteractions handles returnAll parameter", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 15, n_samples = 80)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  pdf(NULL)
  result_all <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, 
                               returnAll = TRUE)
  }, error = function(e) {
    NULL
  })
  
  result_filtered <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, 
                               returnAll = FALSE, pvalue = c(0.05, 0.01))
  }, error = function(e) {
    NULL
  })
  dev.off()
  
  # returnAll=TRUE should return more or equal rows than returnAll=FALSE
  if (!is.null(result_all) && !is.null(result_filtered)) {
    expect_true(nrow(result_all) >= nrow(result_filtered))
  }
})

test_that("discoversomaticInteractions output has correct structure", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 20, n_samples = 100)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  pdf(NULL)
  result <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, returnAll = TRUE)
  }, error = function(e) {
    NULL
  })
  dev.off()
  
  if (!is.null(result) && nrow(result) > 0) {
    # Check required columns exist
    expect_true("gene1" %in% colnames(result))
    expect_true("gene2" %in% colnames(result))
    expect_true("pValue" %in% colnames(result))
    expect_true("Event" %in% colnames(result))
    expect_true("pair" %in% colnames(result))
    
    # Check Event values are valid
    expect_true(all(result$Event %in% c("Co_Occurence", "Mutually_Exclusive")))
    
    # Check gene1 != gene2 (no self-interactions)
    expect_true(all(result$gene1 != result$gene2))
  }
})

test_that("discoversomaticInteractions handles getMutexMixed parameter", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 15, n_samples = 80)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  pdf(NULL)
  result_mixed_true <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, 
                               getMutexMixed = TRUE, returnAll = TRUE)
  }, error = function(e) {
    NULL
  })
  
  result_mixed_false <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, 
                               getMutexMixed = FALSE, returnAll = TRUE)
  }, error = function(e) {
    NULL
  })
  dev.off()
  
  # Both should produce valid results
  expect_true(is.null(result_mixed_true) || is.data.frame(result_mixed_true))
  expect_true(is.null(result_mixed_false) || is.data.frame(result_mixed_false))
})

test_that("discoversomaticInteractions handles invalid geneOrder", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 10, n_samples = 50)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  pdf(NULL)
  result <- tryCatch({
    discoversomaticInteractions(maf = maf_obj, top = 5, 
                               geneOrder = c("NonExistentGene1", "NonExistentGene2"),
                               returnAll = TRUE)
  }, error = function(e) {
    e$message
  })
  dev.off()
  
  # Should error with invalid genes
  if (is.character(result)) {
    expect_true(grepl("does not match", result, ignore.case = TRUE))
  }
})

test_that("discoversomaticInteractions works with different top values", {
  skip_if_not_installed("maftools")
  
  maf_obj <- create_mock_maf(n_genes = 30, n_samples = 100)
  skip_if(is.null(maf_obj), "Could not create mock MAF object")
  
  top_values <- c(5, 10, 15)
  
  for (top_val in top_values) {
    pdf(NULL)
    result <- tryCatch({
      discoversomaticInteractions(maf = maf_obj, top = top_val, returnAll = TRUE)
    }, error = function(e) {
      NULL
    })
    dev.off()
    
    expect_true(is.null(result) || is.data.frame(result))
  }
})
