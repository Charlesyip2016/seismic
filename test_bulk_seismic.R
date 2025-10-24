#!/usr/bin/env Rscript
# Simple validation test for bulk-seismic functions
# This tests the core logic without requiring full package installation

cat("=== Bulk-seismic Validation Test ===\n\n")

# Source all required functions in order
cat("Loading functions...\n")
source("R/sweep_sparse.R")
source("R/calc_bulk_specificity.R")
source("R/run_bulk_seismic.R")

# Check that we have necessary libraries installed
required_pkgs <- c("data.table", "Matrix", "magrittr", "speedglm", "stats")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cat("ERROR: Missing required packages:", paste(missing_pkgs, collapse=", "), "\n")
  cat("Please install them using: install.packages(c('", paste(missing_pkgs, collapse="', '"), "'))\n", sep="")
  quit(status = 1)
}

# Load libraries
library(data.table)
library(Matrix)
library(magrittr)
library(speedglm)

cat("All required packages loaded successfully!\n\n")

# Test 1: Basic function validation
cat("Test 1: Validate function definitions...\n")
stopifnot(is.function(sweep_sparse))
stopifnot(is.function(calc_bulk_specificity))
stopifnot(is.function(run_bulk_seismic))
cat("✓ All functions are defined\n\n")

# Test 2: Check function arguments
cat("Test 2: Check function arguments...\n")
calc_bulk_args <- names(formals(calc_bulk_specificity))
expected_calc_bulk_args <- c("expression_matrix", "subtype_labels", "min_uniq_subtype", 
                              "min_subtype_size", "min_avg_exp_subtype")
stopifnot(all(expected_calc_bulk_args %in% calc_bulk_args))

run_bulk_args <- names(formals(run_bulk_seismic))
expected_run_bulk_args <- c("gene_trait_vector", "expression_matrix", "subtype_labels")
stopifnot(all(expected_run_bulk_args %in% run_bulk_args))
cat("✓ Function signatures are correct\n\n")

# Test 3: Create synthetic data
cat("Test 3: Create synthetic test data...\n")
set.seed(42)
n_genes <- 100
n_samples <- 30
gene_names <- paste0("GENE", 1:n_genes)
sample_names <- paste0("S", 1:n_samples)

# Expression matrix
expr_mat <- matrix(
  rlnorm(n_genes * n_samples, meanlog = 2, sdlog = 1),
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_names, sample_names)
)

# Subtype labels
subtype_labels <- c(rep("TypeA", 10), rep("TypeB", 10), rep("TypeC", 10))

# Add differential expression
expr_mat[1:20, 1:10] <- expr_mat[1:20, 1:10] * 2
expr_mat[21:40, 11:20] <- expr_mat[21:40, 11:20] * 2
expr_mat[41:60, 21:30] <- expr_mat[41:60, 21:30] * 2

cat("✓ Created expression matrix:", n_genes, "genes ×", n_samples, "samples\n")
cat("✓ Created", length(unique(subtype_labels)), "subtypes\n\n")

# Test 4: Run calc_bulk_specificity
cat("Test 4: Run calc_bulk_specificity...\n")
tryCatch({
  spec_scores <- calc_bulk_specificity(
    expression_matrix = expr_mat,
    subtype_labels = subtype_labels,
    min_uniq_subtype = 2,
    min_subtype_size = 3,
    min_avg_exp_subtype = 0.1
  )
  cat("✓ calc_bulk_specificity executed successfully\n")
  cat("  Output dimensions:", nrow(spec_scores), "genes ×", ncol(spec_scores), "subtypes\n")
  
  # Verify output structure
  stopifnot(is.matrix(spec_scores) || inherits(spec_scores, "Matrix"))
  stopifnot(ncol(spec_scores) == length(unique(subtype_labels)))
  
  # Check normalization: each gene's scores should sum to ~1 across subtypes
  row_sums <- rowSums(as.matrix(spec_scores))
  cat("  Row sums (should be ~1):", 
      "min =", round(min(row_sums), 3), 
      "max =", round(max(row_sums), 3), "\n")
  stopifnot(all(abs(row_sums - 1) < 0.01))  # Allow small numerical error
  
  cat("✓ Specificity scores are properly normalized\n\n")
}, error = function(e) {
  cat("✗ ERROR in calc_bulk_specificity:", conditionMessage(e), "\n")
  quit(status = 1)
})

# Test 5: Run full bulk-seismic analysis
cat("Test 5: Run full bulk-seismic analysis...\n")

# Create trait associations
trait_data <- data.frame(
  GENE = gene_names,
  ZSTAT = rnorm(n_genes, mean = 0, sd = 1)
)
# Make TypeA-specific genes have higher z-scores
trait_data$ZSTAT[1:20] <- trait_data$ZSTAT[1:20] + 1.5

tryCatch({
  results <- run_bulk_seismic(
    gene_trait_vector = trait_data,
    expression_matrix = expr_mat,
    subtype_labels = subtype_labels,
    gene_col = "GENE",
    score_col = "ZSTAT"
  )
  cat("✓ run_bulk_seismic executed successfully\n")
  
  # Verify output structure
  stopifnot(is.list(results))
  stopifnot(all(c("subtype_associations", "influential_genes", "specificity_scores") %in% names(results)))
  
  # Check subtype_associations
  stopifnot(is.data.frame(results$subtype_associations))
  stopifnot(all(c("subtype", "pvalue", "FDR") %in% names(results$subtype_associations)))
  stopifnot(nrow(results$subtype_associations) == length(unique(subtype_labels)))
  
  cat("  Found", nrow(results$subtype_associations), "subtype associations\n")
  cat("  P-values range:", 
      round(min(results$subtype_associations$pvalue), 4), "to",
      round(max(results$subtype_associations$pvalue), 4), "\n")
  
  # Check influential_genes
  stopifnot(is.list(results$influential_genes))
  if (length(results$influential_genes) > 0) {
    cat("  Found influential genes for", length(results$influential_genes), "significant subtype(s)\n")
  } else {
    cat("  No significant subtypes at p < 0.05\n")
  }
  
  cat("✓ Output structure is correct\n\n")
}, error = function(e) {
  cat("✗ ERROR in run_bulk_seismic:", conditionMessage(e), "\n")
  quit(status = 1)
})

# Test 6: Test with named vector input
cat("Test 6: Test with named vector input...\n")
trait_vector <- setNames(trait_data$ZSTAT, trait_data$GENE)

tryCatch({
  results2 <- run_bulk_seismic(
    gene_trait_vector = trait_vector,
    expression_matrix = expr_mat,
    subtype_labels = subtype_labels
  )
  cat("✓ Named vector input works correctly\n\n")
}, error = function(e) {
  cat("✗ ERROR with named vector:", conditionMessage(e), "\n")
  quit(status = 1)
})

# All tests passed
cat("=== All validation tests PASSED! ===\n")
cat("\nSummary:\n")
cat("✓ Function definitions are correct\n")
cat("✓ calc_bulk_specificity produces properly normalized specificity scores\n")
cat("✓ run_bulk_seismic performs complete analysis with correct output structure\n")
cat("✓ Both data.frame and named vector inputs work correctly\n")
cat("\nThe bulk-seismic implementation is ready to use!\n")
