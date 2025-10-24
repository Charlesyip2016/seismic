#!/usr/bin/env Rscript
# Example usage of bulk-seismic functions
# This script demonstrates how to use the new bulk-seismic functions

library(seismicGWAS)

# Generate synthetic example data for demonstration
set.seed(123)

# 1. Create synthetic bulk RNA-seq expression matrix
n_genes <- 500
n_samples <- 60
gene_names <- paste0("Gene_", 1:n_genes)
sample_names <- paste0("Sample_", 1:n_samples)

# Create expression matrix (genes x samples)
expression_matrix <- matrix(
  rlnorm(n_genes * n_samples, meanlog = 2, sdlog = 1),
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_names, sample_names)
)

# 2. Create disease subtype labels (3 subtypes)
subtype_labels <- c(
  rep("SubtypeA", 20),
  rep("SubtypeB", 20),
  rep("SubtypeC", 20)
)

# Add some differential expression between subtypes
# Make some genes more expressed in SubtypeA
expression_matrix[1:50, 1:20] <- expression_matrix[1:50, 1:20] * 2

# Make some genes more expressed in SubtypeB
expression_matrix[51:100, 21:40] <- expression_matrix[51:100, 21:40] * 2

# Make some genes more expressed in SubtypeC
expression_matrix[101:150, 41:60] <- expression_matrix[101:150, 41:60] * 2

# 3. Create synthetic gene-trait association data (like MAGMA z-scores)
# Simulate stronger associations for subtype-specific genes
trait_associations <- data.frame(
  GENE = gene_names,
  ZSTAT = rnorm(n_genes, mean = 0, sd = 1)
)

# Make subtype-specific genes have higher z-scores
trait_associations$ZSTAT[1:50] <- trait_associations$ZSTAT[1:50] + 1.5
trait_associations$ZSTAT[51:100] <- trait_associations$ZSTAT[51:100] + 1.2
trait_associations$ZSTAT[101:150] <- trait_associations$ZSTAT[101:150] + 1.0

# 4. Run bulk-seismic analysis
cat("Running bulk-seismic analysis...\n")
results <- run_bulk_seismic(
  gene_trait_vector = trait_associations,
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels,
  gene_col = "GENE",
  score_col = "ZSTAT",
  significance_threshold = 0.05
)

# 5. Display results
cat("\n=== Subtype-Trait Associations ===\n")
print(results$subtype_associations)

cat("\n=== Specificity Scores (first 10 genes, all subtypes) ===\n")
print(head(results$specificity_scores, 10))

# 6. Display influential genes for significant subtypes
if (length(results$influential_genes) > 0) {
  cat("\n=== Influential Genes for Significant Subtypes ===\n")
  for (subtype in names(results$influential_genes)) {
    cat("\nSubtype:", subtype, "\n")
    influential <- results$influential_genes[[subtype]]
    cat("Number of influential genes:", sum(influential$is_influential), "\n")
    cat("Top 10 genes by DFBETAS:\n")
    print(head(influential, 10))
  }
} else {
  cat("\nNo significant subtypes found at p < 0.05\n")
}

# 7. Alternative usage: Using only calc_bulk_specificity
cat("\n=== Alternative: Using calc_bulk_specificity directly ===\n")
specificity_scores <- calc_bulk_specificity(
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels
)
cat("Computed specificity scores for", nrow(specificity_scores), "genes and", 
    ncol(specificity_scores), "subtypes\n")
print(head(specificity_scores))

# 8. Example with named vector for gene-trait associations
cat("\n=== Alternative: Using named vector for gene-trait associations ===\n")
trait_vector <- setNames(trait_associations$ZSTAT, trait_associations$GENE)
results2 <- run_bulk_seismic(
  gene_trait_vector = trait_vector,
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels
)
cat("Analysis completed successfully with named vector input\n")

cat("\n=== Example completed successfully! ===\n")
