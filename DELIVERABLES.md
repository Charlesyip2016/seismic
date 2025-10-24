# Bulk-seismic Implementation - Deliverables

This document summarizes the deliverables for the bulk-seismic algorithm implementation, as requested in the problem statement.

## 1. Algorithm Workflow Description

The bulk-seismic algorithm adapts the original seismic framework from single-cell RNA-seq to bulk RNA-seq data. Here's the workflow:

### Workflow Steps:

1. **Data Input**:
   - Bulk RNA-seq expression matrix (genes × patient samples)
   - Disease subtype labels (one per sample)
   - Gene-trait association vector (e.g., MAGMA z-scores, pQTL, mQTL, or imaging GWAS scores)

2. **Specificity Score Calculation** (via `calc_bulk_specificity`):
   - Group patient samples by disease subtype
   - Filter subtypes with sufficient sample size
   - Calculate mean and variance of gene expression per subtype
   - Calculate out-group statistics (all other subtypes combined)
   - Compute Z-statistic comparing in-group vs out-group expression
   - Apply CDF to get p_i^(c) = Φ(Z-statistic)
   - Normalize across subtypes: s_bulk,i^(c) = p_i^(c) / Σ_c p_i^(c)

3. **Regression Analysis** (via `run_bulk_seismic`):
   - For each subtype, perform linear regression: trait_score ~ specificity_score
   - Extract p-values using one-sided hypothesis test (for positive associations)
   - Apply FDR correction for multiple testing

4. **Influential Gene Identification**:
   - For significant subtypes (p < threshold), calculate DFBETAS for each gene
   - Identify influential genes with DFBETAS > 2/√n
   - Return ordered list of genes with their contribution to the association

### Key Algorithmic Differences from Original seismic:

| Component | Original seismic | bulk-seismic |
|-----------|-----------------|--------------|
| Input data | scRNA-seq (SingleCellExperiment) | Bulk RNA-seq (matrix) |
| Grouping variable | Cell types | Disease subtypes |
| Expression ratio (r_i^(c)) | ✓ Used (% cells expressing gene) | ✗ Removed (not meaningful for bulk) |
| Expression consistency (p_i^(c)) | ✓ Used | ✓ Used (adapted for bulk samples) |
| Final score | p_i^(c) × r_i^(c), normalized | p_i^(c) only, normalized |

## 2. Core Function: `calc_bulk_specificity()`

**File**: `R/calc_bulk_specificity.R`

### Function Signature:
```r
calc_bulk_specificity(
  expression_matrix,      # Matrix: genes × samples
  subtype_labels,         # Vector: subtype per sample
  min_uniq_subtype = 2,   # Minimum number of subtypes
  min_subtype_size = 3,   # Minimum samples per subtype
  min_avg_exp_subtype = 0.1  # Minimum mean expression
)
```

### Algorithm Implementation:

The function implements the modified seismic specificity calculation:

1. **Data Validation and Grouping**:
   ```r
   subtype_groups <- data.table(sample = colnames(expression_matrix),
                                 subtype = subtype_labels)
   ```

2. **Mean Expression Calculation** (per subtype):
   ```r
   factor_mat <- Matrix::fac2sparse(factor(subtype_groups$subtype, ...))
   sum_mat <- Matrix::t(expression_matrix_sparse %*% Matrix::t(factor_mat))
   mean_mat <- sweep_sparse(sum_mat, margin = 1, stats = subtype_groups_n$N, fun = "/")
   ```

3. **Variance Calculation** (per subtype):
   ```r
   var_mat <- (Matrix::t(expression_matrix_sparse^2 %*% Matrix::t(factor_mat)) - 
               2 * mean_mat * sum_mat +
               sweep_sparse(mean_mat^2, margin = 1, stats = subtype_groups_n$N, fun = "*")) %>%
             sweep_sparse(margin = 1, stats = subtype_groups_n$N - 1, fun = "/")
   ```

4. **Out-Group Statistics**:
   ```r
   out_group_mat <- matrix(1, nrow = n_subtypes, ncol = n_subtypes) - diag(n_subtypes)
   out_mean <- (out_group_mat %*% sum_mat) %>%
               sweep_sparse(margin = 1, stats = out_group_mat %*% subtype_groups_n$N, fun = "/")
   ```

5. **Expression Consistency (p_i^(c))**:
   ```r
   rel_exp <- (mean_mat - out_mean) /
              sqrt(sweep_sparse(var_mat, margin = 1, stats = subtype_groups_n$N - 1, fun = "/") +
                   sweep_sparse(out_variance, margin = 1, stats = out_group_mat %*% 
                                subtype_groups_n$N - 1, fun = "/"))
   p_mat <- stats::pnorm(as.matrix(rel_exp))
   ```

6. **Normalization** (KEY MODIFICATION - no ratio_mat multiplication):
   ```r
   # Original seismic: spec_score <- pnorm(rel_exp) * ratio_mat
   # Bulk-seismic: spec_score <- pnorm(rel_exp) only
   spec_score <- sweep_sparse(x = p_mat, margin = 2, 
                             stats = Matrix::colSums(p_mat), fun = "/")
   ```

### Returns:
Matrix of specificity scores (genes × subtypes), where each row sums to 1.

## 3. High-Level Function: `run_bulk_seismic()`

**File**: `R/run_bulk_seismic.R`

### Function Signature:
```r
run_bulk_seismic(
  gene_trait_vector,       # Data.frame or named vector
  expression_matrix,       # Matrix: genes × samples
  subtype_labels,          # Vector: subtype per sample
  gene_col = "GENE",       # Column name for genes
  score_col = "ZSTAT",     # Column name for scores
  significance_threshold = 0.05,  # P-value cutoff
  min_uniq_subtype = 2,
  min_subtype_size = 3,
  min_avg_exp_subtype = 0.1
)
```

### Implementation Steps:

1. **Calculate Specificity Scores**:
   ```r
   sscore <- calc_bulk_specificity(expression_matrix, subtype_labels, ...)
   ```

2. **Process Gene-Trait Vector** (flexible input):
   ```r
   # Handles both data.frame and named vector
   if (is.vector(gene_trait_vector) && !is.null(names(gene_trait_vector))) {
     trait_dt <- data.table(gene = names(gene_trait_vector), 
                           score = gene_trait_vector)
   } else if (is.data.frame(gene_trait_vector)) {
     trait_dt <- gene_trait_vector[, c(gene_col, score_col)]
   }
   ```

3. **Linear Regression Analysis**:
   ```r
   for each subtype:
     slm <- speedglm::speedlm(trait_score ~ specificity_score)
     pval <- one_sided_test(slm)  # Test for positive association
   ```

4. **FDR Correction**:
   ```r
   subtype_results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
   ```

5. **Influential Gene Analysis** (for significant subtypes):
   ```r
   lm_out <- stats::lm(score ~ specificity, data = subtype_data)
   dfbetas <- stats::dfbetas(lm_out)[, 2]
   is_influential <- dfbetas > 2 / sqrt(n_genes)
   ```

### Returns:
```r
list(
  subtype_associations = data.frame(
    subtype,   # Subtype name
    pvalue,    # Association p-value
    FDR        # FDR-corrected p-value
  ),
  influential_genes = list(
    subtype1 = data.frame(
      gene,          # Gene name
      specificity,   # Specificity score
      score,         # Trait score
      dfbetas,       # Influence metric
      is_influential # Boolean flag
    ),
    ...
  ),
  specificity_scores = matrix(genes × subtypes)
)
```

## 4. Usage Examples

### Example 1: Basic MAGMA Analysis
```r
library(seismicGWAS)

# Load your data
bulk_expr <- read.csv("bulk_expression.csv", row.names = 1)  # genes × samples
subtypes <- read.csv("patient_subtypes.csv")$subtype
magma_data <- read.table("gwas.genes.out", header = TRUE)

# Run analysis
results <- run_bulk_seismic(
  gene_trait_vector = magma_data,
  expression_matrix = as.matrix(bulk_expr),
  subtype_labels = subtypes
)

# View subtype associations
print(results$subtype_associations)

# View influential genes for most significant subtype
top_subtype <- results$subtype_associations$subtype[1]
influential <- results$influential_genes[[top_subtype]]
print(influential[influential$is_influential, ])
```

### Example 2: Custom Gene-Level Scores (pQTL)
```r
# pQTL data
pqtl_scores <- data.frame(
  gene_symbol = gene_names,
  pqtl_zscore = z_scores
)

# Run analysis
results <- run_bulk_seismic(
  gene_trait_vector = pqtl_scores,
  expression_matrix = bulk_expr,
  subtype_labels = subtypes,
  gene_col = "gene_symbol",
  score_col = "pqtl_zscore"
)
```

### Example 3: Named Vector Input
```r
# Named vector of trait scores
trait_scores <- setNames(magma_data$ZSTAT, magma_data$GENE)

results <- run_bulk_seismic(
  gene_trait_vector = trait_scores,
  expression_matrix = bulk_expr,
  subtype_labels = subtypes
)
```

### Example 4: Specificity Scores Only
```r
# Calculate specificity scores separately for custom analysis
spec_scores <- calc_bulk_specificity(
  expression_matrix = bulk_expr,
  subtype_labels = subtypes,
  min_subtype_size = 5
)

# spec_scores is a matrix: genes (rows) × subtypes (columns)
# Each row sums to 1
```

## 5. Documentation Files

The implementation includes comprehensive documentation:

1. **BULK_SEISMIC.md**: User guide with:
   - Overview and algorithm modifications
   - Function descriptions and parameters
   - Usage examples for various scenarios
   - Mathematical foundation
   - Comparison with original seismic

2. **IMPLEMENTATION_SUMMARY.md**: Technical summary with:
   - Problem statement compliance checklist
   - Implementation details
   - Mathematical formulas
   - Validation information
   - File structure

3. **README.md**: Updated with bulk-seismic section

4. **examples/bulk_seismic_example.R**: Runnable example script with synthetic data

5. **test_bulk_seismic.R**: Validation test suite

## 6. Key Algorithm Modifications Summary

### Modification 1: Removed r_i^(c) (Expression Ratio)

**Original seismic** (line 126-129 in calc_specificity.R):
```r
ratio_mat <- Matrix::t((data_mat > 0) %*% Matrix::t(factor_mat)) %>%
  sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/")
```

**Bulk-seismic**: Completely removed this component.

**Justification**: In bulk RNA-seq, nearly all genes have non-zero expression across samples, making the ratio meaningless. This metric only makes sense for single-cell data where dropout is common.

### Modification 2: Adapted p_i^(c) Context

**Original seismic**:
- c = cell type
- Units = individual cells
- N = number of cells per type

**Bulk-seismic**:
- c = disease subtype
- Units = patient samples
- N = number of samples per subtype

The calculation remains the same, but the interpretation changes.

### Modification 3: New Specificity Score Formula

**Original seismic** (line 156-157 in calc_specificity.R):
```r
spec_score <- stats::pnorm(as.matrix(rel_exp)) * ratio_mat
spec_score <- sweep_sparse(x = spec_score, margin = 2, 
                          stats = Matrix::colSums(spec_score), fun = "/")
```

**Bulk-seismic** (line 154-157 in calc_bulk_specificity.R):
```r
p_mat <- stats::pnorm(as.matrix(rel_exp))
spec_score <- sweep_sparse(x = p_mat, margin = 2, 
                          stats = Matrix::colSums(p_mat), fun = "/")
```

**Formula**: s_bulk,i^(c) = p_i^(c) / Σ_{c'=1}^C p_i^(c')

This ensures each gene's specificity scores sum to 1 across all subtypes.

## 7. Validation and Testing

The implementation has been validated for:

✅ Correct R syntax (all files parse without errors)
✅ Proper algorithm implementation (matches problem statement exactly)
✅ Mathematical correctness (p_i^(c) calculation verified)
✅ Normalization correctness (row sums = 1)
✅ Flexible input handling (data.frame and named vector)
✅ Comprehensive documentation
✅ Example code and test suite

A validation test script is provided that checks:
- Function definitions and signatures
- Execution with synthetic data
- Output structure and properties
- Normalization correctness
- Multiple input format handling

## 8. Files Added/Modified

**New Files**:
- `R/calc_bulk_specificity.R` - Core specificity calculation (161 lines)
- `R/run_bulk_seismic.R` - Main wrapper function (197 lines)
- `BULK_SEISMIC.md` - User documentation (143 lines)
- `IMPLEMENTATION_SUMMARY.md` - Technical summary (243 lines)
- `examples/bulk_seismic_example.R` - Usage examples (105 lines)
- `test_bulk_seismic.R` - Validation tests (179 lines)

**Modified Files**:
- `NAMESPACE` - Added exports for new functions
- `README.md` - Added bulk-seismic section

**Total**: 1,063 lines of code and documentation added

## Conclusion

This implementation fully addresses all requirements from the problem statement:

1. ✅ Flexible gene-trait association input (MAGMA, pQTL, mQTL, imaging GWAS)
2. ✅ Bulk RNA-seq expression matrix + subtype labels as input
3. ✅ Removed r_i^(c) component completely
4. ✅ Retained and adapted p_i^(c) for disease subtypes
5. ✅ New specificity score: s_bulk,i^(c) = p_i^(c) / Σ_c p_i^(c)
6. ✅ Preserved downstream regression and influential gene analysis
7. ✅ Comprehensive documentation and examples
8. ✅ Working implementation with validation tests

The bulk-seismic algorithm is ready for use in identifying disease subtype-GWAS associations using bulk RNA-seq data.
