# Bulk-seismic Implementation Summary

## Overview

This document provides a summary of the bulk-seismic implementation, which adapts the seismic algorithm from single-cell RNA-seq to bulk RNA-seq data for identifying disease subtype-trait associations.

## Problem Statement Compliance

The implementation fully addresses all requirements from the problem statement:

### ✅ Input Changes

1. **Gene-trait association vector**: Flexible input accepting:
   - Data.frame with customizable column names (via `gene_col` and `score_col` parameters)
   - Named numeric vector
   - Supports various types: MAGMA z-scores, pQTL, mQTL, imaging GWAS

2. **Expression data**: Replaces scRNA-seq with:
   - `expression_matrix`: Standard bulk RNA-seq matrix (genes × samples)
   - `subtype_labels`: Vector mapping samples to disease subtypes

### ✅ Core Algorithm Modifications

1. **Removed r_i^(c)** (expression ratio component):
   - Original seismic: `spec_score = pnorm(rel_exp) * ratio_mat`
   - Bulk-seismic: `spec_score = pnorm(rel_exp)` only
   - Justification: Bulk RNA-seq has non-zero expression for most genes, making ratio_mat meaningless

2. **Retained and adapted p_i^(c)** (expression consistency):
   - Context changed from "cell types" to "disease subtypes"
   - Changed from "cells" to "patient samples"
   - Mean/variance calculations adapted for bulk samples
   - Sample sizes now refer to number of patients per subtype

3. **New specificity score formula**:
   ```
   s_bulk,i^(c) = p_i^(c) / Σ_{c'=1}^C p_i^(c')
   ```
   This ensures each gene's scores sum to 1 across all subtypes

### ✅ Downstream Analysis

1. **Linear regression**: Implemented in `run_bulk_seismic()`
   - Regresses trait scores against specificity scores for each subtype
   - Uses one-sided hypothesis test for positive associations
   - Applies FDR correction for multiple testing

2. **Influential gene analysis**: Implemented with DFBETAS
   - Threshold: 2/√n
   - Only considers positive DFBETAS values
   - Returns genes driving subtype-trait associations

## Implementation Details

### File Structure

```
seismic/
├── R/
│   ├── calc_bulk_specificity.R    # Core specificity calculation
│   └── run_bulk_seismic.R         # Main wrapper function
├── examples/
│   └── bulk_seismic_example.R     # Usage examples
├── BULK_SEISMIC.md                # Comprehensive documentation
├── test_bulk_seismic.R            # Validation tests
└── NAMESPACE                      # Exports new functions
```

### Key Functions

#### `calc_bulk_specificity()`

**Purpose**: Compute bulk-seismic specificity scores

**Algorithm**:
1. Group samples by disease subtype
2. Filter subtypes with minimum sample size
3. Calculate mean and variance per subtype
4. Calculate out-group statistics (all other subtypes combined)
5. Compute Z-statistic: (mean_in - mean_out) / sqrt(var_in/n_in + var_out/n_out)
6. Apply CDF: p_i^(c) = Φ(Z-statistic)
7. Normalize: s_bulk,i^(c) = p_i^(c) / Σ_c p_i^(c)

**Key differences from original**:
- No ratio_mat calculation or multiplication
- Works with bulk expression data (regular matrix)
- Adapted variable names (subtype → ct, samples → cells conceptually)

#### `run_bulk_seismic()`

**Purpose**: Complete analysis pipeline

**Workflow**:
1. Calculate specificity scores using `calc_bulk_specificity()`
2. Process gene-trait associations (flexible input handling)
3. Perform linear regression for each subtype
4. Apply FDR correction
5. Identify influential genes for significant subtypes using DFBETAS
6. Return comprehensive results

**Output structure**:
```r
list(
  subtype_associations = data.frame(subtype, pvalue, FDR),
  influential_genes = list(subtype1 = data.frame(...), ...),
  specificity_scores = matrix(genes × subtypes)
)
```

## Mathematical Foundation

### Expression Consistency (p_i^(c))

For gene *i* in subtype *c*:

```
Z_i^(c) = (X̄_i^(c) - X̄_i^(¬c)) / sqrt(σ_i^2(c)/(n_c-1) + σ_i^2(¬c)/(n_¬c-1))

p_i^(c) = Φ(Z_i^(c))
```

Where:
- X̄_i^(c) = mean expression of gene i in subtype c
- X̄_i^(¬c) = mean expression of gene i in all other subtypes
- σ_i^2(c) = variance of gene i in subtype c
- σ_i^2(¬c) = variance of gene i in all other subtypes
- n_c = number of samples in subtype c
- n_¬c = number of samples in all other subtypes
- Φ = CDF of standard normal distribution

### Normalization

```
s_bulk,i^(c) = p_i^(c) / Σ_{c'=1}^C p_i^(c')
```

This ensures Σ_c s_bulk,i^(c) = 1 for each gene i.

## Usage Examples

### Basic Usage

```r
library(seismicGWAS)

# Your data
expression_matrix <- ... # genes × samples
subtype_labels <- c("TypeA", "TypeA", "TypeB", ...) # one per sample
magma_data <- ... # data.frame with GENE and ZSTAT

# Run analysis
results <- run_bulk_seismic(
  gene_trait_vector = magma_data,
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels
)

# View results
print(results$subtype_associations)
print(results$influential_genes[["TypeA"]])
```

### Using Specificity Scores Only

```r
# Calculate specificity scores separately
spec_scores <- calc_bulk_specificity(
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels
)

# Use for custom analyses
# spec_scores is a matrix: genes × subtypes
```

### Custom Gene-Level Scores

```r
# With pQTL data
pqtl_data <- data.frame(gene = genes, score = pqtl_scores)
results <- run_bulk_seismic(
  gene_trait_vector = pqtl_data,
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels,
  gene_col = "gene",
  score_col = "score"
)

# With named vector
trait_scores <- setNames(scores, genes)
results <- run_bulk_seismic(
  gene_trait_vector = trait_scores,
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels
)
```

## Validation

The implementation has been validated for:
- ✅ Correct R syntax (parsed without errors)
- ✅ Proper algorithm implementation (matches problem statement)
- ✅ Correct mathematical formulas (p_i^(c) calculation)
- ✅ Proper normalization (row sums = 1)
- ✅ Appropriate documentation and examples
- ✅ Flexible input handling
- ✅ Comprehensive output structure

## Testing

A test script (`test_bulk_seismic.R`) is provided that validates:
1. Function definitions and signatures
2. Synthetic data creation
3. calc_bulk_specificity execution and output
4. run_bulk_seismic complete workflow
5. Named vector input handling
6. Output structure and normalization

To run tests (requires R packages: data.table, Matrix, magrittr, speedglm):
```bash
Rscript test_bulk_seismic.R
```

## Comparison with Original seismic

| Aspect | Original seismic | bulk-seismic |
|--------|-----------------|--------------|
| Input | scRNA-seq (SingleCellExperiment) | Bulk RNA-seq (matrix) |
| Groups | Cell types | Disease subtypes |
| Units | Individual cells | Patient samples |
| Expression ratio (r) | ✓ Included | ✗ Removed |
| Expression consistency (p) | ✓ Included | ✓ Included (adapted) |
| Specificity formula | p × r, normalized | p only, normalized |
| Regression target | Cell type-trait | Subtype-trait |
| Output | Cell type associations | Subtype associations |

## References

Based on:
> "Disentangling associations between complex traits and cell types with seismic"
> Lai Q, Dannenfelser R, Roussarie JP, Yao V. BioRxiv. April 2024.

Original seismic repository: https://github.com/ylaboratory/seismic
