# Bulk-seismic: Adapting seismic for Bulk RNA-seq Data

## Overview

The `bulk-seismic` functions adapt the original seismic algorithm to work with bulk RNA-seq data instead of single-cell RNA-seq data. The key algorithmic modifications are:

1. **Input Change**: Instead of a single-cell RNA-seq matrix with cell-level data, bulk-seismic takes:
   - A bulk RNA-seq expression matrix (genes × patient samples)
   - A vector of disease subtype labels (one per sample)

2. **Algorithm Modification**: 
   - **Removed**: The expression ratio component r_i^(c) (proportion of cells expressing a gene), which is not meaningful for bulk RNA-seq data
   - **Retained**: The expression consistency component p_i^(c), adapted for patient samples within disease subtypes rather than cells within cell types
   - **New specificity score**: s_bulk,i^(c) = p_i^(c) / Σ_c p_i^(c)

3. **Downstream Analysis**: Linear regression and influential gene analysis remain the same, but now identify disease subtype-trait associations instead of cell type-trait associations.

## Core Functions

### `calc_bulk_specificity()`

Computes bulk-seismic specificity scores for genes across disease subtypes.

**Parameters:**
- `expression_matrix`: A matrix with genes as rows and patient samples as columns
- `subtype_labels`: A character vector indicating disease subtype for each sample
- `min_uniq_subtype`: Minimum number of unique subtypes (default: 2)
- `min_subtype_size`: Minimum samples per subtype (default: 3)
- `min_avg_exp_subtype`: Minimum mean expression threshold (default: 0.1)

**Returns:** A matrix of specificity scores (genes × subtypes)

### `run_bulk_seismic()`

Main wrapper function that performs the complete bulk-seismic analysis.

**Parameters:**
- `gene_trait_vector`: Gene-level trait association scores (e.g., MAGMA z-scores, pQTL, mQTL, imaging GWAS)
- `expression_matrix`: Bulk RNA-seq expression matrix
- `subtype_labels`: Disease subtype labels for samples
- `gene_col`: Column name for gene identifiers (default: 'GENE')
- `score_col`: Column name for association scores (default: 'ZSTAT')
- `significance_threshold`: P-value threshold for significance (default: 0.05)

**Returns:** A list containing:
- `subtype_associations`: Statistical associations between subtypes and trait
- `influential_genes`: List of influential genes for each significant subtype
- `specificity_scores`: The computed specificity score matrix

## Example Usage

```r
library(seismicGWAS)

# Prepare your data
# expression_matrix: genes (rows) x samples (columns)
# subtype_labels: vector of disease subtypes
# magma_data: data.frame with GENE and ZSTAT columns

# Run bulk-seismic analysis
results <- run_bulk_seismic(
  gene_trait_vector = magma_data,
  expression_matrix = bulk_rnaseq_matrix,
  subtype_labels = patient_subtypes
)

# View subtype-trait associations
print(results$subtype_associations)

# View influential genes for the most significant subtype
top_subtype <- results$subtype_associations$subtype[1]
print(results$influential_genes[[top_subtype]])

# Access the specificity scores
spec_scores <- results$specificity_scores
```

## Using Custom Gene-Level Scores

The `gene_trait_vector` parameter is flexible and can accept various types of gene-level association scores:

```r
# Example with pQTL data
pqtl_data <- data.frame(
  gene_id = gene_names,
  pqtl_score = pqtl_scores
)

results <- run_bulk_seismic(
  gene_trait_vector = pqtl_data,
  expression_matrix = bulk_rnaseq_matrix,
  subtype_labels = patient_subtypes,
  gene_col = 'gene_id',
  score_col = 'pqtl_score'
)

# Example with named vector
trait_scores <- setNames(z_scores, gene_names)
results <- run_bulk_seismic(
  gene_trait_vector = trait_scores,
  expression_matrix = bulk_rnaseq_matrix,
  subtype_labels = patient_subtypes
)
```

## Algorithm Details

### Expression Consistency (p_i^(c))

For each gene i in each disease subtype c, we calculate the probability that the gene's expression is consistently higher in subtype c compared to all other subtypes. This is computed using a Z-test based on:

- Mean expression: X̄_i^(c) - mean expression of gene i in subtype c
- Variance: σ_i^2(c) - variance of gene i's expression in subtype c  
- Sample size: |L^(c)| - number of patient samples in subtype c

The Z-statistic compares the mean expression in subtype c versus the mean expression in all other subtypes combined, accounting for their respective variances and sample sizes. The cumulative distribution function (CDF) of the standard normal distribution gives us p_i^(c).

### Specificity Score Normalization

The final bulk-seismic specificity score for gene i in subtype c is:

s_bulk,i^(c) = p_i^(c) / Σ_{c'=1}^C p_i^(c')

This normalization ensures that for each gene, the specificity scores across all subtypes sum to 1, making them interpretable as relative specificities.

### Linear Regression and Influential Genes

Following the calculation of specificity scores, we:

1. Perform linear regression: Z^(c) ~ S_bulk^(c) for each subtype
2. Extract p-values using a one-sided test (testing for positive associations)
3. Apply FDR correction for multiple testing
4. For significant subtypes, identify influential genes using DFBETAS > 2/√n

## Differences from Original seismic

| Aspect | Original seismic | bulk-seismic |
|--------|-----------------|--------------|
| Input data | Single-cell RNA-seq | Bulk RNA-seq |
| Groups | Cell types | Disease subtypes (patient groups) |
| Expression ratio (r) | Included (proportion of cells expressing gene) | Removed (not meaningful for bulk) |
| Specificity score | p × r, normalized | p only, normalized |
| Interpretation | Cell type-trait associations | Disease subtype-trait associations |
