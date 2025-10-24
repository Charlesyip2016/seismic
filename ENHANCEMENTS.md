# Seismic Enhancement Features

This document provides an overview of the major enhancements added to the seismic package for comprehensive single-cell GWAS analysis.

## Overview of New Features

The enhanced seismic package now includes 13 major innovation directions:

1. **Non-linear and Interaction Modeling**
2. **Multi-omics Integration**
3. **Hierarchical Modeling**
4. **Causal Inference (MR)**
5. **Enhanced Influence Analysis**
6. **Uncertainty Quantification**
7. **Rare Cell Type Support**
8. **Simulation Framework**

## 1. Non-linear and Interaction Modeling

### Purpose
Capture non-linear relationships between cell type specificity and trait risk that linear models might miss.

### Functions
- `get_ct_trait_associations_nonlinear()`: Fit GAM or kernel models
- `test_specificity_interactions()`: Test for interaction effects between specificity components

### Example
```r
# Use GAM to detect non-linear relationships
results_gam <- get_ct_trait_associations_nonlinear(
  sscore, magma_data, 
  model_type = "gam",
  gam_k = 10
)

# Test for interactions between proportion expressing and relative expression
interaction_results <- test_specificity_interactions(
  sce, magma_data,
  ct_label_col = "cell_type"
)
```

## 2. Multi-omics Integration

### Purpose
Integrate multiple omics layers (expression, ATAC-seq, eQTL, colocalization) for enhanced specificity scoring.

### Functions
- `calc_composite_specificity()`: Combine multiple omics scores
- `calc_gene_activity_from_atac()`: Convert ATAC peaks to gene scores
- `prepare_eqtl_scores()`: Prepare eQTL data for integration
- `prepare_coloc_scores()`: Prepare colocalization results

### Example
```r
# Calculate gene activity from ATAC-seq
atac_scores <- calc_gene_activity_from_atac(
  atac_sce, peak_to_gene,
  ct_label_col = "cell_type"
)

# Prepare eQTL scores
eqtl_scores <- prepare_eqtl_scores(
  eqtl_results,
  score_type = "pvalue"
)

# Combine into composite scores
composite_scores <- calc_composite_specificity(
  expr_specificity = sscore,
  atac_specificity = atac_scores,
  eqtl_scores = eqtl_scores,
  weights = c(expr = 1.0, atac = 0.5, eqtl = 0.5)
)

# Use composite scores for association
results <- get_ct_trait_associations(composite_scores, magma_data)
```

## 3. Hierarchical Modeling

### Purpose
Borrow information across related cell types in a taxonomy, improving estimates for rare cell types.

### Functions
- `get_hierarchical_associations()`: Hierarchical Bayesian modeling
- `calc_specificity_with_shrinkage()`: Empirical Bayes shrinkage for rare cells

### Example
```r
# Define cell type hierarchy
cell_hierarchy <- data.frame(
  cell_type = c("B_cell", "T_cell", "NK_cell", "Monocyte"),
  parent_type = c("Lymphoid", "Lymphoid", "Lymphoid", "Myeloid")
)

# Hierarchical association analysis
hier_results <- get_hierarchical_associations(
  sscore, magma_data, cell_hierarchy,
  tau_prior = 1.0
)

# Apply shrinkage for rare cell types
sscore_shrunk <- calc_specificity_with_shrinkage(
  sce,
  ct_label_col = "cell_type",
  cell_hierarchy = cell_hierarchy,
  min_cells_shrinkage = 50
)
```

## 4. Causal Inference

### Purpose
Use Mendelian Randomization to infer causal relationships between gene expression and traits.

### Functions
- `perform_mendelian_randomization()`: Two-sample MR analysis
- `test_pleiotropy_egger()`: Test for horizontal pleiotropy
- `perform_mediation_analysis()`: Mediation analysis
- `calc_causal_specificity_scores()`: Combine MR with specificity

### Example
```r
# Perform MR analysis
mr_results <- perform_mendelian_randomization(
  eqtl_data, gwas_data, 
  cell_type = "B_cell",
  f_threshold = 10
)

# Test for pleiotropy
pleiotropy <- test_pleiotropy_egger(
  eqtl_data, gwas_data,
  cell_type = "B_cell"
)

# Calculate causal-specificity scores
causal_scores <- calc_causal_specificity_scores(
  mr_results, sscore,
  mr_pval_threshold = 0.05
)
```

## 5. Enhanced Influence Analysis

### Purpose
Robust identification of influential genes using stability selection and leave-k-out validation.

### Functions
- `find_inf_genes_stable()`: DFBETAS with stability selection
- `leave_k_out_analysis()`: Leave-k-out validation
- `group_influence_analysis()`: Pathway-level influence

### Example
```r
# Find influential genes with stability assessment
inf_genes_stable <- find_inf_genes_stable(
  "B_cell", sscore, magma_data,
  n_bootstrap = 100,
  sampling_fraction = 0.8
)

# Leave-k-out analysis
lko_results <- leave_k_out_analysis(
  "B_cell", sscore, magma_data,
  k_percent = 5,
  n_iterations = 50
)

# Pathway-level influence
gene_groups <- list(
  immune_response = c("CD19", "CD79A", "MS4A1"),
  cell_cycle = c("MKI67", "PCNA", "TOP2A")
)

pathway_influence <- group_influence_analysis(
  "B_cell", sscore, magma_data, gene_groups
)
```

## 6. Uncertainty Quantification

### Purpose
Provide confidence intervals and assess calibration of p-values.

### Functions
- `get_associations_with_uncertainty()`: Bootstrap confidence intervals
- `assess_pvalue_calibration()`: QQ plot data and genomic inflation

### Example
```r
# Get associations with bootstrap CI
results_ci <- get_associations_with_uncertainty(
  sscore, magma_data,
  n_bootstrap = 1000,
  ci_level = 0.95
)

# Assess p-value calibration
calibration <- assess_pvalue_calibration(results)

# Check genomic inflation factor
print(calibration$lambda)
```

## 7. Simulation Framework

### Purpose
Benchmark methods on data with known ground truth.

### Functions
- `simulate_sc_data()`: Generate synthetic scRNA-seq data
- `simulate_gwas_data()`: Generate synthetic GWAS summary stats
- `benchmark_seismic()`: Test performance on simulated data
- `generate_simulation_suite()`: Create multiple scenarios

### Example
```r
# Generate simulated data
sim_sc <- simulate_sc_data(
  n_genes = 5000,
  n_cells_per_type = c(100, 200, 150, 80),
  n_driver_genes = 50,
  driver_effect_size = 2.0,
  seed = 42
)

# Calculate specificity
sscore_sim <- calc_specificity(sim_sc$sce, ct_label_col = "cell_type")

# Generate GWAS data
gwas_sim <- simulate_gwas_data(
  rownames(sscore_sim),
  sim_sc$driver_genes,
  sscore_sim,
  effect_size = 0.5,
  seed = 42
)

# Benchmark performance
benchmark_results <- benchmark_seismic(
  sim_sc, gwas_sim,
  methods = c("linear", "gam", "hierarchical"),
  fdr_threshold = 0.05
)

# View performance metrics
print(benchmark_results$performance)
```

## Integration Example: Complete Workflow

Here's how to use multiple enhancements together:

```r
# 1. Calculate enhanced specificity with shrinkage for rare cells
sscore <- calc_specificity_with_shrinkage(
  sce,
  ct_label_col = "cell_type",
  cell_hierarchy = hierarchy,
  min_cells_shrinkage = 50
)

# 2. Integrate multi-omics data
composite_scores <- calc_composite_specificity(
  expr_specificity = sscore,
  atac_specificity = atac_scores,
  eqtl_scores = eqtl_scores,
  coloc_scores = coloc_scores
)

# 3. Run hierarchical association analysis with uncertainty
results <- get_hierarchical_associations(
  composite_scores, magma_data, hierarchy
)

# 4. Test non-linearity for significant associations
sig_cts <- results[FDR < 0.05, cell_type]
nonlinear_tests <- get_ct_trait_associations_nonlinear(
  composite_scores[, sig_cts], magma_data,
  model_type = "gam"
)

# 5. Find influential genes with stability
for (ct in sig_cts) {
  inf_genes <- find_inf_genes_stable(
    ct, composite_scores, magma_data,
    n_bootstrap = 100
  )
  
  # Save results
  write.csv(inf_genes, paste0("influential_genes_", ct, ".csv"))
}

# 6. Perform causal inference for top associations
mr_results <- perform_mendelian_randomization(
  eqtl_data, gwas_data,
  cell_type = sig_cts[1],
  f_threshold = 10
)

# 7. Assess calibration
calibration <- assess_pvalue_calibration(results)
plot(calibration$expected, calibration$observed, 
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
abline(0, 1, col = "red")
```

## Performance Considerations

1. **Non-linear models**: GAM fitting can be slow for large datasets. Consider using on filtered or top cell types.

2. **Bootstrap methods**: Set `n_bootstrap` appropriately based on available compute time. 100-1000 is usually sufficient.

3. **Multi-omics**: Missing data in some omics layers is handled automatically through masking.

4. **Hierarchical models**: Variational inference is used for computational efficiency.

## Interpreting Results

### Model Comparison
- Compare AIC/BIC between linear and non-linear models
- If ΔAIC > 2, non-linearity is likely present

### Stability Scores
- High (>0.75): Core driver genes
- Moderate (0.5-0.75): Potential drivers, needs validation
- Low (<0.5): Unstable/noise

### MR Results
- Check F-statistic > 10 for instrument strength
- Egger intercept p < 0.05 indicates pleiotropy
- Use multiple cell types to assess specificity

### Calibration
- Lambda (genomic inflation) should be ~1.0
- Lambda > 1.1 may indicate population stratification or model misspecification

## References

See the main seismic paper and vignette for basic usage. These enhancements implement methods from:

- GAM: Wood (2006) Generalized Additive Models
- MR: Davey Smith & Hemani (2014) Mendelian randomization
- Hierarchical modeling: Gelman & Hill (2006) Data Analysis Using Regression
- Stability selection: Meinshausen & Bühlmann (2010) Stability selection
