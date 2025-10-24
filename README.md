# _seismic_: Single-cell Expression Investigation for Study of Molecular Interactions and Connections
This repository contains code for the R package _seismicGWAS_ a method for
calculating cell type-trait associations
given GWAS and single cell RNA-sequencing data.

## Citation

> Disentangling associations between complex traits and cell types with _seismic_.
> Lai Q, Dannenfelser R, Roussarie JP, Yao V. BioRxiv. April 2024.

## About

Integrating single-cell RNA sequencing (scRNA-seq) with Genome-Wide Association
Studies (GWAS) can help reveal GWAS-associated cell types furthering our
understanding of the cell-type-specific biological processes underlying complex
traits and disease. In order to rapidly and accurately pinpoint associations, we
develop a novel framework, _seismic_, which characterizes cell types using a new
specificity score. As part of the _seismic_ framework, the specific genes driving
cell type-trait associations can easily be accessed and analyzed, enabling further
biological insights. The following figure depicts a high level overview of
this process. 

![method overview](man/figures/seismic_overview.png)

## Installation
To install the seismic package first clone the seismic repo and then 
use devtools within R to point to seismic and install.

```R
devtools::install(path_to_seismic_folder)
library('seismicGWAS')
```

## Usage
Below we quickly show how to use `seismicGWAS` to calcuate cell
type-trait associations for the sample data included in the package. 
Full usage instructions, including a walk through of all major functions
can be found in the vignette.

```R
# calculate cell type specificity scores using included sample data
tmfacs_sscore <- calc_specificity(tmfacs_sce_small, ct_label_col='cluster_name')

# convert mouse gene identifiers to human ones that match data in GWAS summary data
# from MAGMA
tmfacs_sscore_hsa <- translate_gene_ids(tmfacs_sscore, from='mmu_symbol')

# calculate cell type-trait associations for type 2 diabetes
get_ct_trait_associations(tmfacs_sscore_hsa, t2d_magma)

# find the influential genes for a significant cell type-trait association
# in type 2 diabetes
find_inf_genes("Pancreas.beta cell", tmfacs_sscore_hsa, t2d_magma)
```

## Bulk-seismic: Adaptation for Bulk RNA-seq Data

In addition to the original seismic algorithm for single-cell RNA-seq, this package now includes **bulk-seismic**, an adaptation for bulk RNA-seq data. Instead of identifying cell type-trait associations, bulk-seismic identifies disease subtype-trait associations using bulk RNA-seq expression data from patient cohorts.

### Key Modifications

- **Input**: Bulk RNA-seq expression matrix (genes × samples) + disease subtype labels
- **Algorithm**: Removes the expression ratio component (r_i^(c)) and uses only expression consistency (p_i^(c))
- **Output**: Disease subtype-trait associations and influential genes

### Quick Start with Bulk-seismic

```R
# Prepare bulk RNA-seq data
expression_matrix <- ... # matrix with genes as rows, samples as columns
subtype_labels <- c("SubtypeA", "SubtypeA", "SubtypeB", ...) # one label per sample

# Run bulk-seismic analysis
results <- run_bulk_seismic(
  gene_trait_vector = magma_data,  # GWAS gene-level scores
  expression_matrix = expression_matrix,
  subtype_labels = subtype_labels
)

# View subtype-trait associations
print(results$subtype_associations)

# View influential genes for significant subtypes
print(results$influential_genes)
```

For detailed documentation on bulk-seismic, see [BULK_SEISMIC.md](BULK_SEISMIC.md) and [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md).
