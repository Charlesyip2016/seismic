#' Adjust specificity scores for batch effects
#'
#' Applies batch effect correction to specificity scores or expression data
#' before calculating associations.
#'
#' @param sscore A matrix of specificity scores (genes x cell types).
#' @param batch_info A data.frame with columns: cell_type, batch.
#' @param method Correction method: "mean_center", "quantile", or "scaling". 
#' Default: "mean_center".
#'
#' @return A matrix of batch-corrected specificity scores.
#' @export
correct_batch_effects <- function(sscore,
                                  batch_info,
                                  method = "mean_center") {
  
  sscore_mat <- as.matrix(sscore)
  batch_info <- as.data.table(batch_info)
  
  # Check that all cell types in sscore are in batch_info
  ct_in_sscore <- colnames(sscore_mat)
  ct_in_batch <- unique(batch_info$cell_type)
  
  if (!all(ct_in_sscore %in% ct_in_batch)) {
    warning("Some cell types in sscore not found in batch_info. ",
           "These will not be corrected.")
  }
  
  # Apply correction method
  corrected_mat <- sscore_mat
  
  if (method == "mean_center") {
    # For each batch, center the scores
    for (batch in unique(batch_info$batch)) {
      batch_cts <- batch_info[batch == !!batch, cell_type]
      batch_cts <- intersect(batch_cts, ct_in_sscore)
      
      if (length(batch_cts) > 0) {
        # Calculate batch mean
        batch_means <- rowMeans(sscore_mat[, batch_cts, drop = FALSE])
        overall_means <- rowMeans(sscore_mat)
        
        # Center
        for (ct in batch_cts) {
          corrected_mat[, ct] <- sscore_mat[, ct] - batch_means + overall_means
        }
      }
    }
    
  } else if (method == "quantile") {
    # Quantile normalization across batches
    for (gene in rownames(sscore_mat)) {
      gene_values <- sscore_mat[gene, ]
      
      # Rank within each batch
      for (batch in unique(batch_info$batch)) {
        batch_cts <- batch_info[batch == !!batch, cell_type]
        batch_cts <- intersect(batch_cts, ct_in_sscore)
        
        if (length(batch_cts) > 1) {
          batch_vals <- gene_values[batch_cts]
          ranks <- rank(batch_vals)
          # Map to overall quantiles
          overall_quantiles <- quantile(gene_values, 
                                       probs = (ranks - 0.5) / length(ranks))
          corrected_mat[gene, batch_cts] <- overall_quantiles
        }
      }
    }
    
  } else if (method == "scaling") {
    # Scale variance to match across batches
    for (batch in unique(batch_info$batch)) {
      batch_cts <- batch_info[batch == !!batch, cell_type]
      batch_cts <- intersect(batch_cts, ct_in_sscore)
      
      if (length(batch_cts) > 0) {
        batch_sds <- apply(sscore_mat[, batch_cts, drop = FALSE], 1, sd)
        overall_sds <- apply(sscore_mat, 1, sd)
        
        for (ct in batch_cts) {
          scaling_factor <- overall_sds / (batch_sds + 1e-10)
          corrected_mat[, ct] <- (sscore_mat[, ct] - rowMeans(sscore_mat[, batch_cts, drop = FALSE])) * 
            scaling_factor + rowMeans(sscore_mat)
        }
      }
    }
  } else {
    stop("Unknown method: ", method)
  }
  
  # Ensure non-negative values
  corrected_mat[corrected_mat < 0] <- 0
  
  # Renormalize columns
  col_sums <- colSums(corrected_mat)
  corrected_mat <- sweep(corrected_mat, 2, col_sums, "/")
  
  return(corrected_mat)
}


#' Perform stratified association analysis by batch
#'
#' Tests associations separately within each batch and then meta-analyzes.
#'
#' @param sscore_list A list of specificity score matrices, one per batch.
#' @param magma A data.frame or file path to MAGMA output.
#' @param batch_names Character vector of batch names. Default: names(sscore_list).
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#' @param meta_method Method for meta-analysis: "fixed" or "random". Default: "fixed".
#'
#' @return A data.frame with meta-analyzed association results.
#' @export
stratified_batch_analysis <- function(sscore_list,
                                     magma,
                                     batch_names = NULL,
                                     magma_gene_col = "GENE",
                                     magma_z_col = "ZSTAT",
                                     meta_method = "fixed") {
  
  pvalue <- FDR <- NULL # due to non-standard evaluation
  
  if (is.null(batch_names)) {
    batch_names <- names(sscore_list)
    if (is.null(batch_names)) {
      batch_names <- paste0("Batch", seq_along(sscore_list))
    }
  }
  
  # Load MAGMA
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # Run association for each batch
  batch_results <- list()
  
  for (i in seq_along(sscore_list)) {
    batch_name <- batch_names[i]
    sscore <- sscore_list[[i]]
    
    # Get associations
    assoc <- get_ct_trait_associations(sscore, magma, magma_gene_col, magma_z_col)
    assoc$batch <- batch_name
    
    batch_results[[i]] <- assoc
  }
  
  # Combine results
  all_results <- rbindlist(batch_results)
  
  # Get common cell types across batches
  ct_counts <- all_results[, .N, by = cell_type]
  common_cts <- ct_counts[N == length(sscore_list), cell_type]
  
  # Meta-analyze for common cell types
  meta_results <- rbindlist(lapply(common_cts, function(ct) {
    ct_results <- all_results[cell_type == ct]
    
    # Extract effect sizes (convert p-values to z-scores)
    z_scores <- stats::qnorm(1 - ct_results$pvalue)
    
    if (meta_method == "fixed") {
      # Fixed effects: inverse variance weighted
      # Assume equal variance for simplicity
      weights <- rep(1, length(z_scores))
      meta_z <- sum(z_scores * weights) / sqrt(sum(weights^2))
      meta_pval <- stats::pnorm(-abs(meta_z))
      
      # Heterogeneity test
      Q <- sum((z_scores - mean(z_scores))^2)
      Q_pval <- 1 - stats::pchisq(Q, df = length(z_scores) - 1)
      
    } else if (meta_method == "random") {
      # Random effects (DerSimonian-Laird)
      mean_z <- mean(z_scores)
      Q <- sum((z_scores - mean_z)^2)
      tau_sq <- max(0, (Q - (length(z_scores) - 1)) / length(z_scores))
      
      weights <- 1 / (1 + tau_sq)
      meta_z <- sum(z_scores * weights) / sqrt(sum(weights^2))
      meta_pval <- stats::pnorm(-abs(meta_z))
      Q_pval <- 1 - stats::pchisq(Q, df = length(z_scores) - 1)
    }
    
    return(data.table(
      cell_type = ct,
      meta_z = meta_z,
      pvalue = meta_pval,
      heterogeneity_Q = Q,
      heterogeneity_pval = Q_pval,
      n_batches = length(z_scores)
    ))
  }))
  
  meta_results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  meta_results <- meta_results[order(pvalue)]
  
  return(list(
    meta_results = meta_results,
    batch_results = all_results
  ))
}


#' Add covariates to association model
#'
#' Extends association testing to include covariates such as batch,
#' sex, age, or other confounders.
#'
#' @param sscore A matrix of specificity scores (genes x cell types).
#' @param magma A data.frame or file path to MAGMA output.
#' @param covariates A data.frame with gene-level covariates.
#' @param covariate_cols Character vector of covariate column names.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#'
#' @return A data.frame with covariate-adjusted associations.
#' @export
get_associations_with_covariates <- function(sscore,
                                            magma,
                                            covariates,
                                            covariate_cols,
                                            magma_gene_col = "GENE",
                                            magma_z_col = "ZSTAT") {
  
  pvalue <- FDR <- NULL # due to non-standard evaluation
  
  # Load MAGMA
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # Convert to data.table
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = TRUE)
  setnames(sscore, "rn", "gene")
  sscore <- melt(sscore, id.vars = "gene", variable.name = "cell_type",
                value.name = "specificity")
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = FALSE]
  names(magma) <- c("gene", "zstat")
  
  covariates <- as.data.table(covariates)
  
  # Merge data
  dt <- merge(sscore, magma, by = "gene")
  dt <- merge(dt, covariates, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  # Get cell types
  cts <- unique(dt$cell_type)
  
  # Test association with covariates for each cell type
  results <- rbindlist(lapply(cts, function(ct) {
    ct_data <- dt[cell_type == ct]
    
    if (nrow(ct_data) < 20) {
      return(data.table(
        cell_type = ct,
        beta_specificity = NA,
        pvalue = NA,
        n_genes = nrow(ct_data)
      ))
    }
    
    # Build formula with covariates
    formula_str <- paste("zstat ~ specificity +", 
                        paste(covariate_cols, collapse = " + "))
    formula_obj <- stats::as.formula(formula_str)
    
    # Fit model
    lm_out <- stats::lm(formula_obj, data = ct_data)
    lm_summ <- summary(lm_out)$coefficients
    
    # Extract specificity term
    beta <- lm_summ["specificity", 1]
    se <- lm_summ["specificity", 2]
    pval <- lm_summ["specificity", 4]
    
    # One-sided test
    pval_onesided <- if (beta > 0) pval / 2 else 1 - pval / 2
    
    # Also get covariate effects for reporting
    covar_effects <- lapply(covariate_cols, function(cov) {
      list(
        covariate = cov,
        beta = lm_summ[cov, 1],
        pval = lm_summ[cov, 4]
      )
    })
    
    return(data.table(
      cell_type = ct,
      beta_specificity = beta,
      se_specificity = se,
      pvalue = pval_onesided,
      n_genes = nrow(ct_data)
    ))
  }))
  
  results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  results <- results[order(pvalue)]
  
  return(results)
}


#' Calculate cross-dataset replicability
#'
#' Tests whether associations replicate across independent datasets.
#'
#' @param sscore1 Specificity scores from dataset 1.
#' @param sscore2 Specificity scores from dataset 2.
#' @param magma A data.frame or file path to MAGMA output.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#' @param fdr_threshold FDR threshold for significance. Default: 0.05.
#'
#' @return A list with replication statistics.
#' @export
assess_cross_dataset_replicability <- function(sscore1,
                                              sscore2,
                                              magma,
                                              magma_gene_col = "GENE",
                                              magma_z_col = "ZSTAT",
                                              fdr_threshold = 0.05) {
  
  # Get associations for both datasets
  assoc1 <- get_ct_trait_associations(sscore1, magma, magma_gene_col, magma_z_col)
  assoc2 <- get_ct_trait_associations(sscore2, magma, magma_gene_col, magma_z_col)
  
  # Get common cell types
  common_cts <- intersect(assoc1$cell_type, assoc2$cell_type)
  
  assoc1 <- assoc1[cell_type %in% common_cts]
  assoc2 <- assoc2[cell_type %in% common_cts]
  
  # Merge on cell type
  merged <- merge(assoc1, assoc2, by = "cell_type", suffixes = c("_ds1", "_ds2"))
  
  # Calculate correlation of p-values (on -log10 scale)
  cor_pval <- stats::cor(-log10(merged$pvalue_ds1), 
                        -log10(merged$pvalue_ds2),
                        method = "spearman")
  
  # Identify significant in dataset 1
  sig_ds1 <- merged[FDR_ds1 < fdr_threshold, cell_type]
  
  # Calculate replication rate
  if (length(sig_ds1) > 0) {
    replicated <- merged[cell_type %in% sig_ds1 & FDR_ds2 < fdr_threshold, cell_type]
    replication_rate <- length(replicated) / length(sig_ds1)
  } else {
    replication_rate <- NA
    replicated <- character()
  }
  
  # Concordance: both significant and same direction
  both_sig <- merged[FDR_ds1 < fdr_threshold & FDR_ds2 < fdr_threshold]
  
  result <- list(
    correlation_pvalues = cor_pval,
    n_common_celltypes = length(common_cts),
    n_significant_ds1 = length(sig_ds1),
    n_replicated = length(replicated),
    replication_rate = replication_rate,
    replicated_celltypes = replicated,
    both_significant = both_sig$cell_type,
    merged_results = merged
  )
  
  return(result)
}


#' Filter genes by expression consistency across batches
#'
#' Removes genes with highly variable expression patterns across batches
#' that might confound association analysis.
#'
#' @param sce_list A list of SingleCellExperiment objects, one per batch.
#' @param ct_label_col Column name for cell type labels.
#' @param assay_name Assay name. Default: 'logcounts'.
#' @param cv_threshold Coefficient of variation threshold. Default: 0.5.
#'
#' @return A character vector of genes passing the filter.
#' @export
filter_genes_by_batch_consistency <- function(sce_list,
                                             ct_label_col = "idents",
                                             assay_name = "logcounts",
                                             cv_threshold = 0.5) {
  
  # Get common genes across batches
  all_genes <- lapply(sce_list, function(sce) {
    rownames(SummarizedExperiment::assay(sce, assay_name))
  })
  common_genes <- Reduce(intersect, all_genes)
  
  # Get common cell types
  all_cts <- lapply(sce_list, function(sce) {
    unique(SummarizedExperiment::colData(sce)[[ct_label_col]])
  })
  common_cts <- Reduce(intersect, all_cts)
  
  # Calculate mean expression per cell type per batch
  batch_means <- list()
  
  for (i in seq_along(sce_list)) {
    sce <- sce_list[[i]]
    expr_mat <- SummarizedExperiment::assay(sce, assay_name)
    cell_meta <- SummarizedExperiment::colData(sce)
    
    ct_means <- matrix(0, nrow = length(common_genes), ncol = length(common_cts))
    rownames(ct_means) <- common_genes
    colnames(ct_means) <- common_cts
    
    for (ct in common_cts) {
      ct_cells <- which(cell_meta[[ct_label_col]] == ct)
      if (length(ct_cells) > 0) {
        ct_means[, ct] <- Matrix::rowMeans(expr_mat[common_genes, ct_cells, drop = FALSE])
      }
    }
    
    batch_means[[i]] <- ct_means
  }
  
  # Calculate CV across batches for each gene-celltype pair
  consistent_genes <- character()
  
  for (gene in common_genes) {
    gene_consistent <- TRUE
    
    for (ct in common_cts) {
      # Get values across batches
      batch_vals <- sapply(batch_means, function(bm) bm[gene, ct])
      
      # Calculate CV
      if (mean(batch_vals) > 0.1) {  # Only for reasonably expressed genes
        cv <- sd(batch_vals) / mean(batch_vals)
        
        if (cv > cv_threshold) {
          gene_consistent <- FALSE
          break
        }
      }
    }
    
    if (gene_consistent) {
      consistent_genes <- c(consistent_genes, gene)
    }
  }
  
  message("Retained ", length(consistent_genes), " out of ", 
         length(common_genes), " genes")
  
  return(consistent_genes)
}
