#' Calculate composite multi-omics specificity scores
#'
#' Integrates multiple omics layers (expression, accessibility, eQTL, colocalization)
#' into a unified specificity score for enhanced cell type-trait association.
#'
#' @param expr_specificity A matrix of expression specificity scores (genes x cell types).
#' @param atac_specificity Optional matrix of ATAC-seq accessibility scores (genes x cell types).
#' @param eqtl_scores Optional matrix of eQTL evidence scores (genes x cell types).
#' @param coloc_scores Optional matrix of colocalization PP4 values (genes x cell types).
#' @param weights Named vector of weights for each omics layer. 
#' Default: c(expr=1.0, atac=0.5, eqtl=0.5, coloc=0.8).
#' @param gamma_powers Named vector of exponents for each component.
#' Default: all 1.0 (linear combination in log space).
#' @param normalize Logical, whether to normalize final scores. Default: TRUE.
#'
#' @return A matrix of composite specificity scores (genes x cell types).
#' @export
calc_composite_specificity <- function(expr_specificity,
                                      atac_specificity = NULL,
                                      eqtl_scores = NULL,
                                      coloc_scores = NULL,
                                      weights = c(expr = 1.0, atac = 0.5, 
                                                 eqtl = 0.5, coloc = 0.8),
                                      gamma_powers = c(expr = 1.0, atac = 1.0,
                                                      eqtl = 1.0, coloc = 1.0),
                                      normalize = TRUE) {
  
  # Convert to matrix if needed
  if (inherits(expr_specificity, "dgCMatrix") || inherits(expr_specificity, "dgeMatrix")) {
    expr_specificity <- as.matrix(expr_specificity)
  }
  
  genes <- rownames(expr_specificity)
  cell_types <- colnames(expr_specificity)
  
  # Initialize composite score with expression component
  # Formula: F = w_E * E^γ1 * (1 + w_A * A)^γ2 * (1 + w_Q * Q)^γ3 * (1 + w_C * C)^γ4
  # In log space: log F = log w_E + γ1*log E + γ2*log(1 + w_A*A) + ...
  
  # Start with expression (primary component)
  log_composite <- log(weights["expr"] + 1e-10) + 
                   gamma_powers["expr"] * log(expr_specificity + 1e-10)
  
  # Add ATAC component if provided
  if (!is.null(atac_specificity)) {
    if (inherits(atac_specificity, "dgCMatrix") || inherits(atac_specificity, "dgeMatrix")) {
      atac_specificity <- as.matrix(atac_specificity)
    }
    
    # Align genes and cell types
    common_genes <- intersect(genes, rownames(atac_specificity))
    common_cts <- intersect(cell_types, colnames(atac_specificity))
    
    if (length(common_genes) > 0 && length(common_cts) > 0) {
      atac_aligned <- matrix(0, nrow = length(genes), ncol = length(cell_types))
      rownames(atac_aligned) <- genes
      colnames(atac_aligned) <- cell_types
      atac_aligned[common_genes, common_cts] <- atac_specificity[common_genes, common_cts]
      
      log_composite <- log_composite + 
        gamma_powers["atac"] * log(1 + weights["atac"] * atac_aligned)
    }
  }
  
  # Add eQTL component if provided
  if (!is.null(eqtl_scores)) {
    if (inherits(eqtl_scores, "dgCMatrix") || inherits(eqtl_scores, "dgeMatrix")) {
      eqtl_scores <- as.matrix(eqtl_scores)
    }
    
    common_genes <- intersect(genes, rownames(eqtl_scores))
    common_cts <- intersect(cell_types, colnames(eqtl_scores))
    
    if (length(common_genes) > 0 && length(common_cts) > 0) {
      eqtl_aligned <- matrix(0, nrow = length(genes), ncol = length(cell_types))
      rownames(eqtl_aligned) <- genes
      colnames(eqtl_aligned) <- cell_types
      eqtl_aligned[common_genes, common_cts] <- eqtl_scores[common_genes, common_cts]
      
      log_composite <- log_composite + 
        gamma_powers["eqtl"] * log(1 + weights["eqtl"] * eqtl_aligned)
    }
  }
  
  # Add colocalization component if provided
  if (!is.null(coloc_scores)) {
    if (inherits(coloc_scores, "dgCMatrix") || inherits(coloc_scores, "dgeMatrix")) {
      coloc_scores <- as.matrix(coloc_scores)
    }
    
    common_genes <- intersect(genes, rownames(coloc_scores))
    common_cts <- intersect(cell_types, colnames(coloc_scores))
    
    if (length(common_genes) > 0 && length(common_cts) > 0) {
      coloc_aligned <- matrix(0, nrow = length(genes), ncol = length(cell_types))
      rownames(coloc_aligned) <- genes
      colnames(coloc_aligned) <- cell_types
      coloc_aligned[common_genes, common_cts] <- coloc_scores[common_genes, common_cts]
      
      log_composite <- log_composite + 
        gamma_powers["coloc"] * log(1 + weights["coloc"] * coloc_aligned)
    }
  }
  
  # Convert back from log space
  composite_scores <- exp(log_composite)
  
  # Normalize if requested
  if (normalize) {
    # Normalize each column (cell type) to sum to 1
    col_sums <- colSums(composite_scores)
    composite_scores <- sweep(composite_scores, 2, col_sums, "/")
  }
  
  return(composite_scores)
}


#' Calculate gene activity scores from ATAC-seq data
#'
#' Converts peak-level ATAC-seq accessibility into gene-level activity scores
#' based on peak-to-gene linkage.
#'
#' @param atac_sce A SingleCellExperiment object with ATAC-seq peak counts.
#' @param peak_to_gene A data.frame with columns: peak_id, gene_id, distance, correlation.
#' @param ct_label_col Column name in colData for cell type labels.
#' @param assay_name Assay name in atac_sce. Default: "counts".
#' @param distance_weight Logical, whether to weight by distance. Default: TRUE.
#' @param max_distance Maximum peak-gene distance to consider. Default: 100000.
#'
#' @return A matrix of gene activity scores (genes x cell types).
#' @export
calc_gene_activity_from_atac <- function(atac_sce, 
                                        peak_to_gene,
                                        ct_label_col = "idents",
                                        assay_name = "counts",
                                        distance_weight = TRUE,
                                        max_distance = 100000) {
  
  # Extract ATAC data
  atac_mat <- SummarizedExperiment::assay(atac_sce, assay_name)
  cell_meta <- SummarizedExperiment::colData(atac_sce)
  
  # Get cell type groupings
  ct_groups <- data.table(
    cell = rownames(cell_meta),
    ct = cell_meta[[ct_label_col]], key = "ct"
  )
  
  ct_groups_n <- ct_groups[, .N, by = ct]
  
  # Filter peak-to-gene linkages by distance
  peak_to_gene <- as.data.table(peak_to_gene)
  if ("distance" %in% names(peak_to_gene)) {
    peak_to_gene <- peak_to_gene[abs(distance) <= max_distance]
  }
  
  # Calculate mean accessibility per cell type for each peak
  factor_mat <- Matrix::fac2sparse(factor(ct_groups$ct, levels = unique(ct_groups$ct)))
  sum_mat <- Matrix::t(atac_mat %*% Matrix::t(factor_mat))
  mean_mat <- sum_mat %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/")
  
  # Map peaks to genes
  unique_genes <- unique(peak_to_gene$gene_id)
  unique_cts <- unique(ct_groups$ct)
  
  gene_activity <- matrix(0, nrow = length(unique_genes), ncol = length(unique_cts))
  rownames(gene_activity) <- unique_genes
  colnames(gene_activity) <- unique_cts
  
  # Aggregate peak accessibility to gene level
  for (gene in unique_genes) {
    linked_peaks <- peak_to_gene[gene_id == gene]
    
    # Find peaks that exist in our data
    peak_idx <- match(linked_peaks$peak_id, rownames(atac_mat))
    peak_idx <- peak_idx[!is.na(peak_idx)]
    
    if (length(peak_idx) > 0) {
      # Weight by distance if requested
      if (distance_weight && "distance" %in% names(linked_peaks)) {
        weights <- 1 / (1 + abs(linked_peaks$distance[!is.na(match(linked_peaks$peak_id, 
                                                                    rownames(atac_mat)))]) / 10000)
      } else {
        weights <- rep(1, length(peak_idx))
      }
      
      # Weight by correlation if available
      if ("correlation" %in% names(linked_peaks)) {
        corr_vals <- linked_peaks$correlation[!is.na(match(linked_peaks$peak_id, 
                                                           rownames(atac_mat)))]
        corr_vals[is.na(corr_vals)] <- 0
        corr_vals[corr_vals < 0] <- 0  # Only positive correlations
        weights <- weights * (corr_vals + 0.1)  # Add small constant to avoid zeros
      }
      
      # Weighted sum of peak accessibilities
      peak_values <- as.matrix(mean_mat[, peak_idx, drop = FALSE])
      gene_activity[gene, ] <- as.vector(peak_values %*% weights / sum(weights))
    }
  }
  
  return(gene_activity)
}


#' Prepare eQTL scores for integration
#'
#' Converts eQTL association results into cell-type-specific gene scores
#' suitable for composite specificity calculation.
#'
#' @param eqtl_results A data.frame with columns: gene, cell_type, pvalue (or beta, se).
#' @param score_type Type of score to compute: "pvalue", "pip" (posterior inclusion prob),
#' or "effect". Default: "pvalue".
#' @param genes Vector of gene IDs to include. If NULL, uses all genes in eqtl_results.
#' @param cell_types Vector of cell type names. If NULL, uses all in eqtl_results.
#'
#' @return A matrix of eQTL scores (genes x cell types).
#' @export
prepare_eqtl_scores <- function(eqtl_results, 
                               score_type = "pvalue",
                               genes = NULL,
                               cell_types = NULL) {
  
  eqtl_results <- as.data.table(eqtl_results)
  
  # Get unique genes and cell types
  if (is.null(genes)) {
    genes <- unique(eqtl_results$gene)
  }
  if (is.null(cell_types)) {
    cell_types <- unique(eqtl_results$cell_type)
  }
  
  # Initialize score matrix
  score_matrix <- matrix(0, nrow = length(genes), ncol = length(cell_types))
  rownames(score_matrix) <- genes
  colnames(score_matrix) <- cell_types
  
  # Calculate scores based on type
  if (score_type == "pvalue") {
    # Convert p-values to -log10(p) scores
    if (!"pvalue" %in% names(eqtl_results)) {
      stop("eqtl_results must have 'pvalue' column for score_type='pvalue'")
    }
    
    for (i in 1:nrow(eqtl_results)) {
      gene <- eqtl_results$gene[i]
      ct <- eqtl_results$cell_type[i]
      if (gene %in% genes && ct %in% cell_types) {
        pval <- max(eqtl_results$pvalue[i], 1e-300)  # Avoid log(0)
        score_matrix[gene, ct] <- -log10(pval)
      }
    }
    
  } else if (score_type == "pip") {
    # Use posterior inclusion probability directly
    if (!"pip" %in% names(eqtl_results)) {
      stop("eqtl_results must have 'pip' column for score_type='pip'")
    }
    
    for (i in 1:nrow(eqtl_results)) {
      gene <- eqtl_results$gene[i]
      ct <- eqtl_results$cell_type[i]
      if (gene %in% genes && ct %in% cell_types) {
        score_matrix[gene, ct] <- eqtl_results$pip[i]
      }
    }
    
  } else if (score_type == "effect") {
    # Use absolute effect size weighted by significance
    if (!all(c("beta", "se") %in% names(eqtl_results))) {
      stop("eqtl_results must have 'beta' and 'se' columns for score_type='effect'")
    }
    
    for (i in 1:nrow(eqtl_results)) {
      gene <- eqtl_results$gene[i]
      ct <- eqtl_results$cell_type[i]
      if (gene %in% genes && ct %in% cell_types) {
        z_score <- abs(eqtl_results$beta[i] / eqtl_results$se[i])
        score_matrix[gene, ct] <- z_score
      }
    }
  }
  
  # Normalize scores to [0, 1] range per cell type
  for (j in 1:ncol(score_matrix)) {
    max_val <- max(score_matrix[, j])
    if (max_val > 0) {
      score_matrix[, j] <- score_matrix[, j] / max_val
    }
  }
  
  return(score_matrix)
}


#' Prepare colocalization scores for integration
#'
#' Converts colocalization analysis results (e.g., from coloc package) into
#' gene-cell type scores for composite specificity.
#'
#' @param coloc_results A data.frame with columns: gene, cell_type, PP4 (or PP3, PP4).
#' @param pp_threshold Minimum PP4 to consider significant. Default: 0.5.
#' @param genes Vector of gene IDs. If NULL, uses all in coloc_results.
#' @param cell_types Vector of cell type names. If NULL, uses all in coloc_results.
#'
#' @return A matrix of colocalization scores (genes x cell types).
#' @export
prepare_coloc_scores <- function(coloc_results,
                                pp_threshold = 0.5,
                                genes = NULL,
                                cell_types = NULL) {
  
  coloc_results <- as.data.table(coloc_results)
  
  # Get unique genes and cell types
  if (is.null(genes)) {
    genes <- unique(coloc_results$gene)
  }
  if (is.null(cell_types)) {
    cell_types <- unique(coloc_results$cell_type)
  }
  
  # Initialize score matrix
  score_matrix <- matrix(0, nrow = length(genes), ncol = length(cell_types))
  rownames(score_matrix) <- genes
  colnames(score_matrix) <- cell_types
  
  # Fill in PP4 values
  if (!"PP4" %in% names(coloc_results)) {
    stop("coloc_results must have 'PP4' column")
  }
  
  for (i in 1:nrow(coloc_results)) {
    gene <- coloc_results$gene[i]
    ct <- coloc_results$cell_type[i]
    if (gene %in% genes && ct %in% cell_types) {
      pp4 <- coloc_results$PP4[i]
      # Only include if above threshold
      if (pp4 >= pp_threshold) {
        score_matrix[gene, ct] <- pp4
      }
    }
  }
  
  return(score_matrix)
}
