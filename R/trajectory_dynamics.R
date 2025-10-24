#' Calculate trajectory-specific gene expression patterns
#'
#' Identifies genes with stage-specific expression patterns along a trajectory
#' (pseudotime) and tests association with trait risk.
#'
#' @param sce A SingleCellExperiment object with trajectory information.
#' @param pseudotime_col Column name in colData containing pseudotime values.
#' @param ct_label_col Column name for cell type labels.
#' @param assay_name Assay name. Default: 'logcounts'.
#' @param n_windows Number of pseudotime windows to test. Default: 5.
#' @param window_overlap Fractional overlap between windows. Default: 0.2.
#'
#' @return A matrix of stage-specific expression scores (genes x cell_types x stages).
#' @export
calc_trajectory_specificity <- function(sce,
                                       pseudotime_col,
                                       ct_label_col = "idents",
                                       assay_name = "logcounts",
                                       n_windows = 5,
                                       window_overlap = 0.2) {
  
  # Extract data
  expr_mat <- SummarizedExperiment::assay(sce, assay_name)
  cell_meta <- SummarizedExperiment::colData(sce)
  
  if (!pseudotime_col %in% names(cell_meta)) {
    stop("Pseudotime column '", pseudotime_col, "' not found in colData")
  }
  
  pseudotime <- cell_meta[[pseudotime_col]]
  cell_types <- cell_meta[[ct_label_col]]
  
  # Remove cells with missing pseudotime
  valid_cells <- !is.na(pseudotime)
  expr_mat <- expr_mat[, valid_cells]
  pseudotime <- pseudotime[valid_cells]
  cell_types <- cell_types[valid_cells]
  
  # Define pseudotime windows
  pt_range <- range(pseudotime)
  window_size <- (pt_range[2] - pt_range[1]) / (n_windows - window_overlap * (n_windows - 1))
  
  windows <- list()
  for (i in 1:n_windows) {
    start <- pt_range[1] + (i - 1) * window_size * (1 - window_overlap)
    end <- start + window_size
    windows[[i]] <- list(
      id = paste0("stage", i),
      start = start,
      end = end,
      cells = which(pseudotime >= start & pseudotime <= end)
    )
  }
  
  # Calculate expression for each cell type and stage
  unique_cts <- unique(cell_types)
  gene_names <- rownames(expr_mat)
  
  # Store results as list of matrices (one per stage)
  stage_specificity <- list()
  
  for (w in 1:n_windows) {
    window <- windows[[w]]
    stage_name <- window$id
    
    # Initialize matrix for this stage
    stage_mat <- matrix(0, nrow = length(gene_names), ncol = length(unique_cts))
    rownames(stage_mat) <- gene_names
    colnames(stage_mat) <- unique_cts
    
    for (ct in unique_cts) {
      # Get cells of this type in this window
      ct_cells <- which(cell_types == ct)
      window_ct_cells <- intersect(window$cells, ct_cells)
      
      if (length(window_ct_cells) < 5) {
        # Too few cells, use NA
        stage_mat[, ct] <- NA
        next
      }
      
      # Calculate mean expression in window for this cell type
      in_window_mean <- Matrix::rowMeans(expr_mat[, window_ct_cells, drop = FALSE])
      
      # Calculate mean expression outside window for this cell type
      out_window_ct_cells <- setdiff(ct_cells, window_ct_cells)
      if (length(out_window_ct_cells) > 0) {
        out_window_mean <- Matrix::rowMeans(expr_mat[, out_window_ct_cells, drop = FALSE])
      } else {
        out_window_mean <- rep(0, length(gene_names))
      }
      
      # Stage-specific score: difference between in and out of window
      stage_score <- in_window_mean - out_window_mean
      stage_mat[, ct] <- stage_score
    }
    
    stage_specificity[[stage_name]] <- stage_mat
  }
  
  return(stage_specificity)
}


#' Test trajectory-stage-specific associations with trait
#'
#' Tests whether gene expression at specific trajectory stages drives
#' cell type-trait associations.
#'
#' @param stage_specificity Output from calc_trajectory_specificity.
#' @param magma A data.frame or file path to MAGMA output.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#'
#' @return A data.frame with stage-specific association results.
#' @export
test_trajectory_associations <- function(stage_specificity,
                                        magma,
                                        magma_gene_col = "GENE",
                                        magma_z_col = "ZSTAT") {
  
  pvalue <- FDR <- NULL # due to non-standard evaluation notes
  
  # Load MAGMA data
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  magma <- magma[, c(magma_gene_col, magma_z_col), with = FALSE]
  names(magma) <- c("gene", "zstat")
  
  # Test each stage
  results <- rbindlist(lapply(names(stage_specificity), function(stage) {
    stage_mat <- stage_specificity[[stage]]
    
    # Convert to data.table format
    stage_dt <- as.data.table(stage_mat, keep.rownames = TRUE)
    setnames(stage_dt, "rn", "gene")
    stage_dt <- melt(stage_dt, id.vars = "gene", 
                    variable.name = "cell_type", value.name = "stage_score")
    
    # Merge with MAGMA
    dt <- merge(stage_dt, magma, by = "gene")
    dt <- dt[stats::complete.cases(dt)]
    
    # Get cell types
    cts <- unique(dt$cell_type)
    
    # Test association for each cell type in this stage
    stage_results <- rbindlist(lapply(cts, function(ct) {
      ct_data <- dt[cell_type == ct]
      
      if (nrow(ct_data) < 10) {
        return(data.table(
          stage = stage,
          cell_type = ct,
          beta = NA,
          pvalue = NA,
          n_genes = nrow(ct_data)
        ))
      }
      
      # Fit model
      lm_stage <- stats::lm(zstat ~ stage_score, data = ct_data)
      lm_summ <- summary(lm_stage)$coefficients
      
      beta <- lm_summ[2, 1]
      pval <- lm_summ[2, 4]
      
      # One-sided test
      pval_onesided <- if (beta > 0) pval / 2 else 1 - pval / 2
      
      return(data.table(
        stage = stage,
        cell_type = ct,
        beta = beta,
        pvalue = pval_onesided,
        n_genes = nrow(ct_data)
      ))
    }))
    
    return(stage_results)
  }))
  
  results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  results <- results[order(pvalue)]
  
  return(results)
}


#' Identify stage-specific driver genes
#'
#' Finds genes that are specifically influential at certain trajectory stages.
#'
#' @param cell_type Character string for cell type.
#' @param stage Character string for trajectory stage (e.g., "stage1", "stage2").
#' @param stage_specificity Output from calc_trajectory_specificity.
#' @param magma A data.frame or file path to MAGMA output.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#'
#' @return A data.frame with stage-specific influential genes.
#' @export
find_stage_specific_drivers <- function(cell_type,
                                       stage,
                                       stage_specificity,
                                       magma,
                                       magma_gene_col = "GENE",
                                       magma_z_col = "ZSTAT") {
  
  dfbetas <- is_influential <- NULL # due to non-standard evaluation
  
  if (!stage %in% names(stage_specificity)) {
    stop("Stage '", stage, "' not found in stage_specificity")
  }
  
  stage_mat <- stage_specificity[[stage]]
  
  if (!cell_type %in% colnames(stage_mat)) {
    stop("Cell type '", cell_type, "' not found in stage matrix")
  }
  
  # Load MAGMA data
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  magma <- magma[, c(magma_gene_col, magma_z_col), with = FALSE]
  names(magma) <- c("gene", "zstat")
  
  # Prepare data
  dt <- data.table(
    gene = rownames(stage_mat),
    stage_score = stage_mat[, cell_type]
  )
  
  dt <- merge(dt, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  if (nrow(dt) < 10) {
    warning("Too few genes for analysis")
    return(data.table())
  }
  
  # Fit model and calculate DFBETAS
  lm_out <- stats::lm(zstat ~ stage_score, data = dt)
  dt[, dfbetas := stats::dfbetas(lm_out)[, 2]]
  
  # Mark influential
  thresh <- 2 / sqrt(nrow(dt))
  dt[, is_influential := dfbetas > thresh]
  
  dt <- dt[order(-dfbetas)]
  
  return(dt)
}


#' Compare associations across trajectory stages
#'
#' Tests whether association strength varies significantly across
#' different trajectory stages.
#'
#' @param cell_type Character string for cell type.
#' @param stage_specificity Output from calc_trajectory_specificity.
#' @param magma A data.frame or file path to MAGMA output.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#'
#' @return A data.frame comparing associations across stages.
#' @export
compare_stage_associations <- function(cell_type,
                                      stage_specificity,
                                      magma,
                                      magma_gene_col = "GENE",
                                      magma_z_col = "ZSTAT") {
  
  # Load MAGMA
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  magma <- magma[, c(magma_gene_col, magma_z_col), with = FALSE]
  names(magma) <- c("gene", "zstat")
  
  # Get betas for each stage
  stages <- names(stage_specificity)
  stage_betas <- numeric(length(stages))
  stage_ses <- numeric(length(stages))
  stage_pvals <- numeric(length(stages))
  stage_ngenes <- numeric(length(stages))
  
  for (i in seq_along(stages)) {
    stage <- stages[i]
    stage_mat <- stage_specificity[[stage]]
    
    if (!cell_type %in% colnames(stage_mat)) {
      stage_betas[i] <- NA
      stage_ses[i] <- NA
      stage_pvals[i] <- NA
      stage_ngenes[i] <- 0
      next
    }
    
    dt <- data.table(
      gene = rownames(stage_mat),
      stage_score = stage_mat[, cell_type]
    )
    
    dt <- merge(dt, magma, by = "gene")
    dt <- dt[stats::complete.cases(dt)]
    
    if (nrow(dt) < 10) {
      stage_betas[i] <- NA
      stage_ses[i] <- NA
      stage_pvals[i] <- NA
      stage_ngenes[i] <- nrow(dt)
      next
    }
    
    lm_stage <- stats::lm(zstat ~ stage_score, data = dt)
    lm_summ <- summary(lm_stage)$coefficients
    
    stage_betas[i] <- lm_summ[2, 1]
    stage_ses[i] <- lm_summ[2, 2]
    stage_pvals[i] <- lm_summ[2, 4]
    stage_ngenes[i] <- nrow(dt)
  }
  
  # Create result
  result <- data.table(
    cell_type = cell_type,
    stage = stages,
    beta = stage_betas,
    se = stage_ses,
    pvalue = stage_pvals,
    n_genes = stage_ngenes
  )
  
  # Test for heterogeneity across stages
  # Use Cochran's Q test if we have multiple valid estimates
  valid_stages <- which(!is.na(result$beta))
  
  if (length(valid_stages) >= 2) {
    betas <- result$beta[valid_stages]
    ses <- result$se[valid_stages]
    
    # Cochran's Q
    weights <- 1 / (ses^2)
    weighted_mean <- sum(betas * weights) / sum(weights)
    Q <- sum(weights * (betas - weighted_mean)^2)
    df <- length(valid_stages) - 1
    Q_pval <- 1 - stats::pchisq(Q, df)
    
    result$heterogeneity_Q <- Q
    result$heterogeneity_pval <- Q_pval
    result$heterogeneity_df <- df
  }
  
  return(result)
}


#' Fit dynamic expression model with time-varying effects
#'
#' Models how gene expression changes over pseudotime using smoothing splines
#' or GAM, then tests association with trait.
#'
#' @param sce A SingleCellExperiment object.
#' @param pseudotime_col Column name for pseudotime values.
#' @param ct_label_col Column name for cell type labels.
#' @param magma A data.frame or file path to MAGMA output.
#' @param assay_name Assay name. Default: 'logcounts'.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#' @param smooth_method Smoothing method: "loess" or "gam". Default: "loess".
#'
#' @return A data.frame with dynamic association results.
#' @export
test_dynamic_associations <- function(sce,
                                     pseudotime_col,
                                     ct_label_col,
                                     magma,
                                     assay_name = "logcounts",
                                     magma_gene_col = "GENE",
                                     magma_z_col = "ZSTAT",
                                     smooth_method = "loess") {
  
  # Load MAGMA
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # Extract data
  expr_mat <- SummarizedExperiment::assay(sce, assay_name)
  cell_meta <- SummarizedExperiment::colData(sce)
  
  pseudotime <- cell_meta[[pseudotime_col]]
  cell_types <- cell_meta[[ct_label_col]]
  
  # Filter valid cells
  valid_cells <- !is.na(pseudotime)
  expr_mat <- expr_mat[, valid_cells]
  pseudotime <- pseudotime[valid_cells]
  cell_types <- cell_types[valid_cells]
  
  unique_cts <- unique(cell_types)
  genes_to_test <- intersect(rownames(expr_mat), magma[[magma_gene_col]])
  
  # For each cell type, fit dynamic models
  results <- rbindlist(lapply(unique_cts, function(ct) {
    ct_cells <- which(cell_types == ct)
    
    if (length(ct_cells) < 30) {
      return(data.table(
        cell_type = ct,
        dynamic_score = NA,
        pvalue = NA
      ))
    }
    
    ct_expr <- expr_mat[genes_to_test, ct_cells]
    ct_pt <- pseudotime[ct_cells]
    
    # Calculate dynamic variance (how much expression varies over time)
    dynamic_scores <- apply(ct_expr, 1, function(gene_expr) {
      if (smooth_method == "loess") {
        fit <- stats::loess(gene_expr ~ ct_pt, span = 0.3)
        fitted_vals <- stats::fitted(fit)
      } else if (smooth_method == "gam" && requireNamespace("mgcv", quietly = TRUE)) {
        fit <- mgcv::gam(gene_expr ~ s(ct_pt))
        fitted_vals <- stats::fitted(fit)
      } else {
        # Fallback to linear
        fit <- stats::lm(gene_expr ~ ct_pt)
        fitted_vals <- stats::fitted(fit)
      }
      
      # Dynamic score: variance of fitted trajectory
      return(stats::var(fitted_vals))
    })
    
    # Test association between dynamic scores and GWAS
    dt <- data.table(
      gene = genes_to_test,
      dynamic_score = dynamic_scores
    )
    
    magma_sub <- magma[, c(magma_gene_col, magma_z_col), with = FALSE]
    names(magma_sub) <- c("gene", "zstat")
    dt <- merge(dt, magma_sub, by = "gene")
    dt <- dt[stats::complete.cases(dt)]
    
    if (nrow(dt) < 10) {
      return(data.table(
        cell_type = ct,
        beta = NA,
        pvalue = NA,
        n_genes = nrow(dt)
      ))
    }
    
    # Test association
    lm_dyn <- stats::lm(zstat ~ dynamic_score, data = dt)
    lm_summ <- summary(lm_dyn)$coefficients
    
    beta <- lm_summ[2, 1]
    pval <- lm_summ[2, 4]
    pval_onesided <- if (beta > 0) pval / 2 else 1 - pval / 2
    
    return(data.table(
      cell_type = ct,
      beta = beta,
      pvalue = pval_onesided,
      n_genes = nrow(dt)
    ))
  }))
  
  results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  results <- results[order(pvalue)]
  
  return(results)
}
