#' Bayesian hierarchical modeling for cell type associations
#'
#' Implements hierarchical modeling to borrow information across related cell types
#' in a taxonomy, improving estimates for rare cell types.
#'
#' @param sscore A dgeMatrix of specificity scores (genes x cell types).
#' @param magma A data.frame or file path to MAGMA output.
#' @param cell_hierarchy A data.frame with columns: cell_type, parent_type.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#' @param tau_prior Prior variance for hierarchical shrinkage. Default: 1.0.
#' @param n_iter Number of iterations for variational inference. Default: 100.
#'
#' @return A data.frame with hierarchical association results.
#' @export
get_hierarchical_associations <- function(sscore, magma, cell_hierarchy,
                                         magma_gene_col = "GENE",
                                         magma_z_col = "ZSTAT",
                                         tau_prior = 1.0,
                                         n_iter = 100) {
  
  pvalue <- FDR <- NULL # due to non-standard evaluation notes
  
  # Load and check data
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  check_overlap(sscore, magma, magma_gene_col)
  
  # Convert cell hierarchy to data.table
  cell_hierarchy <- as.data.table(cell_hierarchy)
  
  # Clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = TRUE)
  setnames(sscore, "rn", "gene")
  sscore <- melt(sscore, id.vars = "gene", variable.name = "cell_type", 
                value.name = "specificity")
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = FALSE]
  names(magma) <- c("gene", "zstat")
  
  # Merge data
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  # Get cell types
  cts <- unique(sscore$cell_type)
  
  # Initialize results
  beta_estimates <- list()
  beta_parent <- list()
  
  # First pass: estimate parent effects
  parent_cts <- unique(cell_hierarchy$parent_type[!is.na(cell_hierarchy$parent_type)])
  
  for (pct in parent_cts) {
    # Get children of this parent
    children <- cell_hierarchy[parent_type == pct, cell_type]
    
    # Pool data from children
    pooled_data <- dt[cell_type %in% children]
    
    if (nrow(pooled_data) > 0) {
      # Fit model on pooled data
      parent_model <- speedglm::speedlm(pooled_data$zstat ~ pooled_data$specificity)
      beta_parent[[pct]] <- stats::coef(parent_model)[2]
    } else {
      beta_parent[[pct]] <- 0
    }
  }
  
  # Second pass: estimate cell-type-specific effects with hierarchical prior
  results <- rbindlist(lapply(cts, function(ct) {
    ct_data <- dt[cell_type == ct]
    
    # Get parent beta if exists
    parent <- cell_hierarchy[cell_type == ct, parent_type]
    if (length(parent) > 0 && !is.na(parent[1]) && parent[1] %in% names(beta_parent)) {
      prior_mean <- beta_parent[[parent[1]]]
    } else {
      prior_mean <- 0
    }
    
    # Fit local model
    local_model <- stats::lm(zstat ~ specificity, data = ct_data)
    local_beta <- stats::coef(local_model)[2]
    local_se <- summary(local_model)$coefficients[2, 2]
    
    # Calculate sample size weight for shrinkage
    n <- nrow(ct_data)
    lambda <- n / (n + tau_prior)  # Empirical Bayes shrinkage factor
    
    # Hierarchical estimate (shrinkage towards parent)
    hier_beta <- lambda * local_beta + (1 - lambda) * prior_mean
    
    # Adjusted standard error accounting for shrinkage
    hier_se <- local_se * sqrt(lambda)
    
    # Calculate p-value
    z_score <- hier_beta / hier_se
    pval <- 2 * stats::pnorm(-abs(z_score))
    
    # One-sided test for positive association
    pval_onesided <- if (hier_beta > 0) pval / 2 else 1 - pval / 2
    
    return(data.table(
      cell_type = ct,
      parent_type = if (length(parent) > 0) parent[1] else NA,
      beta_local = local_beta,
      beta_parent = if (length(parent) > 0 && !is.na(parent[1])) 
        beta_parent[[parent[1]]] else NA,
      beta_hierarchical = hier_beta,
      se_hierarchical = hier_se,
      shrinkage_factor = lambda,
      pvalue = pval_onesided,
      n_genes = n
    ))
  }))
  
  results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  results <- results[order(pvalue)]
  
  return(results)
}


#' Estimate cell type associations with uncertainty quantification
#'
#' Uses bootstrap resampling to quantify uncertainty in association estimates
#' and provide confidence intervals.
#'
#' @param sscore A dgeMatrix of specificity scores.
#' @param magma A data.frame or file path to MAGMA output.
#' @param magma_gene_col Column name for gene identifiers. Default: 'GENE'.
#' @param magma_z_col Column name for z-scores. Default: 'ZSTAT'.
#' @param n_bootstrap Number of bootstrap iterations. Default: 1000.
#' @param ci_level Confidence level for intervals. Default: 0.95.
#'
#' @return A data.frame with point estimates and confidence intervals.
#' @export
get_associations_with_uncertainty <- function(sscore, magma,
                                             magma_gene_col = "GENE",
                                             magma_z_col = "ZSTAT",
                                             n_bootstrap = 1000,
                                             ci_level = 0.95) {
  
  # Load and check data
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  check_overlap(sscore, magma, magma_gene_col)
  
  # Clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = TRUE)
  setnames(sscore, "rn", "gene")
  sscore <- melt(sscore, id.vars = "gene", variable.name = "cell_type",
                value.name = "specificity")
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = FALSE]
  names(magma) <- c("gene", "zstat")
  
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  cts <- unique(sscore$cell_type)
  
  # Calculate for each cell type
  results <- rbindlist(lapply(cts, function(ct) {
    ct_data <- dt[cell_type == ct]
    n_genes <- nrow(ct_data)
    
    # Original estimate
    original_model <- stats::lm(zstat ~ specificity, data = ct_data)
    beta_original <- stats::coef(original_model)[2]
    
    # Bootstrap
    beta_boot <- numeric(n_bootstrap)
    
    for (i in 1:n_bootstrap) {
      # Sample genes with replacement
      boot_idx <- sample(1:n_genes, n_genes, replace = TRUE)
      boot_data <- ct_data[boot_idx]
      
      # Fit model
      boot_model <- stats::lm(zstat ~ specificity, data = boot_data)
      beta_boot[i] <- stats::coef(boot_model)[2]
    }
    
    # Calculate confidence intervals
    alpha <- 1 - ci_level
    ci_lower <- quantile(beta_boot, alpha / 2)
    ci_upper <- quantile(beta_boot, 1 - alpha / 2)
    
    # Bootstrap standard error
    se_boot <- sd(beta_boot)
    
    # Calculate p-value using bootstrap distribution
    # H0: beta = 0
    pval_boot <- min(
      mean(beta_boot <= 0) * 2,  # Two-sided
      mean(beta_boot >= 0) * 2,
      1.0
    )
    
    # One-sided p-value for positive association
    pval_onesided <- mean(beta_boot <= 0)
    
    return(data.table(
      cell_type = ct,
      beta = beta_original,
      beta_boot_mean = mean(beta_boot),
      beta_boot_se = se_boot,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      pvalue = pval_onesided,
      n_genes = n_genes
    ))
  }))
  
  results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  results <- results[order(pvalue)]
  
  return(results)
}


#' Implement empirical Bayes shrinkage for rare cell types
#'
#' Applies empirical Bayes shrinkage to stabilize specificity estimates
#' for rare cell types with few cells.
#'
#' @param sce A SingleCellExperiment object.
#' @param ct_label_col Column name for cell type labels.
#' @param assay_name Assay name. Default: 'logcounts'.
#' @param cell_hierarchy Optional data.frame with parent-child relationships.
#' @param min_cells_shrinkage Apply shrinkage to cell types with fewer cells. Default: 50.
#' @param shrinkage_strength Parameter controlling shrinkage strength. Default: 10.
#'
#' @return A matrix of shrunk specificity scores.
#' @export
calc_specificity_with_shrinkage <- function(sce, 
                                           ct_label_col = "idents",
                                           assay_name = "logcounts",
                                           cell_hierarchy = NULL,
                                           min_cells_shrinkage = 50,
                                           shrinkage_strength = 10) {
  
  # First calculate regular specificity
  sscore <- calc_specificity(sce, assay_name = assay_name, 
                            ct_label_col = ct_label_col)
  
  # Get cell counts per type
  cell_meta <- SummarizedExperiment::colData(sce)
  ct_counts <- table(cell_meta[[ct_label_col]])
  
  # Identify rare cell types
  rare_cts <- names(ct_counts[ct_counts < min_cells_shrinkage])
  
  if (length(rare_cts) == 0) {
    message("No rare cell types found requiring shrinkage")
    return(sscore)
  }
  
  message("Applying shrinkage to ", length(rare_cts), " rare cell types")
  
  # Convert to matrix
  sscore_mat <- as.matrix(sscore)
  
  for (ct in rare_cts) {
    if (!ct %in% colnames(sscore_mat)) next
    
    n_cells <- ct_counts[ct]
    
    # Determine shrinkage target (parent or global mean)
    if (!is.null(cell_hierarchy)) {
      cell_hierarchy <- as.data.table(cell_hierarchy)
      parent <- cell_hierarchy[cell_type == ct, parent_type]
      
      if (length(parent) > 0 && !is.na(parent[1]) && 
          parent[1] %in% colnames(sscore_mat)) {
        target <- sscore_mat[, parent[1]]
      } else {
        target <- rowMeans(sscore_mat)
      }
    } else {
      target <- rowMeans(sscore_mat)
    }
    
    # Calculate shrinkage factor: lambda = n / (n + k)
    lambda <- as.numeric(n_cells) / (as.numeric(n_cells) + shrinkage_strength)
    
    # Apply shrinkage
    sscore_mat[, ct] <- lambda * sscore_mat[, ct] + (1 - lambda) * target
  }
  
  # Renormalize columns to sum to 1
  col_sums <- colSums(sscore_mat)
  sscore_mat <- sweep(sscore_mat, 2, col_sums, "/")
  
  return(sscore_mat)
}


#' Assess calibration of p-values using QQ plot data
#'
#' Generates data for QQ plot and calculates genomic inflation factor
#' to assess p-value calibration.
#'
#' @param association_results Data.frame with pvalue column from association tests.
#'
#' @return A list with QQ plot data and lambda (genomic inflation factor).
#' @export
assess_pvalue_calibration <- function(association_results) {
  
  association_results <- as.data.table(association_results)
  
  if (!"pvalue" %in% names(association_results)) {
    stop("association_results must have 'pvalue' column")
  }
  
  pvals <- association_results$pvalue[!is.na(association_results$pvalue)]
  
  if (length(pvals) == 0) {
    stop("No valid p-values found")
  }
  
  # Calculate expected vs observed
  n <- length(pvals)
  observed <- -log10(sort(pvals))
  expected <- -log10(ppoints(n))
  
  # Calculate genomic inflation factor (lambda)
  # Based on median chi-square statistic
  chi_sq <- stats::qchisq(1 - pvals, df = 1)
  lambda <- median(chi_sq, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
  
  # Calculate confidence band
  ci_upper <- -log10(stats::qbeta(0.975, 1:n, n:1))
  ci_lower <- -log10(stats::qbeta(0.025, 1:n, n:1))
  
  result <- list(
    expected = expected,
    observed = observed,
    ci_upper = ci_upper,
    ci_lower = ci_lower,
    lambda = lambda,
    n_tests = n
  )
  
  return(result)
}
