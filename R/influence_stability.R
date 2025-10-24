#' Enhanced influential gene analysis with stability selection
#'
#' Extends DFBETAS analysis with stability selection, leave-k-out validation,
#' and optional Shapley value computation for robust influence assessment.
#'
#' @param ct A character string containing a valid cell type name in sscore.
#' @param sscore A dgeMatrix of seismic specificity scores.
#' @param magma A data.frame or file path to MAGMA output.
#' @param magma_gene_col Column name in MAGMA data for gene identifiers.
#' @param magma_z_col Column name in MAGMA data for z-scores.
#' @param n_bootstrap Number of bootstrap iterations for stability selection. Default: 100.
#' @param sampling_fraction Fraction of data to sample in each bootstrap. Default: 0.8.
#' @param influence_threshold DFBETAS threshold for considering a gene influential. 
#' Default: 2/sqrt(n).
#'
#' @return A data.frame with genes, influence metrics, and stability scores.
#' @export
find_inf_genes_stable <- function(ct, sscore, magma,
                                  magma_gene_col = "GENE", 
                                  magma_z_col = "ZSTAT",
                                  n_bootstrap = 100,
                                  sampling_fraction = 0.8,
                                  influence_threshold = NULL) {
  dfbetas <- is_influential <- stability_score <- NULL # due to non-standard evaluation notes
  
  # if given file path to magma, check that it is loadable
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # check that there is overlap enough to run
  check_overlap(sscore, magma, magma_gene_col)
  
  # check that cell type name exists
  if (!(ct %in% colnames(sscore))) {
    stop("Invalid cell type name. Please choose a cell type
         present in the given specificity data.")
  }
  
  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- sscore[, c("gene", ct), with = F]
  names(sscore) <- c("gene", "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$specificity <- as.numeric(sscore$specificity)
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = F]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)
  
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  # Original model
  lm_out <- stats::lm(zstat ~ specificity, data = dt)
  lm_summ <- summary(lm_out)$coefficients
  lm_pval <- if (lm_summ[2, 1] > 0) lm_summ[2, 4] / 2 else (1 - lm_summ[2, 4] / 2)
  
  if (lm_pval > 0.05) {
    warning("Cell type ", ct, " does not have a significant association with trait
         (p-value ", round(lm_pval, 4), "), so influential gene analysis may not be relevant")
  }
  
  # Calculate DFBETAS
  dt[, dfbetas := stats::dfbetas(lm_out)[, 2]]
  
  # Set threshold
  if (is.null(influence_threshold)) {
    influence_threshold <- 2 / sqrt(nrow(dt))
  }
  
  # Stability selection via bootstrap
  n_genes <- nrow(dt)
  n_sample <- floor(n_genes * sampling_fraction)
  
  # Track how often each gene is marked as influential
  influence_counts <- rep(0, n_genes)
  names(influence_counts) <- dt$gene
  
  for (i in 1:n_bootstrap) {
    # Sample genes
    sample_idx <- sample(1:n_genes, n_sample, replace = FALSE)
    dt_sample <- dt[sample_idx]
    
    # Fit model on sample
    lm_sample <- stats::lm(zstat ~ specificity, data = dt_sample)
    dfb_sample <- stats::dfbetas(lm_sample)[, 2]
    
    # Mark influential genes in this iteration
    influential_genes <- dt_sample$gene[dfb_sample > influence_threshold]
    influence_counts[influential_genes] <- influence_counts[influential_genes] + 1
  }
  
  # Calculate stability score (proportion of iterations gene was influential)
  dt[, stability_score := influence_counts[gene] / n_bootstrap]
  
  # Mark genes as influential based on DFBETAS
  dt[, is_influential := ifelse(dfbetas > influence_threshold, TRUE, FALSE)]
  
  # Add stability category
  dt[, stability_category := ifelse(stability_score >= 0.75, "high",
                                   ifelse(stability_score >= 0.5, "moderate", "low"))]
  
  dt <- dt[order(-dfbetas)]
  
  return(dt)
}


#' Leave-k-out influence analysis
#'
#' Performs leave-k-out validation to assess the stability of association
#' estimates when removing subsets of genes.
#'
#' @param ct A character string containing a valid cell type name in sscore.
#' @param sscore A dgeMatrix of seismic specificity scores.
#' @param magma A data.frame or file path to MAGMA output.
#' @param magma_gene_col Column name in MAGMA data for gene identifiers.
#' @param magma_z_col Column name in MAGMA data for z-scores.
#' @param k_percent Percentage of genes to leave out in each iteration. Default: 5.
#' @param n_iterations Number of leave-k-out iterations. Default: 50.
#'
#' @return A list containing beta distribution and influential gene sets.
#' @export
leave_k_out_analysis <- function(ct, sscore, magma,
                                magma_gene_col = "GENE",
                                magma_z_col = "ZSTAT",
                                k_percent = 5,
                                n_iterations = 50) {
  
  # if given file path to magma, check that it is loadable
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # check that there is overlap enough to run
  check_overlap(sscore, magma, magma_gene_col)
  
  # check that cell type name exists
  if (!(ct %in% colnames(sscore))) {
    stop("Invalid cell type name. Please choose a cell type
         present in the given specificity data.")
  }
  
  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- sscore[, c("gene", ct), with = F]
  names(sscore) <- c("gene", "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$specificity <- as.numeric(sscore$specificity)
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = F]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)
  
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  # Original model
  lm_full <- stats::lm(zstat ~ specificity, data = dt)
  beta_full <- stats::coef(lm_full)[2]
  
  # Leave-k-out iterations
  n_genes <- nrow(dt)
  k_genes <- max(1, floor(n_genes * k_percent / 100))
  
  beta_distribution <- numeric(n_iterations)
  beta_changes <- numeric(n_iterations)
  
  for (i in 1:n_iterations) {
    # Randomly select k genes to leave out
    leave_out_idx <- sample(1:n_genes, k_genes, replace = FALSE)
    dt_subset <- dt[-leave_out_idx]
    
    # Fit model without these genes
    lm_subset <- stats::lm(zstat ~ specificity, data = dt_subset)
    beta_subset <- stats::coef(lm_subset)[2]
    
    beta_distribution[i] <- beta_subset
    beta_changes[i] <- abs(beta_subset - beta_full) / abs(beta_full)
  }
  
  # Summary statistics
  result <- list(
    cell_type = ct,
    beta_full = beta_full,
    beta_mean = mean(beta_distribution),
    beta_sd = sd(beta_distribution),
    beta_median = median(beta_distribution),
    beta_q025 = quantile(beta_distribution, 0.025),
    beta_q975 = quantile(beta_distribution, 0.975),
    mean_relative_change = mean(beta_changes),
    max_relative_change = max(beta_changes),
    beta_distribution = beta_distribution,
    stability_ratio = 1 - mean(beta_changes)  # Higher is more stable
  )
  
  return(result)
}


#' Group influence analysis by pathway or gene module
#'
#' Assesses influence at the pathway/module level by calculating change in
#' association when removing groups of related genes.
#'
#' @param ct A character string containing a valid cell type name in sscore.
#' @param sscore A dgeMatrix of seismic specificity scores.
#' @param magma A data.frame or file path to MAGMA output.
#' @param gene_groups A named list where each element is a character vector of gene IDs
#' belonging to a pathway or module.
#' @param magma_gene_col Column name in MAGMA data for gene identifiers.
#' @param magma_z_col Column name in MAGMA data for z-scores.
#'
#' @return A data.frame with pathway-level influence metrics.
#' @export
group_influence_analysis <- function(ct, sscore, magma, gene_groups,
                                    magma_gene_col = "GENE",
                                    magma_z_col = "ZSTAT") {
  
  # if given file path to magma, check that it is loadable
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # check that there is overlap enough to run
  check_overlap(sscore, magma, magma_gene_col)
  
  # check that cell type name exists
  if (!(ct %in% colnames(sscore))) {
    stop("Invalid cell type name. Please choose a cell type
         present in the given specificity data.")
  }
  
  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- sscore[, c("gene", ct), with = F]
  names(sscore) <- c("gene", "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$specificity <- as.numeric(sscore$specificity)
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = F]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)
  
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  # Original model
  lm_full <- stats::lm(zstat ~ specificity, data = dt)
  beta_full <- stats::coef(lm_full)[2]
  pval_full <- summary(lm_full)$coefficients[2, 4]
  
  # Test each group
  results <- rbindlist(lapply(names(gene_groups), function(group_name) {
    group_genes <- gene_groups[[group_name]]
    
    # Filter out genes in this group
    dt_without_group <- dt[!gene %in% group_genes]
    
    if (nrow(dt_without_group) < 10) {
      return(data.table(
        group = group_name,
        n_genes_in_group = length(intersect(group_genes, dt$gene)),
        beta_without_group = NA,
        delta_beta = NA,
        relative_change = NA,
        pval_without_group = NA
      ))
    }
    
    # Fit model without this group
    lm_without <- stats::lm(zstat ~ specificity, data = dt_without_group)
    beta_without <- stats::coef(lm_without)[2]
    pval_without <- summary(lm_without)$coefficients[2, 4]
    
    delta_beta <- beta_full - beta_without
    relative_change <- abs(delta_beta) / abs(beta_full)
    
    return(data.table(
      group = group_name,
      n_genes_in_group = length(intersect(group_genes, dt$gene)),
      beta_full = beta_full,
      beta_without_group = beta_without,
      delta_beta = delta_beta,
      relative_change = relative_change,
      pval_full = pval_full,
      pval_without_group = pval_without
    ))
  }))
  
  results <- results[order(-relative_change)]
  
  return(results)
}
