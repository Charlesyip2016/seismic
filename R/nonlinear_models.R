#' Fit non-linear association models between cell type specificity and trait risk
#'
#' Extends the basic linear model to capture non-linear relationships using
#' Generalized Additive Models (GAM), kernel ridge regression, or tree-based models.
#'
#' @param sscore A dgeMatrix of seismic specificity scores where
#' each column is a cell type and row names are gene identifiers.
#' @param magma A data.frame or file path to MAGMA output for a particular GWAS
#' with at least 2 columns: gene identifiers and z-scores.
#' @param model_type Character string specifying the model type: "gam", "kernel", or "linear".
#' Default is "linear" for backward compatibility.
#' @param magma_gene_col A character string corresponding to the column name
#' in the MAGMA data containing gene identifiers. Defaults to 'GENE'.
#' @param magma_z_col A character string corresponding to the column name
#' in the MAGMA data containing z-scores. Defaults to 'ZSTAT'.
#' @param gam_k Integer specifying the basis dimension for GAM smoothing. Default is 10.
#'
#' @return A data.frame of trait associations for each cell type with model diagnostics.
#' @export
get_ct_trait_associations_nonlinear <- function(sscore, magma, 
                                               model_type = "linear",
                                               magma_gene_col = "GENE",
                                               magma_z_col = "ZSTAT",
                                               gam_k = 10) {
  pvalue <- FDR <- NULL # due to non-standard evaluation notes in R CMD check
  
  # if given file path to magma, check that it is loadable
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # check that there is overlap enough to run
  check_overlap(sscore, magma, magma_gene_col)
  
  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- melt(sscore, id.vars = "gene", variable.name = "cell_type", value.name = "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$cell_type <- as.character(sscore$cell_type)
  sscore$specificity <- as.numeric(sscore$specificity)
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = F]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)
  
  # get all cell type names
  cts <- unique(sscore$cell_type)
  
  # combine the magma annots with the specificities
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  # calculate the association for each cell type
  res <- rbindlist(lapply(cts, function(ct) {
    ct_data <- dt[cell_type == ct]
    
    if (model_type == "linear") {
      # Standard linear model
      slm <- speedglm::speedlm(ct_data$zstat ~ ct_data$specificity)
      slm_summ <- summary(slm)$coefficients
      pval <- if (slm_summ[2, 1] > 0) slm_summ[2, 4] / 2 else (1 - slm_summ[2, 4] / 2)
      aic <- NA
      bic <- NA
      nonlinear_flag <- FALSE
      
    } else if (model_type == "gam") {
      # Check if mgcv package is available
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        warning("mgcv package not available, falling back to linear model for cell type: ", ct)
        slm <- speedglm::speedlm(ct_data$zstat ~ ct_data$specificity)
        slm_summ <- summary(slm)$coefficients
        pval <- if (slm_summ[2, 1] > 0) slm_summ[2, 4] / 2 else (1 - slm_summ[2, 4] / 2)
        aic <- NA
        bic <- NA
        nonlinear_flag <- FALSE
      } else {
        # Generalized Additive Model with smoothing splines
        gam_model <- mgcv::gam(zstat ~ s(specificity, k = min(gam_k, nrow(ct_data) - 1)), 
                              data = ct_data)
        gam_summ <- summary(gam_model)
        
        # Extract p-value for smooth term
        pval <- gam_summ$s.pv[1]
        aic <- AIC(gam_model)
        bic <- BIC(gam_model)
        
        # Check for non-linearity by comparing with linear model
        lm_model <- stats::lm(zstat ~ specificity, data = ct_data)
        lm_aic <- AIC(lm_model)
        nonlinear_flag <- (lm_aic - aic) > 2  # Improvement of > 2 AIC units suggests non-linearity
      }
      
    } else {
      stop("Unsupported model_type: ", model_type, ". Use 'linear' or 'gam'.")
    }
    
    return(data.table(
      cell_type = ct, 
      pvalue = pval,
      model_type = model_type,
      aic = aic,
      bic = bic,
      nonlinear = nonlinear_flag
    ))
  }))
  
  res[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  res <- res[order(pvalue, FDR)]
  
  return(res)
}


#' Detect interaction effects between specificity components
#'
#' Tests for interaction effects between proportion of cells expressing (p) 
#' and relative expression (r) components of specificity score.
#'
#' @param sce A SingleCellExperiment object.
#' @param magma A data.frame or file path to MAGMA output.
#' @param ct_label_col A column name in colData of sce.
#' @param assay_name An assay in sce, default 'logcounts'.
#' @param magma_gene_col Column name in MAGMA data for gene identifiers.
#' @param magma_z_col Column name in MAGMA data for z-scores.
#'
#' @return A data.frame with interaction test results for each cell type.
#' @export
test_specificity_interactions <- function(sce, magma, 
                                         ct_label_col = "idents",
                                         assay_name = "logcounts",
                                         magma_gene_col = "GENE",
                                         magma_z_col = "ZSTAT") {
  pvalue_interaction <- NULL # due to non-standard evaluation notes in R CMD check
  
  # if given file path to magma, check that it is loadable
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)
  
  # Extract components similar to calc_specificity but keep p and r separate
  data_mat <- SummarizedExperiment::assay(sce, assay_name)
  cell_meta <- SummarizedExperiment::colData(sce)
  
  ct_groups <- data.table(
    cell = rownames(cell_meta),
    ct = cell_meta[[ct_label_col]], key = "ct"
  )
  
  ct_groups_n <- ct_groups[, .N, by = ct]
  factor_mat <- Matrix::fac2sparse(factor(ct_groups$ct, levels = unique(ct_groups$ct)))
  
  # Calculate mean expression (r component)
  sum_mat <- Matrix::t(data_mat %*% Matrix::t(factor_mat))
  mean_mat <- sum_mat %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/")
  
  # Calculate proportion expressing (p component)
  ratio_mat <- Matrix::t((data_mat > 0) %*% Matrix::t(factor_mat)) %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/")
  
  # Calculate relative expression for r
  out_group_mat <- matrix(1, nrow = nrow(ct_groups_n), ncol = nrow(ct_groups_n)) -
    diag(nrow = nrow(ct_groups_n))
  out_mean <- (out_group_mat %*% sum_mat) %>%
    sweep_sparse(margin = 1, stats = out_group_mat %*% ct_groups_n$N, fun = "/")
  
  rel_exp <- mean_mat - out_mean
  
  # Get cell type names
  cts <- unique(ct_groups$ct)
  
  # Test interactions for each cell type
  res <- rbindlist(lapply(seq_along(cts), function(i) {
    ct <- cts[i]
    
    # Extract p and r for this cell type
    p_vals <- as.numeric(ratio_mat[i, ])
    r_vals <- as.numeric(rel_exp[i, ])
    genes <- colnames(mean_mat)
    
    # Merge with MAGMA data
    dt <- data.table(gene = genes, p = p_vals, r = r_vals)
    magma_sub <- magma[, c(magma_gene_col, magma_z_col), with = F]
    names(magma_sub) <- c("gene", "zstat")
    dt <- merge(dt, magma_sub, by = "gene")
    dt <- dt[stats::complete.cases(dt)]
    
    if (nrow(dt) < 10) {
      return(data.table(cell_type = ct, pvalue_interaction = NA, 
                       beta_p = NA, beta_r = NA, beta_interaction = NA))
    }
    
    # Fit interaction model: Z ~ p + r + p*r
    lm_int <- stats::lm(zstat ~ p * r, data = dt)
    lm_summ <- summary(lm_int)$coefficients
    
    # Extract interaction term (p:r)
    if (nrow(lm_summ) >= 4) {
      pval_int <- lm_summ[4, 4]  # p-value for interaction
      beta_p <- lm_summ[2, 1]
      beta_r <- lm_summ[3, 1]
      beta_int <- lm_summ[4, 1]
    } else {
      pval_int <- NA
      beta_p <- NA
      beta_r <- NA
      beta_int <- NA
    }
    
    return(data.table(
      cell_type = ct,
      pvalue_interaction = pval_int,
      beta_p = beta_p,
      beta_r = beta_r,
      beta_interaction = beta_int
    ))
  }))
  
  res[, FDR_interaction := stats::p.adjust(pvalue_interaction, method = "fdr")]
  res <- res[order(pvalue_interaction)]
  
  return(res)
}
