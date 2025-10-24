#' Perform Mendelian Randomization analysis
#'
#' Conducts two-sample MR using cell-type-specific eQTLs as instruments to
#' infer causal relationships between gene expression and trait.
#'
#' @param eqtl_data A data.frame with columns: gene, cell_type, snp, beta_eqtl, se_eqtl, pval_eqtl.
#' @param gwas_data A data.frame with columns: snp, beta_gwas, se_gwas, pval_gwas.
#' @param cell_type Character string specifying cell type to analyze.
#' @param f_threshold F-statistic threshold for instrument strength. Default: 10.
#' @param pval_threshold P-value threshold for selecting eQTL instruments. Default: 5e-8.
#'
#' @return A data.frame with MR results for each gene.
#' @export
perform_mendelian_randomization <- function(eqtl_data, gwas_data, cell_type,
                                           f_threshold = 10,
                                           pval_threshold = 5e-8) {
  
  eqtl_data <- as.data.table(eqtl_data)
  gwas_data <- as.data.table(gwas_data)
  
  # Filter eQTL data for specified cell type
  eqtl_ct <- eqtl_data[cell_type == !!cell_type]
  
  # Select significant eQTLs as instruments
  eqtl_ct[, f_stat := (beta_eqtl / se_eqtl)^2]
  instruments <- eqtl_ct[pval_eqtl < pval_threshold & f_stat > f_threshold]
  
  if (nrow(instruments) == 0) {
    warning("No valid instruments found for cell type: ", cell_type)
    return(data.table(gene = character(), 
                     mr_beta = numeric(),
                     mr_se = numeric(),
                     mr_pval = numeric(),
                     n_instruments = numeric()))
  }
  
  # Merge with GWAS data
  mr_data <- merge(instruments, gwas_data, by = "snp", suffixes = c("_eqtl", "_gwas"))
  
  # Get unique genes
  genes <- unique(mr_data$gene)
  
  # Perform MR for each gene
  results <- rbindlist(lapply(genes, function(g) {
    gene_data <- mr_data[gene == g]
    
    if (nrow(gene_data) < 3) {
      # Need at least 3 instruments for robust MR
      return(data.table(
        gene = g,
        cell_type = cell_type,
        mr_beta = NA,
        mr_se = NA,
        mr_pval = NA,
        n_instruments = nrow(gene_data),
        method = "insufficient_instruments"
      ))
    }
    
    # IVW (Inverse Variance Weighted) MR
    ivw_result <- perform_ivw_mr(gene_data$beta_eqtl, gene_data$se_eqtl,
                                 gene_data$beta_gwas, gene_data$se_gwas)
    
    return(data.table(
      gene = g,
      cell_type = cell_type,
      mr_beta = ivw_result$beta,
      mr_se = ivw_result$se,
      mr_pval = ivw_result$pval,
      n_instruments = nrow(gene_data),
      method = "IVW"
    ))
  }))
  
  # Calculate FDR
  results[, mr_fdr := stats::p.adjust(mr_pval, method = "fdr")]
  results <- results[order(mr_pval)]
  
  return(results)
}


#' Inverse Variance Weighted MR estimation
#'
#' Internal function to perform IVW MR.
#'
#' @param beta_exposure Vector of SNP-exposure effects.
#' @param se_exposure Vector of SNP-exposure standard errors.
#' @param beta_outcome Vector of SNP-outcome effects.
#' @param se_outcome Vector of SNP-outcome standard errors.
#'
#' @return List with beta, se, and pval.
#' @keywords internal
perform_ivw_mr <- function(beta_exposure, se_exposure, beta_outcome, se_outcome) {
  
  # Calculate ratio estimates (Wald ratios)
  ratio <- beta_outcome / beta_exposure
  
  # Calculate weights (inverse variance of ratio)
  # Var(ratio) ≈ (se_outcome / beta_exposure)^2
  weights <- (beta_exposure / se_outcome)^2
  
  # IVW estimate
  beta_ivw <- sum(ratio * weights) / sum(weights)
  
  # Standard error
  se_ivw <- sqrt(1 / sum(weights))
  
  # P-value
  z_score <- beta_ivw / se_ivw
  pval <- 2 * stats::pnorm(-abs(z_score))
  
  return(list(beta = beta_ivw, se = se_ivw, pval = pval))
}


#' Test for horizontal pleiotropy using MR-Egger
#'
#' Performs MR-Egger regression to detect and correct for horizontal pleiotropy.
#'
#' @param eqtl_data A data.frame with columns: gene, cell_type, snp, beta_eqtl, se_eqtl.
#' @param gwas_data A data.frame with columns: snp, beta_gwas, se_gwas.
#' @param cell_type Character string specifying cell type to analyze.
#' @param pval_threshold P-value threshold for selecting instruments. Default: 5e-6.
#'
#' @return A data.frame with MR-Egger results including intercept test.
#' @export
test_pleiotropy_egger <- function(eqtl_data, gwas_data, cell_type,
                                 pval_threshold = 5e-6) {
  
  eqtl_data <- as.data.table(eqtl_data)
  gwas_data <- as.data.table(gwas_data)
  
  # Filter eQTL data
  eqtl_ct <- eqtl_data[cell_type == !!cell_type & pval_eqtl < pval_threshold]
  
  # Merge with GWAS
  mr_data <- merge(eqtl_ct, gwas_data, by = "snp")
  
  genes <- unique(mr_data$gene)
  
  results <- rbindlist(lapply(genes, function(g) {
    gene_data <- mr_data[gene == g]
    
    if (nrow(gene_data) < 4) {
      return(data.table(
        gene = g,
        cell_type = cell_type,
        egger_beta = NA,
        egger_se = NA,
        egger_pval = NA,
        intercept = NA,
        intercept_pval = NA,
        pleiotropy_detected = NA
      ))
    }
    
    # MR-Egger regression: beta_outcome ~ beta_exposure
    # Weighted by inverse variance
    weights <- 1 / (gene_data$se_gwas^2)
    
    fit <- stats::lm(beta_gwas ~ beta_eqtl, data = gene_data, weights = weights)
    fit_summ <- summary(fit)$coefficients
    
    intercept <- fit_summ[1, 1]
    intercept_pval <- fit_summ[1, 4]
    egger_beta <- fit_summ[2, 1]
    egger_se <- fit_summ[2, 2]
    egger_pval <- fit_summ[2, 4]
    
    # Pleiotropy detected if intercept significantly different from 0
    pleiotropy <- intercept_pval < 0.05
    
    return(data.table(
      gene = g,
      cell_type = cell_type,
      egger_beta = egger_beta,
      egger_se = egger_se,
      egger_pval = egger_pval,
      intercept = intercept,
      intercept_pval = intercept_pval,
      pleiotropy_detected = pleiotropy
    ))
  }))
  
  return(results)
}


#' Perform mediation analysis
#'
#' Tests whether gene expression mediates the relationship between genetic
#' variants and trait using mediation analysis framework.
#'
#' @param snp_exposure Vector of SNP dosages or effects.
#' @param gene_expression Vector of gene expression values (mediator).
#' @param trait_outcome Vector of trait values (outcome).
#' @param covariates Optional matrix of covariates.
#' @param bootstrap Logical, whether to use bootstrap for confidence intervals. Default: TRUE.
#' @param n_bootstrap Number of bootstrap iterations. Default: 1000.
#'
#' @return A list with direct effect, indirect effect, total effect, and proportion mediated.
#' @export
perform_mediation_analysis <- function(snp_exposure, gene_expression, trait_outcome,
                                      covariates = NULL,
                                      bootstrap = TRUE,
                                      n_bootstrap = 1000) {
  
  # Create data frame
  if (is.null(covariates)) {
    med_data <- data.frame(
      snp = snp_exposure,
      expr = gene_expression,
      trait = trait_outcome
    )
  } else {
    med_data <- data.frame(
      snp = snp_exposure,
      expr = gene_expression,
      trait = trait_outcome,
      covariates
    )
  }
  
  # Remove missing values
  med_data <- med_data[stats::complete.cases(med_data), ]
  
  if (nrow(med_data) < 50) {
    warning("Insufficient data for mediation analysis (n < 50)")
    return(NULL)
  }
  
  # Path a: SNP -> Expression
  if (is.null(covariates)) {
    model_a <- stats::lm(expr ~ snp, data = med_data)
  } else {
    model_a <- stats::lm(expr ~ ., data = med_data[, c("expr", "snp", 
                                                       names(covariates))])
  }
  alpha <- stats::coef(model_a)["snp"]
  
  # Path b: Expression -> Trait (controlling for SNP)
  if (is.null(covariates)) {
    model_b <- stats::lm(trait ~ snp + expr, data = med_data)
  } else {
    model_b <- stats::lm(trait ~ ., data = med_data[, c("trait", "snp", "expr",
                                                        names(covariates))])
  }
  beta <- stats::coef(model_b)["expr"]
  
  # Path c: Total effect (SNP -> Trait)
  if (is.null(covariates)) {
    model_c <- stats::lm(trait ~ snp, data = med_data)
  } else {
    model_c <- stats::lm(trait ~ ., data = med_data[, c("trait", "snp",
                                                        names(covariates))])
  }
  tau <- stats::coef(model_c)["snp"]
  
  # Path c': Direct effect (SNP -> Trait controlling for Expression)
  tau_prime <- stats::coef(model_b)["snp"]
  
  # Indirect effect (mediated effect)
  indirect <- alpha * beta
  
  # Proportion mediated
  prop_mediated <- indirect / tau
  
  # Bootstrap confidence intervals if requested
  if (bootstrap) {
    boot_indirect <- numeric(n_bootstrap)
    
    for (i in 1:n_bootstrap) {
      # Sample with replacement
      boot_idx <- sample(1:nrow(med_data), nrow(med_data), replace = TRUE)
      boot_data <- med_data[boot_idx, ]
      
      # Refit models
      if (is.null(covariates)) {
        boot_a <- stats::lm(expr ~ snp, data = boot_data)
        boot_b <- stats::lm(trait ~ snp + expr, data = boot_data)
      } else {
        boot_a <- stats::lm(expr ~ ., data = boot_data[, c("expr", "snp",
                                                           names(covariates))])
        boot_b <- stats::lm(trait ~ ., data = boot_data[, c("trait", "snp", "expr",
                                                            names(covariates))])
      }
      
      boot_alpha <- stats::coef(boot_a)["snp"]
      boot_beta <- stats::coef(boot_b)["expr"]
      boot_indirect[i] <- boot_alpha * boot_beta
    }
    
    # Calculate CI
    ci_lower <- quantile(boot_indirect, 0.025)
    ci_upper <- quantile(boot_indirect, 0.975)
    
    # Significance: CI does not include 0
    significant <- (ci_lower > 0 && ci_upper > 0) || (ci_lower < 0 && ci_upper < 0)
    
  } else {
    ci_lower <- NA
    ci_upper <- NA
    significant <- NA
  }
  
  result <- list(
    total_effect = tau,
    direct_effect = tau_prime,
    indirect_effect = indirect,
    proportion_mediated = prop_mediated,
    path_a = alpha,
    path_b = beta,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    significant = significant
  )
  
  return(result)
}


#' Calculate cell-type-specific causal scores
#'
#' Integrates MR results with cell type specificity to prioritize causal
#' gene-cell type pairs.
#'
#' @param mr_results Data.frame from perform_mendelian_randomization.
#' @param specificity_scores Matrix of specificity scores (genes x cell types).
#' @param mr_pval_threshold P-value threshold for significant MR. Default: 0.05.
#'
#' @return A data.frame with integrated causal-specificity scores.
#' @export
calc_causal_specificity_scores <- function(mr_results, specificity_scores,
                                          mr_pval_threshold = 0.05) {
  
  mr_results <- as.data.table(mr_results)
  
  # Filter for significant MR results
  mr_sig <- mr_results[mr_pval < mr_pval_threshold & !is.na(mr_beta)]
  
  if (nrow(mr_sig) == 0) {
    warning("No significant MR results found")
    return(data.table())
  }
  
  # Convert specificity to data.table
  spec_dt <- as.data.table(as.matrix(specificity_scores), keep.rownames = TRUE)
  setnames(spec_dt, "rn", "gene")
  spec_dt <- melt(spec_dt, id.vars = "gene", variable.name = "cell_type",
                 value.name = "specificity")
  
  # Merge MR results with specificity
  integrated <- merge(mr_sig, spec_dt, by = c("gene", "cell_type"))
  
  # Calculate combined causal-specificity score
  # Score = |MR_beta| * specificity * -log10(MR_pval)
  integrated[, causal_spec_score := abs(mr_beta) * specificity * (-log10(mr_pval))]
  
  # Normalize scores
  max_score <- max(integrated$causal_spec_score, na.rm = TRUE)
  if (max_score > 0) {
    integrated[, causal_spec_score := causal_spec_score / max_score]
  }
  
  integrated <- integrated[order(-causal_spec_score)]
  
  return(integrated)
}
