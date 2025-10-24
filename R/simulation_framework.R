#' Generate synthetic single-cell data with known ground truth
#'
#' Creates simulated scRNA-seq data with specified cell type hierarchy,
#' gene expression patterns, and known driver genes for benchmarking.
#'
#' @param n_genes Number of genes to simulate. Default: 5000.
#' @param n_cells_per_type Vector of cell counts for each cell type.
#' @param cell_type_names Character vector of cell type names.
#' @param n_driver_genes Number of true driver genes per cell type. Default: 50.
#' @param driver_effect_size Effect size for driver genes. Default: 2.0.
#' @param noise_sd Standard deviation of noise. Default: 0.5.
#' @param seed Random seed for reproducibility. Default: NULL.
#'
#' @return A list with sce object, true drivers, and simulation parameters.
#' @export
simulate_sc_data <- function(n_genes = 5000,
                            n_cells_per_type = c(100, 200, 150, 80),
                            cell_type_names = c("TypeA", "TypeB", "TypeC", "TypeD"),
                            n_driver_genes = 50,
                            driver_effect_size = 2.0,
                            noise_sd = 0.5,
                            seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_cell_types <- length(cell_type_names)
  total_cells <- sum(n_cells_per_type)
  
  # Generate gene names
  gene_names <- paste0("Gene", 1:n_genes)
  
  # Generate baseline expression (log-normal)
  baseline_expr <- matrix(
    stats::rlnorm(n_genes * total_cells, meanlog = 0, sdlog = 1),
    nrow = n_genes, ncol = total_cells
  )
  rownames(baseline_expr) <- gene_names
  
  # Add cell type specific expression patterns
  cell_assignments <- rep(cell_type_names, times = n_cells_per_type)
  
  # Select driver genes for each cell type
  driver_genes <- list()
  for (i in 1:n_cell_types) {
    driver_genes[[cell_type_names[i]]] <- sample(gene_names, n_driver_genes)
  }
  
  # Add driver gene effects
  expr_mat <- baseline_expr
  cell_idx <- 1
  
  for (i in 1:n_cell_types) {
    ct <- cell_type_names[i]
    n_cells <- n_cells_per_type[i]
    drivers <- driver_genes[[ct]]
    
    # Increase expression of driver genes in this cell type
    for (gene in drivers) {
      expr_mat[gene, cell_idx:(cell_idx + n_cells - 1)] <- 
        expr_mat[gene, cell_idx:(cell_idx + n_cells - 1)] * 
        exp(stats::rnorm(n_cells, mean = driver_effect_size, sd = 0.2))
    }
    
    cell_idx <- cell_idx + n_cells
  }
  
  # Add noise
  expr_mat <- expr_mat + matrix(
    stats::rnorm(n_genes * total_cells, mean = 0, sd = noise_sd),
    nrow = n_genes, ncol = total_cells
  )
  
  # Convert to log scale
  expr_mat[expr_mat < 0] <- 0  # Remove negative values
  logcounts <- log1p(expr_mat)
  
  # Create SingleCellExperiment object
  coldata <- data.frame(
    cell_id = paste0("Cell", 1:total_cells),
    cell_type = cell_assignments,
    row.names = paste0("Cell", 1:total_cells)
  )
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = logcounts),
    colData = coldata
  )
  
  result <- list(
    sce = sce,
    driver_genes = driver_genes,
    n_genes = n_genes,
    n_cells_per_type = n_cells_per_type,
    cell_type_names = cell_type_names,
    driver_effect_size = driver_effect_size,
    noise_sd = noise_sd
  )
  
  return(result)
}


#' Generate synthetic GWAS data with known driver genes
#'
#' Creates simulated GWAS summary statistics where trait risk is driven
#' by specified genes with cell-type-specific effects.
#'
#' @param gene_names Character vector of gene names.
#' @param driver_genes Named list of driver genes per cell type.
#' @param specificity_scores Matrix of specificity scores (genes x cell types).
#' @param effect_size Effect size for driver genes. Default: 0.5.
#' @param noise_sd Standard deviation of noise in Z-scores. Default: 1.0.
#' @param seed Random seed. Default: NULL.
#'
#' @return A data.frame with gene-level GWAS results (MAGMA format).
#' @export
simulate_gwas_data <- function(gene_names,
                              driver_genes,
                              specificity_scores,
                              effect_size = 0.5,
                              noise_sd = 1.0,
                              seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_genes <- length(gene_names)
  
  # Initialize Z-scores with noise
  z_scores <- stats::rnorm(n_genes, mean = 0, sd = noise_sd)
  names(z_scores) <- gene_names
  
  # Add signal from driver genes based on specificity
  spec_mat <- as.matrix(specificity_scores)
  
  for (ct in names(driver_genes)) {
    if (ct %in% colnames(spec_mat)) {
      drivers <- driver_genes[[ct]]
      drivers <- drivers[drivers %in% gene_names]
      
      for (gene in drivers) {
        # Effect proportional to specificity
        gene_spec <- spec_mat[gene, ct]
        z_scores[gene] <- z_scores[gene] + effect_size * gene_spec * 
          stats::rnorm(1, mean = 1, sd = 0.1)
      }
    }
  }
  
  # Convert to p-values
  p_values <- 2 * stats::pnorm(-abs(z_scores))
  
  # Create MAGMA-style output
  gwas_data <- data.frame(
    GENE = gene_names,
    ZSTAT = z_scores,
    P = p_values,
    stringsAsFactors = FALSE
  )
  
  return(gwas_data)
}


#' Benchmark seismic performance on simulated data
#'
#' Tests seismic's ability to recover known driver genes and cell types
#' using simulated data with ground truth.
#'
#' @param sim_sc_data Output from simulate_sc_data.
#' @param sim_gwas_data Output from simulate_gwas_data.
#' @param methods Character vector of methods to test. 
#' Default: c("linear", "gam", "hierarchical").
#' @param fdr_threshold FDR threshold for significance. Default: 0.05.
#'
#' @return A list with performance metrics (power, FDR, AUC, etc.).
#' @export
benchmark_seismic <- function(sim_sc_data,
                             sim_gwas_data,
                             methods = c("linear", "gam", "hierarchical"),
                             fdr_threshold = 0.05) {
  
  # Calculate specificity scores
  sscore <- calc_specificity(sim_sc_data$sce, ct_label_col = "cell_type")
  
  # True driver genes
  true_drivers <- sim_sc_data$driver_genes
  
  # Results list
  results <- list()
  
  # Test each method
  for (method in methods) {
    message("Testing method: ", method)
    
    # Get associations based on method
    if (method == "linear") {
      assoc <- get_ct_trait_associations(sscore, sim_gwas_data)
    } else if (method == "gam") {
      assoc <- get_ct_trait_associations_nonlinear(
        sscore, sim_gwas_data, model_type = "gam"
      )
    } else if (method == "hierarchical") {
      # Create simple hierarchy (for demo, assume flat - in real use provide structure)
      cell_hierarchy <- data.frame(
        cell_type = sim_sc_data$cell_type_names,
        parent_type = rep(NA, length(sim_sc_data$cell_type_names))
      )
      assoc <- get_hierarchical_associations(sscore, sim_gwas_data, cell_hierarchy)
    }
    
    # Calculate power: proportion of true driver cell types identified
    sig_cts <- assoc[FDR < fdr_threshold, cell_type]
    true_cts <- names(true_drivers)
    
    true_positives <- sum(sig_cts %in% true_cts)
    false_positives <- sum(!sig_cts %in% true_cts)
    false_negatives <- sum(!true_cts %in% sig_cts)
    true_negatives <- length(unique(assoc$cell_type)) - 
      true_positives - false_positives - false_negatives
    
    power <- true_positives / length(true_cts)
    fdr_empirical <- false_positives / max(length(sig_cts), 1)
    
    # Calculate AUC
    if (requireNamespace("pROC", quietly = TRUE)) {
      true_labels <- ifelse(assoc$cell_type %in% true_cts, 1, 0)
      pred_scores <- -log10(assoc$pvalue + 1e-300)
      auc <- pROC::auc(pROC::roc(true_labels, pred_scores, quiet = TRUE))
    } else {
      auc <- NA
    }
    
    # Test influential genes for a known driver cell type
    if (length(sig_cts) > 0 && sig_cts[1] %in% true_cts) {
      inf_genes <- find_inf_genes(sig_cts[1], sscore, sim_gwas_data)
      
      # Calculate precision and recall for influential genes
      top_n <- min(50, nrow(inf_genes))
      pred_drivers <- inf_genes$gene[1:top_n]
      true_ct_drivers <- true_drivers[[sig_cts[1]]]
      
      gene_tp <- sum(pred_drivers %in% true_ct_drivers)
      gene_precision <- gene_tp / top_n
      gene_recall <- gene_tp / length(true_ct_drivers)
    } else {
      gene_precision <- NA
      gene_recall <- NA
    }
    
    results[[method]] <- list(
      method = method,
      power = power,
      empirical_fdr = fdr_empirical,
      auc = as.numeric(auc),
      true_positives = true_positives,
      false_positives = false_positives,
      false_negatives = false_negatives,
      gene_precision = gene_precision,
      gene_recall = gene_recall
    )
  }
  
  # Combine results
  results_df <- rbindlist(lapply(results, function(x) as.data.table(x)))
  
  return(list(
    performance = results_df,
    detailed_results = results
  ))
}


#' Generate comprehensive simulation suite
#'
#' Creates multiple simulation scenarios varying key parameters for
#' systematic benchmarking.
#'
#' @param n_scenarios Number of simulation scenarios. Default: 10.
#' @param vary_params Parameters to vary: "sample_size", "effect_size", "noise". 
#' Default: all.
#' @param seed Random seed. Default: NULL.
#'
#' @return A list of simulation scenarios.
#' @export
generate_simulation_suite <- function(n_scenarios = 10,
                                     vary_params = c("sample_size", "effect_size", "noise"),
                                     seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  scenarios <- list()
  
  for (i in 1:n_scenarios) {
    # Base parameters
    n_cells_per_type <- c(100, 200, 150, 80)
    driver_effect_size <- 2.0
    noise_sd <- 0.5
    gwas_effect_size <- 0.5
    
    # Vary parameters
    if ("sample_size" %in% vary_params) {
      multiplier <- runif(1, 0.5, 2.0)
      n_cells_per_type <- round(n_cells_per_type * multiplier)
    }
    
    if ("effect_size" %in% vary_params) {
      driver_effect_size <- runif(1, 1.0, 3.0)
      gwas_effect_size <- runif(1, 0.3, 0.8)
    }
    
    if ("noise" %in% vary_params) {
      noise_sd <- runif(1, 0.2, 1.0)
    }
    
    # Generate data
    sc_data <- simulate_sc_data(
      n_genes = 5000,
      n_cells_per_type = n_cells_per_type,
      driver_effect_size = driver_effect_size,
      noise_sd = noise_sd,
      seed = seed + i
    )
    
    sscore <- calc_specificity(sc_data$sce, ct_label_col = "cell_type")
    
    gwas_data <- simulate_gwas_data(
      gene_names = rownames(sscore),
      driver_genes = sc_data$driver_genes,
      specificity_scores = sscore,
      effect_size = gwas_effect_size,
      seed = seed + i + 1000
    )
    
    scenarios[[i]] <- list(
      scenario_id = i,
      sc_data = sc_data,
      gwas_data = gwas_data,
      sscore = sscore,
      params = list(
        n_cells_per_type = n_cells_per_type,
        driver_effect_size = driver_effect_size,
        noise_sd = noise_sd,
        gwas_effect_size = gwas_effect_size
      )
    )
  }
  
  return(scenarios)
}
