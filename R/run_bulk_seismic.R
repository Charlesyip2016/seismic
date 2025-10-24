#' Run bulk-seismic analysis to identify disease subtype-trait associations
#'
#' This is the main wrapper function for bulk-seismic analysis. It computes 
#' specificity scores for disease subtypes using bulk RNA-seq data, performs 
#' linear regression analysis to identify subtype-trait associations, and 
#' identifies influential genes driving significant associations.
#'
#' @param gene_trait_vector A data.frame or named numeric vector containing 
#' gene-level trait association scores. If a data.frame, it should have at least 
#' two columns: one for gene identifiers and one for association scores (e.g., 
#' MAGMA z-scores, but can also be pQTL, mQTL, or imaging GWAS scores). If a 
#' named numeric vector, names should be gene identifiers.
#' @param expression_matrix A numeric matrix of gene expression values where rows 
#' are genes and columns are patient samples. Row names should be gene identifiers 
#' (matching those in gene_trait_vector).
#' @param subtype_labels A character vector indicating the disease subtype for 
#' each sample in expression_matrix.
#' @param gene_col If gene_trait_vector is a data.frame, the column name containing 
#' gene identifiers. Default: 'GENE' (as in MAGMA output).
#' @param score_col If gene_trait_vector is a data.frame, the column name containing 
#' trait association scores. Default: 'ZSTAT' (as in MAGMA output).
#' @param significance_threshold The p-value threshold for considering a subtype 
#' association as significant. Default: 0.05.
#' @param min_uniq_subtype Minimum number of unique disease subtypes. Default: 2.
#' @param min_subtype_size Minimum number of samples per subtype. Default: 3.
#' @param min_avg_exp_subtype Minimum mean expression threshold for genes. Default: 0.1.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{subtype_associations}{A data.frame of trait associations for each subtype, 
#'         including p-values and FDR-corrected p-values, ordered by significance.}
#'   \item{influential_genes}{A named list where each element corresponds to a 
#'         significant subtype (p < significance_threshold) and contains a data.frame 
#'         of genes with their specificity scores, trait scores, DFBETAS values, 
#'         and a Boolean indicating if the gene is influential.}
#'   \item{specificity_scores}{The matrix of bulk-seismic specificity scores 
#'         (genes as rows, subtypes as columns) used in the analysis.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with MAGMA z-scores
#' result <- run_bulk_seismic(
#'   gene_trait_vector = magma_data,  # data.frame with GENE and ZSTAT columns
#'   expression_matrix = bulk_rnaseq_matrix,  # genes x samples
#'   subtype_labels = patient_subtypes  # vector of subtype labels
#' )
#' 
#' # View significant subtype associations
#' print(result$subtype_associations)
#' 
#' # View influential genes for the most significant subtype
#' top_subtype <- result$subtype_associations$subtype[1]
#' print(result$influential_genes[[top_subtype]])
#' 
#' # Example with custom gene-level scores (e.g., from pQTL)
#' pqtl_scores <- data.frame(gene = gene_names, pqtl_score = scores)
#' result <- run_bulk_seismic(
#'   gene_trait_vector = pqtl_scores,
#'   expression_matrix = bulk_rnaseq_matrix,
#'   subtype_labels = patient_subtypes,
#'   gene_col = 'gene',
#'   score_col = 'pqtl_score'
#' )
#' }
run_bulk_seismic <- function(gene_trait_vector, expression_matrix, subtype_labels,
                             gene_col = "GENE", score_col = "ZSTAT",
                             significance_threshold = 0.05,
                             min_uniq_subtype = 2,
                             min_subtype_size = 3,
                             min_avg_exp_subtype = 0.1) {
  pvalue <- FDR <- dfbetas <- is_influential <- NULL # due to non-standard evaluation notes in R CMD check

  # Step 1: Calculate bulk-seismic specificity scores
  message("Calculating bulk-seismic specificity scores...")
  sscore <- calc_bulk_specificity(
    expression_matrix = expression_matrix,
    subtype_labels = subtype_labels,
    min_uniq_subtype = min_uniq_subtype,
    min_subtype_size = min_subtype_size,
    min_avg_exp_subtype = min_avg_exp_subtype
  )

  # Step 2: Process gene-trait vector
  if (is.vector(gene_trait_vector) && !is.null(names(gene_trait_vector))) {
    # Convert named vector to data.table
    trait_dt <- data.table(
      gene = names(gene_trait_vector),
      score = as.numeric(gene_trait_vector)
    )
  } else if (is.data.frame(gene_trait_vector)) {
    # Extract relevant columns from data.frame
    if (!gene_col %in% names(gene_trait_vector)) {
      stop("Column '", gene_col, "' not found in gene_trait_vector.")
    }
    if (!score_col %in% names(gene_trait_vector)) {
      stop("Column '", score_col, "' not found in gene_trait_vector.")
    }
    trait_dt <- as.data.table(gene_trait_vector)[, c(gene_col, score_col), with = FALSE]
    setnames(trait_dt, c(gene_col, score_col), c("gene", "score"))
  } else {
    stop("gene_trait_vector must be either a named numeric vector or a data.frame.")
  }

  trait_dt$gene <- as.character(trait_dt$gene)
  trait_dt$score <- as.numeric(trait_dt$score)

  # Check for overlap between genes
  gene_overlap <- intersect(rownames(sscore), trait_dt$gene)
  if (length(gene_overlap) < 100) {
    warning("Only ", length(gene_overlap), " genes overlap between expression_matrix and gene_trait_vector. ",
            "Consider checking that gene identifiers match between the two datasets.")
  }
  if (length(gene_overlap) == 0) {
    stop("No overlapping genes found between expression_matrix and gene_trait_vector. ",
         "Please ensure gene identifiers match.")
  }

  # Step 3: Calculate subtype-trait associations via linear regression
  message("Calculating subtype-trait associations...")
  
  # Convert specificity scores to long format
  sscore_dt <- as.data.table(as.matrix(sscore), keep.rownames = TRUE)
  setnames(sscore_dt, "rn", "gene")
  sscore_long <- melt(sscore_dt, id.vars = "gene", variable.name = "subtype", value.name = "specificity")
  sscore_long$gene <- as.character(sscore_long$gene)
  sscore_long$subtype <- as.character(sscore_long$subtype)
  sscore_long$specificity <- as.numeric(sscore_long$specificity)

  # Get all subtype names
  subtypes <- unique(sscore_long$subtype)

  # Merge specificity scores with trait scores
  merged_dt <- merge(sscore_long, trait_dt, by = "gene")
  merged_dt <- merged_dt[stats::complete.cases(merged_dt)]

  # Calculate association for each subtype
  subtype_results <- rbindlist(lapply(subtypes, function(st) {
    subtype_data <- merged_dt[subtype == st]
    
    # Perform linear regression: trait_score ~ specificity
    slm <- speedglm::speedlm(subtype_data$score ~ subtype_data$specificity)
    slm_summ <- summary(slm)$coefficients

    # Extract 1-sided hypothesis test p-value (testing for positive association)
    pval <- if (slm_summ[2, 1] > 0) slm_summ[2, 4] / 2 else (1 - slm_summ[2, 4] / 2)
    
    return(data.table(subtype = st, pvalue = pval))
  }))

  # Apply FDR correction
  subtype_results[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  subtype_results <- subtype_results[order(pvalue, FDR)]

  # Step 4: Identify influential genes for significant subtypes
  message("Identifying influential genes for significant subtypes...")
  
  significant_subtypes <- subtype_results[pvalue < significance_threshold]$subtype
  
  influential_genes_list <- list()
  
  if (length(significant_subtypes) > 0) {
    for (st in significant_subtypes) {
      # Extract data for this subtype
      subtype_data <- merged_dt[subtype == st]
      
      # Perform linear regression
      lm_out <- stats::lm(score ~ specificity, data = subtype_data)
      lm_summ <- summary(lm_out)$coefficients
      lm_pval <- if (lm_summ[2, 1] > 0) lm_summ[2, 4] / 2 else (1 - lm_summ[2, 4] / 2)
      
      # Calculate DFBETAS for influential gene analysis
      subtype_data[, dfbetas := stats::dfbetas(lm_out)[, 2]]
      
      # Threshold for influential genes: 2 / sqrt(n)
      thresh <- 2 / sqrt(nrow(subtype_data))
      
      # Only consider positive dfbetas values (positive relationships)
      subtype_data[, is_influential := ifelse(dfbetas > thresh, TRUE, FALSE)]
      subtype_data <- subtype_data[order(-dfbetas)]
      
      # Store results
      influential_genes_list[[st]] <- subtype_data[, .(gene, specificity, score, dfbetas, is_influential)]
    }
    message("Found ", length(significant_subtypes), " significant subtype(s) at p < ", significance_threshold)
  } else {
    message("No significant subtypes found at p < ", significance_threshold)
  }

  # Return results
  return(list(
    subtype_associations = subtype_results,
    influential_genes = influential_genes_list,
    specificity_scores = sscore
  ))
}
