#' Compute bulk-seismic specificity score for each gene and disease subtype.
#'
#' This function adapts the seismic algorithm for bulk RNA-seq data instead of 
#' single-cell RNA-seq data. It computes specificity scores for genes across 
#' disease subtypes (patient groups) rather than cell types. The key difference 
#' from the original seismic algorithm is that the expression ratio component 
#' (r_i^(c)) is removed, as bulk RNA-seq data typically has non-zero expression 
#' for most genes.
#'
#' @param expression_matrix A numeric matrix of gene expression values where rows 
#' are genes and columns are patient samples. Row names should be gene identifiers 
#' (matching those in your gene-trait association data), and column names should 
#' be sample identifiers.
#' @param subtype_labels A character vector of the same length as the number of 
#' columns in expression_matrix, indicating the disease subtype (or patient group) 
#' for each sample. The names or order should correspond to the column names/order 
#' of expression_matrix.
#' @param min_uniq_subtype The minimum number of unique disease subtypes needed 
#' in the data to proceed with specificity score calculations. Default: 2.
#' @param min_subtype_size Filter for subtypes to be included in the specificity 
#' score calculation based on the number of samples in that subtype. Default: 3.
#' @param min_avg_exp_subtype Filter for genes to be included in the specificity 
#' score calculation based on the mean of subtype-averaged expressions (mean of 
#' mean subtype expression). Default: 0.1.
#'
#' @return A matrix of specificity scores (genes as rows, disease subtypes as columns).
#' @export
calc_bulk_specificity <- function(expression_matrix, subtype_labels,
                                   min_uniq_subtype = 2,
                                   min_subtype_size = 3,
                                   min_avg_exp_subtype = 0.1) {
  subtype <- N <- ave_exp_subtype <- NULL # due to non-standard evaluation notes in R CMD check

  # data formatting checks
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data.frame with genes as rows and samples as columns.")
  }
  
  # convert to matrix if data.frame
  if (is.data.frame(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }

  if (length(subtype_labels) != ncol(expression_matrix)) {
    stop("Length of subtype_labels must match the number of columns in expression_matrix.")
  }

  if (is.null(rownames(expression_matrix))) {
    stop("expression_matrix must have row names (gene identifiers).")
  }

  # ensure column names exist for samples
  if (is.null(colnames(expression_matrix))) {
    colnames(expression_matrix) <- paste0("sample.", 1:ncol(expression_matrix))
  }

  # create subtype grouping data.table
  subtype_groups <- data.table(
    sample = colnames(expression_matrix),
    subtype = subtype_labels,
    key = "subtype"
  )

  # check that there are at least a few different subtypes
  subtype_groups_n <- subtype_groups[, .N, by = subtype]

  if (nrow(subtype_groups_n) < min_uniq_subtype) {
    stop("There are fewer than ", min_uniq_subtype, " unique subtypes in the data.
         Decrease the min_uniq_subtype threshold or provide data with more unique subtypes.")
  }

  # filter out subtypes that do not have minimum # of samples
  subtype_groups_n <- subtype_groups_n[N >= min_subtype_size]
  subtype_groups <- subtype_groups[subtype %in% subtype_groups_n$subtype]
  expression_matrix <- expression_matrix[, subtype_groups$sample]

  if (nrow(subtype_groups_n) < min_uniq_subtype) {
    stop("There are fewer than ", min_uniq_subtype, " unique subtypes after
          filtering out subtypes with fewer than ", min_subtype_size, " samples.
         Decrease the min_uniq_subtype or min_subtype_size thresholds or provide
         data with more unique subtypes / larger numbers of samples.")
  }

  # calculate mean gene expression per subtype
  # create factor matrix for aggregation
  factor_mat <- Matrix::fac2sparse(factor(subtype_groups$subtype, levels = unique(subtype_groups$subtype)))
  
  # convert expression matrix to sparse if beneficial
  if (!inherits(expression_matrix, "sparseMatrix")) {
    expression_matrix_sparse <- Matrix::Matrix(expression_matrix, sparse = TRUE)
  } else {
    expression_matrix_sparse <- expression_matrix
  }
  
  sum_mat <- Matrix::t(expression_matrix_sparse %*% Matrix::t(factor_mat))
  mean_mat <- sum_mat %>%
    sweep_sparse(margin = 1, stats = subtype_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(expression_matrix)) %>%
    magrittr::set_rownames(unique(subtype_groups$subtype))

  # filter out genes whose subtype-averaged expression does not exceed a baseline level
  stats.dt <- data.table(
    gene = rownames(expression_matrix),
    ave_exp_subtype = Matrix::colMeans(mean_mat)
  )
  stats.dt <- stats.dt[ave_exp_subtype >= min_avg_exp_subtype]

  if (nrow(stats.dt) < 100) {
    warning("There are fewer than 100 genes that have at least mean expression ",
            min_avg_exp_subtype, " across subtypes.
            Consider relaxing the threshold or double check the input data.")
  }

  expression_matrix_sparse <- expression_matrix_sparse[stats.dt$gene, ]
  sum_mat <- sum_mat[, stats.dt$gene]
  mean_mat <- mean_mat[, stats.dt$gene]

  # calculate variance of gene expression per subtype
  var_mat <- (Matrix::t(expression_matrix_sparse^2 %*% Matrix::t(factor_mat)) - 
               2 * mean_mat * sum_mat +
               sweep_sparse(mean_mat^2, margin = 1, stats = subtype_groups_n$N, fun = "*")) %>%
    sweep_sparse(margin = 1, stats = subtype_groups_n$N - 1, fun = "/") %>%
    magrittr::set_colnames(stats.dt$gene) %>%
    magrittr::set_rownames(unique(subtype_groups$subtype))

  # calculate an indicator matrix for all other subtypes (out-group)
  out_group_mat <- matrix(1, nrow = dim(subtype_groups_n)[1], ncol = dim(subtype_groups_n)[1]) -
    diag(nrow = dim(subtype_groups_n)[1], ncol = dim(subtype_groups_n)[1])

  # mean for out group per subtype
  out_mean <- (out_group_mat %*% sum_mat) %>%
    sweep_sparse(margin = 1, stats = out_group_mat %*% subtype_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(stats.dt$gene) %>%
    magrittr::set_rownames(unique(subtype_groups$subtype))

  # variance for out groups per subtype
  tot_var <- sweep_sparse(var_mat, margin = 1, stats = subtype_groups_n$N - 1, fun = "*")
  tot_mean_sq <- sweep_sparse(mean_mat^2, margin = 1, stats = subtype_groups_n$N, fun = "*")
  out_variance <- (out_group_mat %*% tot_var + out_group_mat %*% tot_mean_sq -
                    sweep_sparse(out_mean^2, margin = 1, stats = out_group_mat %*% subtype_groups_n$N, fun = "*")) %>%
    sweep_sparse(margin = 1, stats = out_group_mat %*% subtype_groups_n$N - 1, fun = "/") %>%
    magrittr::set_colnames(stats.dt$gene) %>%
    magrittr::set_rownames(unique(subtype_groups$subtype))

  # probability of relatively higher expression for each gene in each subtype (vs other subtypes)
  # This is p_i^(c) in the original seismic algorithm
  rel_exp <- (mean_mat - out_mean) /
    sqrt(sweep_sparse(var_mat, margin = 1, stats = subtype_groups_n$N - 1, fun = "/") +
           sweep_sparse(out_variance, margin = 1, stats = out_group_mat %*% subtype_groups_n$N - 1, fun = "/"))

  # calculate bulk-seismic specificity score
  # Note: Unlike the original seismic, we do NOT multiply by ratio_mat (r_i^(c))
  # because bulk RNA-seq data typically has non-zero expression for most genes
  p_mat <- stats::pnorm(as.matrix(rel_exp))
  
  # Normalize p_i^(c) across all subtypes to get the final specificity score
  spec_score <- sweep_sparse(x = p_mat, margin = 2, stats = Matrix::colSums(p_mat), fun = "/")

  # Return transposed matrix (genes as rows, subtypes as columns)
  return(Matrix::t(spec_score))
}
