
#' Variance Filter Function
#'
#' Functions that takes tissue matrices and variance percentage and outputs matrices of the most variable genes for each tissue set.
#' @param tissue_matrix_1 A file containing gene-donor data of the first tissue
#' @param tissue_matrix_2 A file containing gene-donor data of the second tissue
#' @param percentage_variance A decimal entry representing the variance percentage for gene analysis.
#' @export
#' @example var_filter(tissue_1, tissue_2, 0.3)
#'
#'
var_filter <- function(tissue_matrix_1,tissue_matrix_2,percentage_variance){

  ## Solve for the variance within genes across donors
  library(iemisc)
  gene_expression_variance <- c()
  for (i in 1:size(as.matrix(final_tissue_1.matrix),1)) {
    gene_expression_variance[i] <- var(as.numeric(final_tissue_1.matrix[i,]),as.numeric(final_tissue_2.matrix[i,]))
  }

  gene_variance_test <- mcmapply(function(x) var(as.numeric(Adipose.matrix[x,]),as.numeric(Liver.matrix[x,])),1:nrow(Adipose.matrix))

  ## Filter for the most variable genes (based on input variance_percentage)
  ind.variance <- which(gene_expression_variance >= quantile(gene_expression_variance, probs = 1 - variance_percentage))

  variable.tissue_1.genes <- final_tissue_1.matrix[ind.variance,]
  variable.tissue_2.genes <- final_tissue_2.matrix[ind.variance,]

  assign("variable.tissue_1.genes",variable.tissue_1.genes, envir = globalenv())
  assign("variable.tissue_2.genes",variable.tissue_2.genes, envir = globalenv())
}
