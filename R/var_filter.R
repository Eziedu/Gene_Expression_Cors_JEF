#' Variance Filter Function
#'
#' Functions that takes tissue matrices and variance percentage and outputs matrices of the most variable genes for each tissue set.
#' This function is primarily utilized with inputted tissue sets with matching donors and genes.
#' Therefore the var_filter function serves as the supplement after the use of the function `gene_names` which outputs matching donor/gene tissue matrices.
#'@details var_filter details
#'@param tissue_matrix_1 A file containing gene-donor data of the first tissue. Recommend using 'final_tissue_1.matrix' acquired from the use of gene_names()
#'@param tissue_matrix_2 A file containing gene-donor data of the second tissue. Recommend using'final_tissue_2.matrix' acquired from the use of gene_names()
#'@param percentage_variance A decimal entry representing the variance percentage for gene analysis.
#'@return Two matricies of the most variable genes according to variance percentage.
#'@export

var_filter <- function(tissue_matrix_1,tissue_matrix_2, percentage_variance){

  ## Solve for the variance within genes across donors
  library(parallel)
  gene_expression_variance <- mcmapply(function(x) var(as.numeric(tissue_matrix_1[x,]),as.numeric(tissue_matrix_2[x,])),1:nrow(tissue_matrix_1))

  ## Filter for the most variable genes (based on input variance_percentage)
  ind.variance <- which(gene_expression_variance >= quantile(gene_expression_variance, probs = 1 - percentage_variance))

  variable.tissue_1.genes <- final_tissue_1.matrix[ind.variance,]
  variable.tissue_2.genes <- final_tissue_2.matrix[ind.variance,]

  assign("variable.tissue_1.genes",variable.tissue_1.genes, envir = globalenv())
  assign("variable.tissue_2.genes",variable.tissue_2.genes, envir = globalenv())
}
