#' Gene Expression Correlations Amongst Tissue Samples.
#'
#'
#' A package that takes two gene data files and conducts a biweight Midcorrelation analysis.
#' The function will output two matrices containing the biweight midcorrelation correlation values and P-values.
#'
#'@param tissue_1_matrix A file containing gene-donor data of the first tissue
#'@param tissue_2_matrix A file containing gene-donor data of the second tissue
#'@return Biweight Mid correlations and P-values of gene comparisons across both tissue sets
#'@export
#'
#'
tissue_pair_gene_expression <- function(tissue_1_matrix,tissue_2_matrix) {

  ## Test for correlations among genes between tissues, in a single donor

  library(WGCNA)
  gene_biweight_midcorrelation <- bicor(t(as.matrix(tissue_1_matrix)),t(as.matrix(tissue_2_matrix)))
  pvalues_biweight <- bicorAndPvalue(t(as.matrix(tissue_1_matrix)),t(as.matrix(tissue_2_matrix)))$p
  assign("gene_biweight_midcorrelation",gene_biweight_midcorrelation, envir = globalenv())
  assign("pvalues_biweight",pvalues_biweight, envir = globalenv())
}
