#' geneExpCor Wrapper Function //WORKING TITLE//
#'
#' This function performs the entire tissue pair correlation analysis in a single function call.
#' It requires 2 gene expression matrices, a variance percentage for filtering, and the desired amount of top correlations to return.
#' @param tissue1 A matrix object of normalized TPM expression data, with genes as rows and sample subject identifiers as columns.
#' @param tissue2 A matrix object of normalized TPM expression data, with genes as rows and sample subject identifiers as columns.
#' @param var_perc A decimal representing the percentage of highest variable genes to keep.
#' @param var_perc2 A decimal representing the percentage of highest most significant pvalues to keep
#' @param var_plot Boolean value for whether you would like histogram visual of gene variances
#' @export
#'
wrapper_func <- function(tissue1, tissue2, var_perc=.05,var_perc2=.0001,var_plot = TRUE) {

  # filter tissue matrices to retain matching genes and subjects
  out.pair <- gene_names(tissue1, tissue2)

  # filter tissue matrices to retain indicated percentage of highest variable genes
  out.var <- var_filter(out.pair$tissue_1, out.pair$tissue_2, var_perc,var_plot)

  # perform biweight midcorrelation on tissue matrices to determine correlation coefficients and respective p-values
  out.corr <- tissue_pair_gene_expression(out.var$tissue_1, out.var$tissue_2,var_perc2)

}
