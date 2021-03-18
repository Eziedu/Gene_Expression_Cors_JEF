#Collective updated Main script
#Workspace edit for push test
#' Gene Expression Correlations Amongst Tissue Samples.
#'
#'
#' A package that takes two gene data files and produces graphics and correlation information based on preferred method and gene variability. Methods include, 'pearson', 'Spearman', and 'Biweight MidCorrelation' where the biweight method is the default parameter. The operation and optimization of the function is still in progress
#'
#'@param tissue_1_matrix A file containing gene-donor data of the first tissue
#'@param tissue_2_matrix A file containing gene-donor data of the second tissue
#'@param variance_percentage A decimal value requesting the percentage of the most variable genes
#'@param correlation_method Type of correlation method used for analysis. Choose amongst 'Pearson', 'Spearman', and 'Biweight Mid' correlation. 'Biweight Mid' is default.
#'@return Graphics of gene correlations and a matrix p-values associated with analysis
#'@export
#'@example
#'tissue_pair_gene_expression(tissue_1, tissue_2, 0.2, 'Spearman')

tissue_pair_gene_expression <- function(tissue_1_matrix,tissue_2_matrix,variance_percentage,correlation_method) {
  ## Create a list of column names for downstream matching of donors
  tissue_1.samples = colnames(tissue_1_matrix)
  tissue_2.samples = colnames(tissue_2_matrix)

  ## Separate strings of column names to match donor tag identifiers
  tissue_1.subjects = sapply(tissue_1.samples,function(x){
    subj = strsplit(x,'[.]')[[1]]
    out = paste0(subj[1],'_',subj[2])
    return(out)
  })
  tissue_2.subjects = sapply(tissue_2.samples,function(x){
    subj = strsplit(x,'[.]')[[1]]
    out = paste0(subj[1],'_',subj[2])
    return(out)
  })
  combined_tissue.subjects = intersect(tissue_1.subjects,tissue_2.subjects)

  ## Index the matching donors so as to remove any excess information
  ind.tissue_1 = sapply(combined_tissue.subjects,function(x){
    out = which(tissue_1.subjects == x)
    return(out)
  })
  ind.tissue_2 = sapply(combined_tissue.subjects,function(x){
    out = which(tissue_2.subjects == x)
    return(out)
  })
  donor_matched_tissue_1.matrix = tissue_1_matrix[,ind.tissue_1]
  donor_matched_tissue_2.matrix = tissue_2_matrix[,ind.tissue_2]

  ## Perform similar data truncation with the gene information
  tissue_1.genes <- rownames(tissue_1_matrix)
  tissue_2.genes <- rownames(tissue_2_matrix)
  combined_tissue.genes <- intersect(tissue_1.genes,tissue_2.genes)
  ind.gene_tissue_1 = sapply(combined_tissue.genes,function(x){
    out = which(tissue_1.genes == x)
    return(out)
  })
  ind.gene_tissue_2 = sapply(combined_tissue.genes,function(x){
    out = which(tissue_2.genes == x)
    return(out)
  })
  final_tissue_1.matrix = donor_matched_tissue_1.matrix[ind.gene_tissue_1,]
  final_tissue_2.matrix = donor_matched_tissue_2.matrix[ind.gene_tissue_2,]
  #head(rownames(final_tissues_1.matrix)) == head(rownames(final_tissue_2.matrix))

  ## Solve for the variance within genes across donors
  library(iemisc)
  gene_expression_variance <- c()
  for (i in 1:size(as.matrix(final_tissue_1.matrix),1)) {
    gene_expression_variance[i] <- var(as.numeric(final_tissue_1.matrix[i,]),as.numeric(final_tissue_2.matrix[i,]))
  }

  ## Filter for the most variable genes (based on input variance_percentage)
  ind.variance <- which(gene_expression_variance >= quantile(gene_expression_variance, probs = 1 - variance_percentage))

  variable.tissue_1.genes <- final_tissue_1.matrix[ind.variance,]
  variable.tissue_2.genes <- final_tissue_2.matrix[ind.variance,]

  ## Test for correlations among genes between tissues, across donors

  if (correlation_method == 'pearson') {
    pvalues_pearson <- c()
    gene_pearson_correlation <- cor(t(as.matrix(variable.tissue_1.genes)),t(as.matrix(variable.tissue_2.genes)), method = "pearson")
    for (i in 1:size(as.matrix(variable.tissue_1.genes),1)){
      pvalues_pearson[i] <- cor.test(as.numeric(variable.tissue_1.genes[,i]),as.numeric(variable.tissue_2.genes[,i]), method = "pearson")$p.value
    }
    assign("gene_pearson_correlation",gene_pearson_correlation, envir = globalenv())
    assign("pvalues_pearson",pvalues_pearson, envir = globalenv())
  } else if (correlation_method == 'spearman') {
    gene_spearman_correlation <- cor(t(as.matrix(variable.tissue_1.genes)),t(as.matrix(variable.tissue_2.genes)), method = "spearman")
    pvalues_spearman <- cor.test(t(as.matrix(variable.tissue_1.genes)),t(as.matrix(variable.tissue_2.genes)), method = "spearman")$p.value
    assign("gene_spearman_correlation",gene_spearman_correlation, envir = globalenv())
    assign("pvalues_spearman",pvalues_spearman, envir = globalenv())
  }  else {
    library(WGCNA)
    gene_biweight_midcorrelation <- bicor(t(as.matrix(variable.tissue_1.genes[i,])),t(as.matrix(variable.tissue_2.genes[i,])))
    pvalues_biweight <- bicorAndPvalue(t(as.matrix(variable.tissue_1.genes[i,])),t(as.matrix(variable.tissue_2.genes[i,])))$p
    assign("gene_biweight_midcorrelation",gene_biweight_midcorrelation, envir = globalenv())
    assign("pvalues_biweight",pvalues_biweight, envir = globalenv())
  }
}
