#'
#' Bicor Analysis Function //Working Title
#'
#'Takes in two tissue data sets and attempts to match the genes within both tissues.
#'This function also attempts to label the genes with the proper gene names.
#'Returns two matrices for each tissue with the matching genes.
#'
#'@param tissue_1_matrix A file containing gene-donor data of the first tissue
#'@param tissue_2_matrix A file containing gene-donor data of the second tissue
#'@param percentage_variance A decimal entry representing the variance percentage for gene analysis.
#'@return A series of plots in a pdf associated with the bicor midcorrelation analysis conducted on the inputted tissue matrices
#'@export


tissue_pair_gene_expression <- function(tissue_1_matrix,tissue_2_matrix,percentage_variance) {

  ## Test for correlations among genes between tissues, in a single donor

  library(WGCNA)
  gene_biweight_midcorrelation <- bicor(t(as.matrix(tissue_1_matrix)),t(as.matrix(tissue_2_matrix)))
  pvalues_biweight <- bicorAndPvalue(t(as.matrix(tissue_1_matrix)),t(as.matrix(tissue_2_matrix)))$p
  assign("gene_biweight_midcorrelation",gene_biweight_midcorrelation, envir = globalenv())
  assign("pvalues_biweight",pvalues_biweight, envir = globalenv())
  #cc_and_pvalues = list(Bicor = gene_biweight_midcorrelation,Pvalues = pvalues_biweight)

  # Write a pdf file to contain data visualizations
  pdf("Analysis_Visualization.pdf")

  # Create a histogram of correlation coefficients
  hist(gene_biweight_midcorrelation,main = 'Histogram of Gene Correlation Coefficients',ylab = 'Frequency',xlab = 'Correlation Magnitude')

  # Create a matrix of inverse log base 10 p-values
  logpvalues <- -log(pvalues_biweight)

  # load library and data sets
  library(GSEAplot)
  data("Kadoki_ligands.db")
  data("Kadoki_receptors.db")

  # assign gene names into objects
  ligand_symbols <- get_genesymbols(Kadoki_ligands.db)$Kadoki_ligand

  # find the indeces for genes encoding for secreted ligands in origin tissue
  ind.encoding_tissue_1 <- match(intersect(rownames(tissue_1_matrix), ligand_symbols), rownames(tissue_1_matrix))

  # filtering to keep only rows for the origin tissue that encode for secreted ligands
  pvalues_ligand_filtered <- logpvalues[ind.encoding_tissue_1,]

  # Plot the distribution of inverse log pvalues for origin tissue genes encoding secreted proteins
  hist(pvalues_ligand_filtered, main = 'Histogram of Inverse Log Pvalues for Origin Tissue Ligand Secreting Genes',xlab = 'Inverse Log Pvalue Magnitude',ylab = 'Frequency', breaks = 100)
  abline(v=quantile(pvalues_ligand_filtered, probs = 1 - percentage_variance), col = 'red')

  # Create an array of scores for each ligand encoding gene in the origin tissue
  scores <- rowSums(pvalues_ligand_filtered)

  # Normalize the scores for the origin tissue ligand encoding genes, add gene name identifiers, arrange in descending order
  library(dplyr)
  Ssec = data.frame(Gene_symbol = names(scores), score = scores) %>%
    mutate(Ssec = score / length(colnames(tissue_2_matrix))) %>%
    arrange(desc(Ssec))
  Ssec = Ssec[,c(1,3)]

  # Provide a table with ligand encoding gene significance scores
  write.table(Ssec,file = 'Ssec_scores_all_encoding_genes', row.names=F, col.names=T, sep='\t', quote=F)

  # Melt the biweight midcorrelation matrix
  library(reshape2)
  bicor.data = melt(gene_biweight_midcorrelation)
  pvalue.data = melt(logpvalues)
  colnames(bicor.data) = c('Gene_symbol_1', 'Gene_symbol_2', 'bicor')
  colnames(pvalue.data) = c('Gene_symbol_1', 'Gene_symbol_2', 'pvalue')

  # Isolate data for the origin gene with the highest Ssec score
  ##MAJOR EDIT: added toString()
  gene_interest= filter(bicor.data, Gene_symbol_1 == toString( Ssec[1,1]) ) %>%
    rename(Gene_symbol = Gene_symbol_2)
  gene_interest = gene_interest[,c(2,3)]

  # Positively correlated pathways - "gene enhancement"
  upreg = gene_interest %>% arrange(desc(bicor)) %>% head(5)
  write.table(upreg, file = 'Positively_correlated_pathways_top_rank_gene', col.names=F, sep='\t', quote=F)

  # Negatively correlated pathways - "gene suppression"
  downreg = gene_interest %>% arrange(bicor) %>% head(5)
  write.table(downreg, file = 'Negatively_correlated_pathways_top_rank_gene', col.names=F, sep='\t', quote=F)

  # All pathways engaged by gene of interest
  totalreg = gene_interest %>% arrange(desc(abs(bicor))) %>% head(10)
  write.table(totalreg, file="Absolute_value_pathways_top_rank_gene", col.names=F, sep='\t', quote=F)

  # Find the gene pair indeces for the ligand encoding genes with the highest inverse log pvalues
  library(iemisc)
  cutoff.ind <- which(pvalues_ligand_filtered>= quantile(pvalues_ligand_filtered, probs = 1 - percentage_variance), arr.ind = T)

  # Create variables containing gene names from our filtered pvalue list and our matched matrices for plotting
  rows.ref = rownames(pvalues_ligand_filtered)
  columns.ref = colnames(pvalues_ligand_filtered)
  rows = rownames(matched_pairs$tissue_1)
  columns = rownames(matched_pairs$tissue_2)

  # Create plots equal to the number of gene pairs with high magnitude pvalues
  for (i in 1:(size(cutoff.ind,1))){
    # Find the gene names corresponding to the indeces in our filtered pvalue matrix
    gene.name1 <- rows.ref[cutoff.ind[i,1]]
    gene.name2 <- columns.ref[cutoff.ind[i,2]]
    # Find the indeces of high magnitude pvalue gene pairs
    gene1.ind = which(rows == gene.name1)
    gene2.ind = which(columns == gene.name2)
    # Add a regression model
    I = lm(formula = as.numeric(matched_pairs$tissue_1[gene1.ind,]) ~ as.numeric(matched_pairs$tissue_2[gene2.ind,]))
    # Plot the expressions of the two genes with high degree of correlation
    plot(as.numeric(matched_pairs$tissue_1[gene1.ind,]),as.numeric(matched_pairs$tissue_2[gene2.ind,]),main = 'Gene Pair Expression',xlab = gene.name1,ylab = gene.name2)
    # Add the linear regression line to the plot
    abline(as.numeric(I$coefficients[1]),as.numeric(I$coefficients[2]), col = 'red')
  }
  dev.off()
}
