gene_names <- function(tissue_matrix_1,tissue_matrix_2){
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
  
  ## Index the matching donors so as to remove non-matching information 
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
  
  assign("final_tissue_1.matrix",final_tissue_1.matrix, envir = globalenv())
  assign("final_tissue_2.matrix",final_tissue_2.matrix, envir = globalenv())
  
  
}