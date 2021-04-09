
#' Ligand Filter Function
#'
#' Takes in matrix of genes and filters the data according to ligand presence.
#' ligand_filter() returns a matrix of genes with matching ligands if the filter.for parameter is TRUE.
#' Otherwise it will return a matrix with the genes removed from the matrix.
#'@param matrix1 A file containing gene-donor data of the first tissue. Recommend using 'final_tissue_1.matrix' acquired from the use of gene_names()
#'@param filter.for boolean entry that determines whether to retian or remove the genes within the inputted matrix. If 'TRUE'
#'function will return the matrix with genes that contain ligands. If 'FALSE' it will remove the genes with ligands. 'Defaulted to FALSE'
#'@return Two matricies of the most variable genes according to variance percentage.
#'export

ligand_filter <- function(matrix1, selection = NULL, filter.for = F) {

  # load library and data sets
  library(GSEAplot)
  data("Kadoki_ligands.db")
  data("Kadoki_receptors.db")

  # assign gene names into objects
  ligand_symbols <- get_genesymbols(Kadoki_ligands.db)$Kadoki_ligand
  receptor_symbols <- get_genesymbols(Kadoki_receptors.db)$Kadoki_receptor

  # construct vector containing both ligand and receptors symbols, excluding overlap
  overlap_id <- match(intersect(ligand_symbols, receptor_symbols), ligand_symbols)
  both_symbols <- c(ligand_symbols[-overlap_id], receptor_symbols)

  # this outputs the indecies where the selected genes exist in the passed matrix
  env <- new.env()

  env[["ligands"]] <- match(intersect(rownames(matrix1), ligand_symbols), rownames(matrix1))

  env[["receptors"]] <- match(intersect(rownames(matrix1), receptor_symbols), rownames(matrix1))

  env[["both"]] <- match(intersect(rownames(matrix1), both_symbols), rownames(matrix1))

  # filtering to keep only rows that contain ligands (NOT DEFAULT)
  if (filter.for) {

    return(matrix1[env[[selection]],])

  } else {

    # return matrix without rows that contain the selected genes
    return(matrix1[-env[[selection]],])

  }
}
