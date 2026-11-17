#' @title Rank Genes Within Samples
#'
#' @description Internal function called by [LundTaxR::int_ratio_score()]. Not meant for out of 
#' package use.
#'
#' @details Returns the same values as `apply(this_data, 2, rank) / nrow(this_data)` (average ranks
#' for ties), but only for the requested genes. Each sample is sorted once and the ranks of the 
#' requested genes are looked up with [findInterval()], which avoids ranking every gene in every 
#' sample. If the data contains missing values, the function falls back to [rank()] so that 
#' missing values are handled exactly as before.
#'
#' @param this_data Required parameter. Matrix (or data frame) with expression values, genes as 
#' rows and samples as columns.
#' @param genes Required parameter. Character vector with the genes (row names) to return ranks for.
#'
#' @return A matrix with genes as rows and samples as columns, holding each gene's rank within the 
#' sample divided by the total number of genes.
#'
#' @examples
#' \dontrun{
#' # No examples provided
#' }
#'
int_rank_genes = function(this_data = NULL,
                          genes = NULL){

  this_data = as.matrix(this_data)
  n_genes = nrow(this_data)

  #fall back to rank() if there are missing values, to keep the original NA handling
  if(anyNA(this_data)){
    rank_data = apply(this_data, 2, rank)
    rownames(rank_data) = rownames(this_data)
    return(rank_data[genes, , drop = FALSE] / n_genes)
  }

  gene_data = this_data[genes, , drop = FALSE]
  rank_data = gene_data

  #average rank = (number of values below + number of values below or equal + 1) / 2
  for(j in seq_len(ncol(this_data))){
    sorted_values = sort.int(this_data[, j], method = "radix")
    n_below = findInterval(gene_data[, j], sorted_values, left.open = TRUE)
    n_below_or_equal = findInterval(gene_data[, j], sorted_values)
    rank_data[, j] = (n_below + n_below_or_equal + 1) / 2
  }

  return(rank_data / n_genes)
}
