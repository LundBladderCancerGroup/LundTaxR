#' @title Calculate Signature Ratio Scores
#'
#' @description Calculate proliferation and progression scores.
#'
#' @details Internal function called by [LundTaxR::int_calc_signatures()].
#' Not meant for out of package use. Takes a data frame of matrix with expression values and
#' calculates scores based on gene expression.
#'
#' @param this_data Required parameter. Data frame or matrix with expression values.
#' @param variable Required parameter. Input should be one of the following; proliferation, or
#' progression.
#' @param gene_id Specify the type of gene identifier used in `this_data`. Accepted values are;
#' hgnc_symbol (default) or ensembl_gene_id.
#' @param verbose A logical value indicating whether processing messages will be printed or not.
#' Default is TRUE.
#'
#' @return A list with two objects. 1, A data frame with scores for the selected variable. 2, A data
#' frame indicating what genes from the incoming data are missing, based on the expected genes for
#' signature calculations.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # No examples provided
#' }
#'
int_ratio_score = function(this_data = NULL,
                           variable = NULL,
                           gene_id = "hgnc_symbol",
                           verbose = TRUE){

  #valid variables
  valid_vars = c("immune", "score141up", "proliferation", "progression")

  #check valid variables
  if(!variable %in% valid_vars){
    stop(paste0(
      "return_this must be one of the following: ",
      paste(valid_vars, collapse = ", ")
    ))
  }

  if(variable == "proliferation"){
    #get all genes with signature LateCellCycle
    up_genes = filter(signatures$proliferation, signature == "LateCellCycle") %>%
      select(signature, !!as.symbol(gene_id))

    #get all genes with signature EarlyCellCycle
    down_genes = filter(signatures$proliferation, signature == "EarlyCellCycle") %>%
      select(signature, !!as.symbol(gene_id))
  }else if(variable == "progression"){
    #get all genes with signature LateCellCycle
    up_genes = filter(signatures$progression, direction_in_prog == "Up") %>%
      select(direction_in_prog, !!as.symbol(gene_id))

    #get all genes with signature EarlyCellCycle
    down_genes = filter(signatures$progression, direction_in_prog == "Down") %>%
      select(direction_in_prog, !!as.symbol(gene_id))
  }

  #intersect with incoming data
  int_up_genes = this_data %>%
    dplyr::filter(rownames(this_data) %in% unique(up_genes[,gene_id]))

  int_down_genes = this_data %>%
    dplyr::filter(rownames(this_data) %in% unique(down_genes[,gene_id]))

  #get genes that are not in the incoming data
  diff_genes_up = setdiff(up_genes[,gene_id], rownames(int_up_genes))
  diff_genes_down = setdiff(down_genes[,gene_id], rownames(int_down_genes))

  #create data frame with missing gene information
  missing_genes = data.frame(genes = as.character(),
                             signature = as.character(),
                             process = as.character())

  #append missing genes information to the data frame
  missing_genes = add_row(missing_genes,
                          genes = diff_genes_up,
                          signature = variable,
                          process = "up")

  missing_genes = add_row(missing_genes,
                          genes = diff_genes_down,
                          signature = variable,
                          process = "down")

  #ratio of ranks
  rank_data = apply(this_data, 2, rank)
  rank_data = rank_data/nrow(rank_data)

  #get vector with relevant genes
  up_vector = rownames(int_up_genes)
  down_vector = rownames(int_down_genes)

  median_up_genes = apply(rank_data[up_vector,,drop=F],2,median)
  median_down_genes = apply(rank_data[down_vector,,drop=F],2,median)

  median_up_down = median_up_genes/median_down_genes

  r_score = data.frame(Score = median_up_down,
                       row.names = colnames(this_data))

  return(list(sig_score = r_score, na_genes = missing_genes))
}
