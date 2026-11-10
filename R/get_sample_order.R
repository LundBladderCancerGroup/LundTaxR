#' @title Get Sample Order based on Proliferation Score.
#'
#' @description Return a sample order based on the late/early ratio (proliferation score).
#'
#' @details Modified version of the heatmap function [LundTaxR::plot_classification_heatmap()], 
#' this function only returns the order of samples to be used for downstream plotting functions. 
#' The order is dictated by the early late ratio. and the sample order can readily be sued 
#' elsewhere, if needed.
#' 
#' @param expr_data Required parameter, should be the output from [LundTaxR::classify_samples()]. 
#' Note, this function required the `classify_samples` to be run with `include_data = TRUE` 
#' argument.
#' @param return_this Set to "late_early" to return late/early ratio. Set to "sample_order" to 
#' return the sample order based on late/early ratio (default).
#' 
#' @return Returns the sample order.
#' 
#' @import dplyr
#' 
#' @export 
#'
#' @examples
#' #run classifier
#' sjodahl_classes = classify_samples(this_data = sjodahl_2017, 
#'                                    log_transform = FALSE, 
#'                                    adjust = TRUE, 
#'                                    impute = TRUE, 
#'                                    include_data = TRUE, 
#'                                    verbose = FALSE)
#'        
#' #get sample order                            
#' get_sample_order(expr_data = sjodahl_classes$data)
#'
get_sample_order = function(expr_data = NULL,
                            return_this = "sample_order"){
  
  #check incoming data and parameter combinations
  this_data = as.matrix(expr_data)
  
  #gene signatures for plotting
  genes_to_plot = list(Early_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "EarlyCellCycle"), 1]),
                       Late_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "LateCellCycle"), 1]),
                       Late_Early = NULL)
  
  #get genes in provided data for downstream filtering steps
  these_genes = row.names(this_data)
  genes_early = intersect(genes_to_plot$Early_CC, these_genes)
  genes_late = intersect(genes_to_plot$Late_CC, these_genes)
  
  #combine genes
  genes_cc = na.omit(c(genes_early, genes_late))
  
  #row split for the heatmap
  row_split = c(rep("Early", length(genes_early)),
                rep("Late", length(genes_late)))
  #row title
  row_title_cc = c("Early Cell Cycle", "Late Cell Cycle")
  
  #late and Early scores
  late_score = apply(this_data[intersect(rownames(this_data),genes_to_plot$Late_CC),], 2, median)
  early_score = apply(this_data[intersect(rownames(this_data),genes_to_plot$Early_CC),], 2, median)
  
  #ratio
  late_early_ratio = late_score - early_score
  
  #add genes to genes_to_plot object
  genes_to_plot$Late_Early = late_early_ratio
  
  #order samples by late_early cell cycle
  sample_order = order(late_early_ratio)
  
  if(return_this == "late_early"){
    message("Returning late/early ratio...")
    return(as.data.frame(late_early_ratio))
  }else if(return_this == "sample_order"){
    message("Returning sample order...")
    return(sample_order)
  }else{
    stop("return_this must be one of the following; sample_order or late_early...")
  }
}
