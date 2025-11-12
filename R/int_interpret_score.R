#' @title Interpret SwitchBox Score
#' 
#' @description Interpret SwitchBox Scores as Binary Classifications
#' 
#' @details Internal function called by [LundTaxR::int_calc_signatures()]. Not meant for out of 
#' package use. Converts continuous prediction scores (0-1) to categorical grade labels using a 
#' threshold cutoff.  
#'
#' @param scores Numeric vector of scores (0-1)
#' @param threshold Cutoff for binary classification (default 0.5)
#' @param labels Character vector of length 2: c("low_class", "high_class")
#'
#' @return Character vector of classifications
#' 
#' @examples
#' \dontrun{
#' # No examples provided
#' }
#' 
int_interpret_score <- function(scores, 
                                threshold = 0.5, 
                                labels = NULL){
  
  ifelse(scores >= threshold, labels[2], labels[1])
}

