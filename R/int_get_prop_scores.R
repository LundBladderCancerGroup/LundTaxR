#' @title Get Proportional Signature Scores
#' 
#' @description Internal function called by `plot_signatures_heatmap` when `proportional_scores = TRUE`.
#' Not meant for out-of-package usage.
#'
#' @details This function calculates the immune/stroma ratio (`immune141_up` minus `stromal141_up`) 
#' and appends it to the scores. The proportional immune and stromal scores used for plotting are 
#' already part of the scores object (the `*_proportion` columns).
#'
#' @param these_predictions A list containing a data frame named `scores` with prediction scores. 
#' The data frame should have row names as sample IDs and columns corresponding to different cell 
#' types and scores.
#'
#' @return A data frame with the original scores and the calculated immune/stroma ratio. The row 
#' names are the sample IDs.
#'
#' @import dplyr tibble
#' 
#' @examples
#' \dontrun{
#' # No examples provided
#' }
#' 
int_get_prop_scores = function(these_predictions = NULL){
  
  #get scores
  updated_scores = as.data.frame(these_predictions$scores)
  
  #calculate stroma/immune ratio
  updated_scores$immune_stroma_ratio = updated_scores$immune141_up - updated_scores$stromal141_up
  
  return(updated_scores)
}
