#' @title Subset a List With Prediction Results
#'
#' @description This function subsets the data frames and named vectors in a list based on a set of 
#' sample IDs.
#' 
#' @details This function is designed to subset various components of a predictions list based on a 
#' set of sample IDs. If `these_sample_ids` is not provided, the function will attempt to extract 
#' sample IDs from `these_samples_metadata`. The function handles different structures within the 
#' list, including data frames and named vectors, ensuring that only the specified samples are 
#' retained. Note, this function expects that `classify_samples` has been run with `include_data = TRUE`, 
#' if not the `data` component will be missing.
#'
#' @param these_predictions Output from `classify_samples` function run with `include_data = TRUE`.
#' @param these_sample_ids A vector of sample IDs to subset the data. If NULL, sample IDs will be 
#' extracted from `these_samples_metadata`.
#' @param these_samples_metadata A data frame containing sample metadata. If provided and 
#' `these_sample_ids` is NULL, sample IDs will be extracted from this data frame.
#' @param samples_rownames A logical value indicating whether the sample IDs are stored in the row 
#' names of `these_samples_metadata`. Default is TRUE.
#'
#' @return A list with the same structure as `these_predictions`, but subset to the specified 
#' sample IDs.
#' 
#' @import dplyr tibble
#' 
#' @export
#'
#' @examples
#' sjodahl_classes = classify_samples(this_data = sjodahl_2017, 
#'                                    log_transform = FALSE, 
#'                                    adjust = TRUE, 
#'                                    impute = TRUE, 
#'                                    include_data = TRUE, 
#'                                    verbose = FALSE)
#'                                    
#' #get some sample subsets
#' my_samples = head(sjodahl_2017_meta$sample_id)
#' my_meta = head(sjodahl_2017_meta)
#'
#' #use sample IDs
#' data_subset <- get_prediction_subset(these_predictions = sjodahl_classes, 
#'                                      these_sample_ids = my_samples)
#'
#' #use a metadata subset 
#' data_subset <- get_prediction_subset(these_predictions = sjodahl_classes, 
#'                                      these_samples_metadata = my_meta,
#'                                      samples_rownames = FALSE)
#' #view data
#' head(data_subset$data)
#' 
#' #view subtype scores
#' head(data_subset$subtype_scores)
#' 
#' #view prediction classes
#' head(data_subset$predictions_5classes)
#'
get_prediction_subset = function(these_predictions = NULL,
                                 these_sample_ids = NULL,
                                 these_samples_metadata = NULL,
                                 samples_rownames = FALSE){

  #check if data is inlcuded in the prediction object
  if(is.null(these_predictions$data)){
    stop("The 'these_predictions' object does not contain a 'data' component. \n  Please ensure that the 'classify_samples' function was run with 'include_data = TRUE'.")
  }

  if(!is.null(these_samples_metadata) && is.null(these_sample_ids)){
    message("Metadata provided, the function will subset to the sample IDs in this object...")
    if(!samples_rownames){
      message("CAUTION: sample_rownames = FALSE, the funciton expects a column in the metadata called `sample_id`")
    }
    if(samples_rownames){
      these_sample_ids = these_samples_metadata %>% 
        rownames_to_column("sample_id") %>% 
        pull(sample_id)
    }else{
      these_sample_ids = these_samples_metadata %>% 
        pull(sample_id)
    }
  }else{
    message("Sample IDs detected...")
  }
  
  #subset the 'data' data frame
  these_predictions$data <- these_predictions$data %>%
    select(all_of(these_sample_ids))
  
  #subset the 'subtype_scores' data frame
  these_predictions$subtype_scores <- these_predictions$subtype_scores[these_sample_ids, , drop = FALSE]
  
  #subset the 'predictions_7classes' named vector
  these_predictions$predictions_7classes <- these_predictions$predictions_7classes[these_sample_ids]
  
  #subset the 'predictions_5classes' named vector
  these_predictions$predictions_5classes <- these_predictions$predictions_5classes[these_sample_ids]
  
  #subset the 'scores' data frame
  these_predictions$scores <- these_predictions$scores %>%
    filter(rownames(these_predictions$scores) %in% these_sample_ids)
  
  return(these_predictions)
}
