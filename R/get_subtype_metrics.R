#' @title Get Summarizing Subtype Metrics
#'
#' @description This function tables a metadata column based on subtype classification.
#' 
#' @details The function performs the following steps:
#' \itemize{
#'   \item Ensures the sample IDs are present in both metadata and predictions.
#'   \item Checks if the desired metadata column is valid.
#'   \item Filters metadata and predictions to include only common samples.
#'   \item Combines metadata and predictions into a single data frame.
#'   \item Counts the number of samples and occurrences of specified level for each subtype.
#' }
#'
#' @param this_metadata A data frame containing metadata with a `sample_id` column.
#' @param this_metadata_variable A string specifying the column name in the metadata to be 
#' summarized.
#' @param these_predictions Output from [LundTaxR::classify_samples()]. The prediction vectors 
#' (`predictions_5classes` / `predictions_7classes`) must have sample IDs as names.
#' @param factor_level The level of the factor/value in `this_metadata_variable` to count 
#' (optional). If not provided, will return counts for all levels.
#' @param subtype_class The classification system, default is 5 class.
#'
#' @return A data frame with the number of samples and counts for each subtype.
#' If `factor_level` is specified, returns counts for that specific level.
#' If `factor_level` is NULL, returns counts for all levels.
#' 
#' @import dplyr
#' @importFrom tidyr pivot_wider
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
#' #check progression events in each subtype                    
#' get_subtype_metrics(these_predictions = sjodahl_classes, 
#'                     this_metadata = sjodahl_2017_meta, 
#'                     this_metadata_variable = "surv_os_event", 
#'                     factor_level = 1)
#'
#' #check number of males in each subtype
#' get_subtype_metrics(these_predictions = sjodahl_classes, 
#'                     this_metadata = sjodahl_2017_meta, 
#'                     this_metadata_variable = "gender", 
#'                     factor_level = "Male")
#'
get_subtype_metrics = function(this_metadata = NULL,
                               this_metadata_variable = NULL,
                               these_predictions = NULL,
                               factor_level = NULL,
                               subtype_class = "5_class") {
  
  #check if sample_id is present in metadata
  if (!"sample_id" %in% colnames(this_metadata)) {
    stop("The metadata must contain a 'sample_id' column.")
  }
  
  #check subtype class
  if (!subtype_class %in% c("5_class", "7_class")) {
    stop("subtype_class must be one of the following; 5_class or 7_class...")
  }
  
  #check if sample_id is present in predictions
  prediction_vector = if (subtype_class == "5_class") these_predictions$predictions_5classes else these_predictions$predictions_7classes
  if (is.null(names(prediction_vector))) {
    stop("The predictions must have sample IDs as names.")
  }
  
  #check if the desired metadata column is valid
  if (is.null(this_metadata_variable)) {
    stop("Please specify a metadata column with `this_metadata_variable`...")
  }
  
  if (!this_metadata_variable %in% colnames(this_metadata)) {
    stop("The specified metadata variable is not a valid column in the metadata.")
  }
  
  if(subtype_class == "5_class"){
    #ensure the sample IDs are the same in both metadata and predictions
    common_samples <- intersect(names(these_predictions$predictions_5classes), this_metadata$sample_id)
    
    #filter metadata and predictions to include only common samples
    filtered_metadata <- this_metadata %>% filter(sample_id %in% common_samples)
    filtered_predictions <- these_predictions$predictions_5classes[common_samples]
  }else if(subtype_class == "7_class"){
    #ensure the sample IDs are the same in both metadata and predictions
    common_samples <- intersect(names(these_predictions$predictions_7classes), this_metadata$sample_id)
    
    #filter metadata and predictions to include only common samples
    filtered_metadata <- this_metadata %>% filter(sample_id %in% common_samples)
    filtered_predictions <- these_predictions$predictions_7classes[common_samples]
  }
  
  #combine metadata and predictions into a single data frame
  combined_data <- filtered_metadata %>%
    mutate(prediction = filtered_predictions[match(sample_id, names(filtered_predictions))])
  
  #count based on whether factor_level is specified
  if (is.null(factor_level)) {
    #return counts for all levels
    result <- combined_data %>%
      group_by(prediction, .data[[this_metadata_variable]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = all_of(this_metadata_variable), 
                  values_from = count, 
                  values_fill = 0) %>%
      mutate(total_samples = rowSums(select(., -prediction)))
  } else {
    #return counts for specific level
    result <- combined_data %>%
      group_by(prediction) %>%
      summarise(
        total_samples = n(),
        !!paste0(this_metadata_variable, "_", factor_level) := sum(.data[[this_metadata_variable]] == factor_level, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  return(result)
}
