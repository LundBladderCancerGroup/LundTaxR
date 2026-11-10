#' @title Subtype Forest Plot
#'
#' @description Create a subtype aware forest plot, with a set clinical end-point.
#'
#' @details This function creates a forest plot for multiple subtypes based on Cox proportional 
#' hazards models. It also includes a table with the number of samples and progression events for 
#' each subtype.
#'
#' @param these_predictions Required parameter if `this_metadata` is not provided. Should be output 
#' from [LundTaxR::classify_samples()].
#' @param these_samples_metadata Required parameter if `these_predictions` is not provided or to 
#' provide survival data. Metadata associated with the prediction output.
#' @param this_subtype A character vector specifying the subtypes to be included in the plot. 
#' Default is c("UroA", "UroB", "UroC", "Uro", "GU", "BaSq", "Mes", "ScNE").
#' @param subtype_class Can be one of the following; 5_class or 7_class. Default is 5_class.
#' @param surv_event A string specifying the column name for the survival event.
#' @param surv_time A string specifying the column name for the survival time.
#' @param sample_id_col Optional parameter. Allows the user to manually specify the name of a column
#' with sample ID.
#' @param row_to_col Optional parameter, set to TRUE to convert row names in metadata to a new 
#' column called sample_id. Default is FALSE.
#' @param plot_title A string specifying the title of the plot. Default is "Subtype Forest Plot".
#' @param plot_subtitle A string specifying the subtitle of the plot.
#' @param plot_caption A string specifying the caption of the plot.
#' @param plot_width A numeric value specifying the width of the plot. Default is 7.
#' @param plot_height A numeric value specifying the height of the plot. Default is 7.
#' @param out_format A string specifying the output format of the plot ("png" or "pdf"). Default is 
#' "png".
#' @param significant_p Numeric parameter for flagging significant p values. Default is 0.05.
#' @param sig_color Color for annotating significant signatures. Default is red.
#' @param out_path A string specifying the output path for saving the plot.
#' @param file_name A string specifying the file name for the saved plot. Default is "subtype".
#' @param return_data Set to TRUE to return plot data, default is FALSE.
#' 
#' @return A combined plot object if `out_path` is NULL. Otherwise, the plot is saved to the 
#' specified path.
#' 
#' @import dplyr survival ggplot2 patchwork
#' @importFrom cowplot plot_grid
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
#' #UroA
#' plot_subtype_forest(these_samples_metadata = sjodahl_2017_meta,
#'                     these_predictions = sjodahl_classes,
#'                     this_subtype = "UroA",
#'                     subtype_class = "7_class",
#'                     surv_event = "surv_css_event",
#'                     surv_time = "surv_css_time")
#'
#' #all 5 class subtypes
#' plot_subtype_forest(these_samples_metadata = sjodahl_2017_meta,
#'                     these_predictions = sjodahl_classes,
#'                     this_subtype = NULL,
#'                     subtype_class = "5_class",
#'                     surv_event = "surv_css_event",
#'                     surv_time = "surv_css_time")
#'
#' #return data
#' forest_data = plot_subtype_forest(these_samples_metadata = sjodahl_2017_meta,
#'                                   these_predictions = sjodahl_classes,
#'                                   subtype_class = "5_class",
#'                                   this_subtype = c("Uro", "GU", "BaSq", "Mes", "ScNE"),
#'                                   return_data = TRUE,
#'                                   surv_event = "surv_css_event",
#'                                   surv_time = "surv_css_time")
#'                                   
#' #view data
#' head(forest_data)
#' 
plot_subtype_forest = function(these_predictions = NULL,
                               these_samples_metadata = NULL,
                               this_subtype = NULL,
                               subtype_class = "5_class",
                               surv_event = NULL,
                               sig_color = "red",
                               significant_p = 0.05,
                               surv_time = NULL,
                               sample_id_col = NULL, 
                               row_to_col = FALSE,
                               plot_title = NULL,
                               plot_subtitle = NULL,
                               plot_caption = NULL,
                               plot_width = 7,
                               plot_height = 7,
                               out_format = "png",
                               out_path = NULL,
                               file_name = NULL,
                               return_data = FALSE){
  
  if(is.null(this_subtype)){
    if(subtype_class == "5_class"){
      this_subtype = c("Uro", "GU", "BaSq", "Mes", "ScNE")
    }else if(subtype_class == "7_class"){
      this_subtype = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")
    }
  }

  #check required parameters
  if (is.null(these_predictions) && is.null(these_samples_metadata)) {
    stop("Either 'these_predictions' or 'these_samples_metadata' must be provided.")
  }
  
  if (is.null(surv_event) || is.null(surv_time)) {
    stop("Both 'surv_event' and 'surv_time' must be specified.")
  }
  
  #process predictions and metadata similar to other functions
  if (!is.null(these_predictions)) {
    #use int_prediction_wrangler or similar logic from other functions
    if (is.null(these_samples_metadata)) {
      stop("'these_samples_metadata' is required when using 'these_predictions'.")
    }
    
    #handle sample ID column
    if (row_to_col && !"sample_id" %in% colnames(these_samples_metadata)) {
      these_samples_metadata <- these_samples_metadata %>%
        tibble::rownames_to_column("sample_id")
    }
    
    if (!is.null(sample_id_col)) {
      these_samples_metadata <- these_samples_metadata %>%
        rename(sample_id = !!sample_id_col)
    }
    
    #ensure sample_id column exists
    if (!"sample_id" %in% colnames(these_samples_metadata)) {
      stop("No sample_id column found. Use 'row_to_col = TRUE' or specify 'sample_id_col'.")
    }
    
    #extract predictions based on subtype_class
    if (subtype_class == "5_class") {
      predictions_vector <- these_predictions$predictions_5classes
      available_subtypes <- names(table(predictions_vector))
    } else if (subtype_class == "7_class") {
      predictions_vector <- these_predictions$predictions_7classes
      available_subtypes <- names(table(predictions_vector))
    } else {
      stop("subtype_class must be either '5_class' or '7_class'")
    }
    
    #filter to common samples
    common_samples <- intersect(names(predictions_vector), these_samples_metadata$sample_id)
    
    if (length(common_samples) == 0) {
      stop("No common samples found between predictions and metadata.")
    }
    
    #create metadata with predictions
    meta_fixed <- these_samples_metadata %>%
      filter(sample_id %in% common_samples) %>%
      mutate(
        predictions = predictions_vector[match(sample_id, names(predictions_vector))]
      )
    
  } else {
    #use provided metadata directly (backward compatibility)
    meta_fixed <- these_samples_metadata
    
    #determine which prediction columns to use
    if (subtype_class == "5_class" && "Predictions_5classes" %in% colnames(meta_fixed)) {
      meta_fixed <- meta_fixed %>% rename(predictions = Predictions_5classes)
    } else if (subtype_class == "7_class" && "Predictions_7classes" %in% colnames(meta_fixed)) {
      meta_fixed <- meta_fixed %>% rename(predictions = Predictions_7classes)
    } else {
      stop("No suitable prediction column found in metadata for the specified subtype_class.")
    }
    
    available_subtypes <- unique(meta_fixed$predictions)
  }
  
  #filter subtypes to only those available in the data
  this_subtype <- intersect(this_subtype, available_subtypes)
  
  if (length(this_subtype) == 0) {
    stop("None of the specified subtypes are available in the data.")
  }
  
  #define helpers (same as original)
  fit_cox_models_sub = function(my_data, 
                                time_column, 
                                event_column, 
                                pred_columns){
    
    if(!all(c(time_column, event_column, pred_columns) %in% colnames(my_data))){
      stop("One or more specified columns do not exist in the data frame.")
    }
    
    my_data[[time_column]] <- as.numeric(my_data[[time_column]])
    my_data[[event_column]] <- as.numeric(my_data[[event_column]])
    
    my_data = my_data %>%
      mutate(across(all_of(pred_columns), ~ factor(., levels = c(0, 1))))
    
    cox_models <- lapply(pred_columns, function(col) {
      formula <- as.formula(paste("Surv(", time_column, ", ", event_column, ") ~", col))
      coxph(formula, data = my_data)
    })
    
    names(cox_models) <- pred_columns
    return(cox_models)
  }
  
  extract_surv_stats = function(cox){
    my_stats = data.frame()
    
    for(name in names(cox)){
      stats = data.frame(score = name)
      
      cox_model = cox[[name]]
      stats$p_value = summary(cox_model)$coefficients[5]
      stats$hazard_ratio = exp(coef(cox_model))
      conf_int = exp(confint(cox_model))
      stats$hazard_conf_2.5 <- conf_int[1]
      stats$hazard_conf_97.5 <- conf_int[2]
      
      my_stats = rbind(my_stats, stats)
    }
    return(my_stats)
  }
  
  #create binary matrix for subtypes
  this_subtype_with_bin <- paste0("bin_", this_subtype)
  
  #create binary columns for each subtype
  for (subtype in this_subtype) {
    bin_col_name <- paste0("bin_", subtype)
    meta_fixed[[bin_col_name]] <- factor(
      ifelse(meta_fixed$predictions == subtype, 1, 0), 
      levels = c(0, 1)
    )
  }
  
  #fit Cox models
  cox_models = fit_cox_models_sub(my_data = meta_fixed, 
                                  time_column = surv_time, 
                                  event_column = surv_event, 
                                  pred_columns = this_subtype_with_bin)
  
  #extract statistics
  stat_results = extract_surv_stats(cox = cox_models)
  
  #process results (same as original)
  stat_results$significant = ifelse(stat_results$p_value < significant_p, "significant", "not significant")
  
  names(stat_results)[3] = "ratio"
  names(stat_results)[4] = "conf_2.5"
  names(stat_results)[5] = "conf_97.5"
  
  stat_results = stat_results %>% 
    rename(subtype = score) %>%
    mutate(subtype = factor(sub("bin_", "", subtype)))
  
  stat_results$subtype <- factor(stat_results$subtype, levels = rev(this_subtype))
  
  if(return_data){
    message("No plot generated, returning data instead...")
    return(stat_results)
  }
  
  #create plot (same as original)
  my_plot = ggplot(stat_results, aes(x = ratio, y = subtype)) +
    geom_point(aes(color = significant), size = 3) +
    geom_errorbar(aes(xmin = conf_2.5, xmax = conf_97.5, color = significant), width = 0.2, orientation = "y") +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "red") +
    labs(title = if(is.null(plot_title)) "Subtype Forest Plot" else plot_title, subtitle = plot_subtitle, caption = plot_caption, x = "Hazard Ratio", y = "") +
    scale_color_manual(values = c("significant" = sig_color, "not significant" = "black")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  #create table
  sample_counts <- sapply(this_subtype_with_bin, function(subtype) sum(meta_fixed[[subtype]] == 1))
  event_counts <- sapply(this_subtype_with_bin, function(subtype) {
    sum(meta_fixed[[subtype]] == 1 & meta_fixed[[surv_event]] == 1, na.rm = TRUE)
  })
  
  table_data <- data.frame(
    Subtype = this_subtype,
    Samples_Events = paste0(event_counts, "/", sample_counts)
  )
  
  table_data$Subtype <- factor(table_data$Subtype, levels = levels(stat_results$subtype))
  
  numbers_plot <- ggplot(table_data, aes(y = Subtype)) +
    geom_text(aes(x = 1, label = Samples_Events), size = 3, hjust = 0, fontface = "italic", family = "sans") +
    labs(y = NULL, x = NULL) +
    theme_void() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    xlim(1, 2)
  
  combined_plot <- cowplot::plot_grid(
    my_plot,
    numbers_plot,
    ncol = 2,
    rel_widths = c(4, 1),
    align = 'h',
    axis = 'tb'
  )
  
  #export or return plot
  if(!is.null(out_path)){
    if(is.null(file_name)){
      file_name = "subtype"
    }
    if(out_format == "pdf"){
      pdf(paste0(out_path, file_name, "_", "hazard_ratio", "_forest.pdf"),
          width = plot_width,
          height = plot_height)
    }else if(out_format == "png"){
      png(paste0(out_path, file_name, "_", "hazard_ratio", "_forest.png"),
          width = plot_width,
          height = plot_height,
          units = "in",
          res = 300,
          pointsize = 10,
          bg = "white")
    }else{
      stop("Enter a valid output format (pdf or png)...")
    }
    print(combined_plot)
    dev.off()
    message(paste0("Plot exported to ", out_path, file_name, "_", "hazard_ratio", "_forest.", out_format))
  }else{
    return(combined_plot) 
  }
}