#' @title Visualize Subtype Prediction Scores in a Sorted/ranked/stacked Barplot
#'
#' @description Visualize subtype prediction score within a set subtype in a stacked barplot.
#'
#' @details Take the output from `classify_samples` and return a stacked barplot plot representing 
#' the distribution of the, for that class, subtype prediction score. Set the subtype of 
#' desire with `this_subtype`. The subtype can be one of the subtypes included in the LundTax 
#' subtype classification nomenclature.
#'
#' @param these_predictions Required parameter, should be the output from 
#' [LundTaxR::classify_samples()].
#' @param this_subtype Required parameter. Should be one of the set subtype classes from the 
#' LundTax nomenclature.
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can further specify the dimensions of the returned plot with `plot_width` and 
#' `plot_height`.
#' @param plot_width This parameter controls the width in inches. Default is 4 (1200 pixels at 300 
#' PPI).
#' @param plot_height This parameter controls the height in inches. Default is 4 (1200 pixels at 300
#' PPI).
#' @param return_data Set to TRUE to return tidy data used by the plotting function. Default is 
#' FALSE.
#'
#' @return Nothing.
#' 
#' @import ggplot2 dplyr reshape2
#' @rawNamespace import(gridExtra, except = combine)
#'
#' @export
#'
#' @examples
#' #run predictor
#' sjodahl_classes = classify_samples(this_data = sjodahl_2017, 
#'                                    log_transform = FALSE, 
#'                                    adjust = TRUE, 
#'                                    impute = TRUE, 
#'                                    include_data = TRUE, 
#'                                    verbose = FALSE)
#' #plot Uro scores
#' plot_subscore_box(these_predictions = sjodahl_classes, 
#'                   this_subtype = "Uro")
#' 
#' #plot UroA scores
#' plot_subscore_box(these_predictions = sjodahl_classes, 
#'                   this_subtype = "UroA")
#' 
#' #return data
#' plot_data = plot_subscore_box(these_predictions = sjodahl_classes, 
#'                               this_subtype = "Uro", 
#'                               return_data = TRUE)
#'
#' #view data
#' head(plot_data)
#'
plot_subscore_box = function(these_predictions, 
                             this_subtype,
                             out_path = NULL,
                             out_format = "png",
                             plot_width = 4,
                             plot_height = 4,
                             return_data = FALSE){
  
  #subset scores
  these_scores = data.frame(these_predictions$subtype_scores)
  
  #get subtypes based on subtype class
  if(this_subtype %in% c("Uro", "GU", "BaSq", "Mes", "ScNE")){
    these_subtypes = data.frame(these_predictions$predictions_5classes)
    subtype_classes = c("Uro", "GU", "BaSq", "Mes", "ScNE")
  }else if(this_subtype %in% c("UroA", "UroB", "UroC")){
    these_subtypes = data.frame(these_predictions$predictions_7classes)
    subtype_classes = c("UroA", "UroB", "UroC")
  }else{
    stop("this_subtype must be one of the following; Uro, GU, BaSq, Mes, ScNE, UroA, UroB or UroC...")
  }
  
  #cnvert rowname to columns
  these_scores = tibble::rownames_to_column(these_scores, "sample_id")
  these_subtypes = tibble::rownames_to_column(these_subtypes, "sample_id")
  
  #rename column
  colnames(these_subtypes)[2] = "subtype"
  
  #subset consensus subtypes
  these_samples = dplyr::filter(these_subtypes, subtype == this_subtype)
  
  #replace NAs with zero
  these_scores[is.na(these_scores)] <- 0
  
  #melt data frame
  melted_scores = melt(data = these_scores,
                       id = "sample_id",
                       variable.name = "subtype",
                       measure.vars = c("Uro",
                                        "UroA",
                                        "UroB",
                                        "UroC",
                                        "GU",
                                        "BaSq",
                                        "Mes",
                                        "ScNE"))
  
  #convert to factor
  melted_scores$subtype = as.factor(melted_scores$subtype) 
  
  #subset on class
  class_melted = melted_scores %>%  
    dplyr::filter(subtype %in% subtype_classes)
  
  #get sample IDs for each subtype
  this_melted = dplyr::filter(class_melted, sample_id %in% these_samples$sample_id)
  
  #filter the data frame to get rows where subtype is the selected subtype
  subtype_df = this_melted[this_melted$subtype == this_subtype, ]
  
  #sort these rows by the value in descending order
  subtype_df = subtype_df[order(-subtype_df$value), ]
  
  #for each sample_id, sort the subtypes by their values in ascending order
  sorted_df = do.call(rbind, lapply(unique(subtype_df$sample_id), function(sample) {
    sample_df = this_melted[this_melted$sample_id == sample, ]
    sample_df = sample_df[order(sample_df$value), ]
    return(sample_df)
  }))
  
  #combine the sorted data frames
  sorted_df = sorted_df[order(match(sorted_df$sample_id, subtype_df$sample_id)), ]
  
  #reorder the sample_id factor levels based on the sorted subtype_df
  sorted_df$sample_id = factor(sorted_df$sample_id, levels = subtype_df$sample_id)
  
  #ensure the selected subtype's values are at the bottom
  sorted_df = sorted_df %>%
    group_by(sample_id) %>%
    arrange(desc(subtype == this_subtype), value, .by_group = TRUE) %>%
    ungroup()
  
  #reorder the subtype factor levels to ensure the selected subtype is plotted first
  sorted_df$subtype = factor(sorted_df$subtype, levels = c(this_subtype, setdiff(levels(sorted_df$subtype), this_subtype)))
  
  #get colour for the subtype to be plotted
  subtype_col = lund_colors$lund_colors[this_subtype]
  
  if(return_data){
    message("No plot generated, returning tidy data instead...")
    return(sorted_df)
  }

  #plot the sorted data frame with stacked bar plots
  my_plot = ggplot(sorted_df, aes(x = sample_id, y = value, fill = subtype)) +
    geom_col(width = 1) +
    scale_fill_manual(values = lund_colors$lund_colors) +
    theme_bw() +
    coord_cartesian(clip = "off") +
    ggtitle(label = paste0("Samples classified as: ", this_subtype)) +
    ylab("Prediction Score") +
    xlab("Samples") +
    labs(fill = "Subtype") +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks = seq(0, 1, by = 0.5)) +
    theme(legend.position = "right", 
          axis.line.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  if(!is.null(out_path)){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, this_subtype, "_subscore_box.pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, this_subtype, "_subscore_box.png"),
          width = plot_width,
          height = plot_height,
          units = "in",
          res = 300,
          pointsize = 10,
          bg = "white")
    }else{
      stop("Enter a valid output format (pdf or png)...")
    }
    print(my_plot)
    dev.off()
    message(paste0("Plot exported to ", out_path, this_subtype, "_subscore_box.", out_format))
  }else{
    return(my_plot) 
  }
}
