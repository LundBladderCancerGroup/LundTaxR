#' @title Create Custom Annotation Tracks from Metadata
#'
#' @description Generate HeatmapAnnotation objects from metadata that can be added to classification heatmaps.
#'
#' @details This function takes a metadata data frame and creates annotation tracks that can be 
#' integrated into the main classification heatmap. It automatically handles different data types 
#' and provides sensible color schemes.
#' 
#' @param metadata A data frame with samples as rows and annotation variables as columns. 
#' Row names should match sample IDs in the expression data.
#' @param annotation_vars Character vector specifying which columns from metadata to include. 
#' If NULL (default), all columns will be used.
#' @param sample_order Optional, vector specifying sample order. Should match the order used 
#' in the main heatmap.
#' @param custom_colors Optional named list of color schemes for specific annotations. 
#' Names should match annotation variable names.
#' @param show_legends Logical vector or single value indicating whether to show legends 
#' for each annotation. Default is FALSE.
#' @param annotation_height Height of annotation tracks in mm. Default is 4.
#' @param font_size Font size for annotation labels. Default is 8.
#' 
#' @return HeatmapAnnotation object that can be passed to plot_classification_heatmap
#' 
#' @import ComplexHeatmap circlize
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats na.omit
#' 
#' @export
#'
#' @examples
#' custom_annot <- get_custom_annotations(
#'   metadata = sjodahl_2017_meta,
#'   annotation_vars = c("gender", "turb_stage"),
#'   show_legends = TRUE)
#'
get_custom_annotations <- function(metadata = NULL,
                                   annotation_vars = NULL,
                                   sample_order = NULL,
                                   custom_colors = NULL,
                                   show_legends = FALSE,
                                   annotation_height = 4,
                                   font_size = 8){
  
  #input validation
  if(is.null(metadata)) {
    stop("metadata is required")
  }
  
  if(!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  
  if(is.null(annotation_vars)) {
    annotation_vars <- colnames(metadata)
  }
  
  #filter metadata to selected variables
  metadata_subset <- metadata[, annotation_vars, drop = FALSE]
  
  #apply sample order if provided
  if(!is.null(sample_order)) {
    #match sample order with metadata rownames
    common_samples <- intersect(sample_order, rownames(metadata_subset))
    if(length(common_samples) == 0) {
      warning("No samples match between sample_order and metadata rownames")
    } else {
      metadata_subset <- metadata_subset[common_samples, , drop = FALSE]
    }
  }
  
  #function to generate colors for different data types
  generate_colors <- function(values, var_name) {
    
    #check if custom colors provided
    if(!is.null(custom_colors) && var_name %in% names(custom_colors)) {
      return(custom_colors[[var_name]])
    }
    
    #remove NA values for color generation
    clean_values <- na.omit(values)
    
    if(is.numeric(clean_values)) {
      #continuous variable - use color ramp
      col_fun <- circlize::colorRamp2(
        c(quantile(clean_values, 0.05, na.rm = TRUE),
          median(clean_values, na.rm = TRUE),
          quantile(clean_values, 0.95, na.rm = TRUE)),
        c("#2166AC", "white", "#B2182B")
      )
      return(col_fun)
      
    } else {
      #categorical variable - use discrete colors
      unique_vals <- unique(as.character(clean_values))
      n_colors <- length(unique_vals)
      
      if(n_colors <= 3) {
        colors <- c("grey", "black", "#B2182B")[1:n_colors]
      } else if(n_colors <= 8) {
        colors <- RColorBrewer::brewer.pal(max(3, n_colors), "Set2")[1:n_colors]
      } else if(n_colors <= 12) {
        colors <- RColorBrewer::brewer.pal(n_colors, "Set3")
      } else {
        #for many categories, use rainbow
        colors <- rainbow(n_colors)
      }
      
      names(colors) <- unique_vals
      return(colors)
    }
  }
  
  #create annotation list and color list
  annotation_list <- list()
  color_list <- list()
  
  for(var in annotation_vars) {
    if(var %in% colnames(metadata_subset)) {
      values <- metadata_subset[[var]]
      annotation_list[[var]] <- values
      color_list[[var]] <- generate_colors(values, var)
    }
  }
  
  #handle show_legends parameter
  if(length(show_legends) == 1) {
    show_legends <- rep(show_legends, length(annotation_vars))
  }
  names(show_legends) <- annotation_vars
  
  #create HeatmapAnnotation
    annotation <- HeatmapAnnotation(
      df = annotation_list,
      col = color_list,
      annotation_name_side = "left",
      annotation_name_rot = 0,
      show_legend = show_legends,
      annotation_name_gp = gpar(fontsize = font_size),
      simple_anno_size = unit(annotation_height, "mm"),
      border = TRUE,
      gap = unit(1, "mm")
    )
  
  return(annotation)
}
