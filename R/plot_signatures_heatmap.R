#' @title Subtype-dependent Visualization of Signature Scores in a Heatmap
#'
#' @description Build a heatmap with signature scores retrieved with the [LundTaxR::classify_samples()] function.
#'
#' @details Construct and export (pdf or png) a highly customizable heatmap visualizing prediction
#' scores for each sample and class, predicted with [LundTaxR::classify_samples()].
#' This function depends on Complexheatmap. It is also possible to return a data frame with
#' prediction scores in a tidy format. To do so, set `return_scores = TRUE`. For a greater explanation
#' on how to use the function, see parameter descriptions and examples.
#'
#' @param these_predictions A list with a data frame object called scores. Returned with [LundTaxR::classify_samples()].
#' @param proportional_scores Set to TRUE to transform signature scores into proportions, default is FALSE.
#' @param out_path Optional, set path to export plot. If not provided, the heatmap is drawn in the active graphics device.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf". The user can further specify the dimensions of the returned plot with `plot_width` and `plot_height`.
#' @param return_scores Set to TRUE to return prediction scores in a tidy format. Default is FALSE.
#' @param to_xlsx Boolean parameter, set to TRUE to export score data frame in xlsx format. Default is FALSE. If set to TRUE, the spreadsheet will be saved to the same path as the heatmap.
#' @param plot_title Required parameter. Heatmap title, will also be pasted to the exported file(s) as well as a new column in the scores data frame under cohort.
#' @param hm_split Optional parameter for controlling how the data is split into different groups. If not provided, the function will split on what is specified within `subtype_annotation`.
#' @param subtype_annotation Can be one of the following; "5_class" (default) or "7_class" annotation.
#' @param hm_cluster Boolean parameter, set to TRUE to cluster the rows (default is FALSE).
#' @param plot_anno_legend Expects a vector with TRUE/FALSE (6 in total, one per annotation track: Predicted Subtype, Proliferation Score, Mol. grade (WHO1999), Mol. grade (WHO2022), Progression Score, Progression Risk), this decides what legends will be on the final heatmap. Default is c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE). If set to NULL, only the subtype legend is shown.
#' @param plot_hm_legend Boolean parameter. Set to TRUE to show heatmap legend. Default is FALSE.
#' @param plot_width This parameter controls the width in inches. Default is 14(4200 pixels at 300 PPI).
#' @param plot_height This parameter controls the height in inches. Default is 6(1800 pixels at 300 PPI)
#' @param plot_font_size Optional parameter to control the size of the font in the generated heatmap. Note, the title of the plot will always be twice that of the set value here (default = 12).
#' @param verbose Set to FALSE to suppress messages. Default is TRUE.
#'
#' @return Data frame with prediction score for each sample and class, if return_scores = TRUE. Otherwise, nothing.
#'
#' @import ComplexHeatmap ggplot2 circlize grid
#' @importFrom stats median quantile
#' @importFrom grDevices dev.off pdf png
#' @importFrom openxlsx write.xlsx
#' @importFrom tibble rownames_to_column
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
#' #plot 5 class
#' plot_signatures_heatmap(these_predictions = sjodahl_classes, 
#'                         plot_title = "Signature Scores (Sjodahl 2017)",
#'                         return_scores = FALSE, 
#'                         plot_anno_legend = FALSE,
#'                         subtype_annotation = "5_class")
#'                         
#' #plot 7 class
#' plot_signatures_heatmap(these_predictions = sjodahl_classes, 
#'                         plot_title = "Signature Scores (Sjodahl 2017)",
#'                         return_scores = FALSE, 
#'                         plot_anno_legend = FALSE,
#'                         subtype_annotation = "7_class")
#'
plot_signatures_heatmap = function(these_predictions = NULL,
                                   proportional_scores = FALSE,
                                   out_path = NULL,
                                   out_format = "png",
                                   return_scores = FALSE,
                                   to_xlsx = FALSE,
                                   plot_title = NULL,
                                   hm_split = NULL,
                                   subtype_annotation = "5_class",
                                   hm_cluster = FALSE,
                                   plot_anno_legend = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
                                   plot_hm_legend = FALSE,
                                   plot_width = 14,
                                   plot_height = 6,
                                   plot_font_size = 12,
                                   verbose = TRUE){
  
  #deal with nonsensical parameter combinations
  if(to_xlsx){
    return_scores = TRUE
  }
  
  if(is.null(plot_title)){
    if(!return_scores){
      stop("Please specify a title for your plot with `plot_title`...")
    }else{
      stop("Please provide a title for your data with `plot_title`...")
    }
  } 
  
  #set the default legend properties
  if(is.null(plot_anno_legend)){
    plot_anno_legend = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  }else{
    plot_anno_legend = plot_anno_legend
  }
  
  #check the type of incoming data
  if(is.null(these_predictions) | !is.list(these_predictions)){
    stop("Input should be the output from classify_samples()...")
  }
  
  #set up output device if path is provided
  save_to_file <- !is.null(out_path)
  
  if(save_to_file){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, plot_title, "_heatmap_scores.pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, plot_title, "_heatmap_scores.png"),
          width = plot_width,
          height = plot_height,
          units = "in",
          res = 300,
          pointsize = 10,
          bg = "white")
    }else{
      stop("Enter a valid output format (pdf or png)...")
    }
    }
  
  #create the heatmap regardless of output path (for display or saving)
  #set sample order
  sample_order = order(these_predictions$scores$proliferation_score)
  
  #set split (top annotation track of heatmap)
  if(!is.null(hm_split)){
    split = hm_split
  }else{
    if(subtype_annotation %in% c("5_class", "7_class")){
      if(subtype_annotation == "5_class"){
        split = factor(these_predictions$predictions_5classes,
                       levels = c("Uro", "GU", "BaSq", "Mes", "ScNE")) 
      }else if(subtype_annotation == "7_class"){
        split = factor(these_predictions$predictions_7classes,
                       levels = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")) 
      }
    }else{
      stop("subtype_annotation must be one of the following; 5_class or 7_class...")
    }
  }
  
  if(proportional_scores){
    message("transforming signature scores into proportions...")
    
    #running helper
    prop_scores = int_get_prop_scores(these_predictions = these_predictions)
    
    #update the these_predictions$scores object
    these_predictions$scores <- prop_scores
    
    #immune + stromal names
    immune_stromal_names <- c("immune141_up", "stromal141_up", "immune_stroma_ratio")
    
    #get immune names
    immune_names <- c("b_cells_proportion", "t_cells_proportion",
                      "t_cells_cd8_proportion", "nk_cells_proportion",
                      "cytotoxicity_score_proportion", "neutrophils_proportion",
                      "monocytic_lineage_proportion", "macrophages_proportion",
                      "m2_macrophage_proportion", "myeloid_dendritic_cells_proportion")
    
    #get stromal names
    stromal_names <- c("endothelial_cells_proportion",
                       "fibroblasts_proportion", "smooth_muscle_proportion")
    
    immune_labels = c("B Cells Prop.", "T Cells Prop.", "CD8+ T Cells Prop.", "NK Cells Prop.", 
                      "Cytotoxicity Score Prop.", "Neutrophils Prop.", "Monocytic Lineage Prop.", 
                      "Macrophages Prop.", "M2 Macrophages Prop.", "Myeloid DCs Prop.")
    
    stromal_labels = c("Endothelial Cells Prop.", "Fibroblasts Prop.", "Smooth Muscle Prop.")

  }else{
    #get immune names
    immune_names <- c("immune141_up", "b_cells", "t_cells",
                      "t_cells_cd8", "nk_cells",
                      "cytotoxicity_score", "neutrophils",
                      "monocytic_lineage", "macrophages",
                      "m2_macrophage", "myeloid_dendritic_cells")
    
    #get stromal names
    stromal_names <- c("stromal141_up", "endothelial_cells",
                       "fibroblasts", "smooth_muscle")  
    
    immune_labels = c("Immune141_UP", "B Cells", "T Cells", "CD8+ T Cells", "NK Cells", 
                      "Cytotoxicity Score", "Neutrophils", "Monocytic Lineage", 
                      "Macrophages", "M2 Macrophages", "Myeloid DCs")

    stromal_labels = c("Stromal 141_UP", "Endothelial Cells", "Fibroblasts", "Smooth Muscle")
  }
  
  #set colours
  #proliferation
  col_fun_proliferation =
    circlize::colorRamp2(c(quantile(these_predictions$scores$proliferation_score, 0.05),
                           median(these_predictions$scores$proliferation_score),
                           quantile(these_predictions$scores$proliferation_score, 0.95)),
                         c("#2166AC","white", "#B2182B"))
  
  #progression
  col_fun_progression =
    circlize::colorRamp2(c(quantile(these_predictions$scores$progression_score, 0.05),
                           median(these_predictions$scores$progression_score),
                           quantile(these_predictions$scores$progression_score, 0.90)),
                         c("#2166AC","white", "#B2182B"))
  
  #create colour object
  colour_obj = list(`Predicted Subtype` = lund_colors$lund_colors,
                    `Proliferation Score` = col_fun_proliferation,
                    `Mol. grade (WHO1999)` = c("G1_2" = "grey", "G3" = "black"),
                    `Mol. grade (WHO2022)` = c("HG" = "black", "LG" = "grey"),
                    `Progression Score` = col_fun_progression,
                    `Progression Risk` = c("HR" = "black", "LR" = "grey"))
  
  # Add custom color functions for immune/stromal signatures
  col_fun_immune = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

  #plotting
  #build heatmap annotation (top)
  hm <- HeatmapAnnotation(`Predicted Subtype` = split,
                          `Proliferation Score` = these_predictions$scores$proliferation_score,
                          `Mol. grade (WHO1999)` = these_predictions$scores$molecular_grade_who_1999,
                          `Mol. grade (WHO2022)` = these_predictions$scores$molecular_grade_who_2022,
                          `Progression Score` = these_predictions$scores$progression_score,
                          `Progression Risk` = these_predictions$scores$progression_risk,
                          annotation_name_side = "left",
                          show_legend = plot_anno_legend,
                          annotation_name_gp = gpar(fontsize = plot_font_size),
                          border = TRUE,
                          col = colour_obj,
                          column_title = NULL,
                          row_title = NULL)
  
  if(proportional_scores){
    #imune + stromal
    hm_immune_stromal <- Heatmap(t(scale(these_predictions$scores[,immune_stromal_names, drop = FALSE])),
                                 top_annotation = hm,
                                 column_order = sample_order,
                                 height = unit(5*ncol(these_predictions$scores[,immune_stromal_names,drop = FALSE]), "mm"),
                                 width = unit(0.5*nrow(these_predictions$scores[,immune_stromal_names,drop = FALSE]), "mm"),
                                 name = "combined",
                                 border = TRUE,
                                 column_split = split,
                                 col = col_fun_immune,
                                 cluster_rows = hm_cluster,
                                 show_heatmap_legend = plot_hm_legend,
                                 row_title = NULL,
                                 row_names_side = "left",
                                 row_names_gp = gpar(fontsize = plot_font_size),
                                 show_column_names = FALSE,
                                 column_title = NULL,
                                 row_labels = c("Immune 141_UP", "Stromal 141_UP", "Immune Stroma Ratio"))
    
    #immune scores heatmap
    hm_immune_scores <- Heatmap(t(scale(these_predictions$scores[,immune_names, drop = FALSE])),
                                column_order = sample_order,
                                height = unit(5*ncol(these_predictions$scores[,immune_names,drop = FALSE]), "mm"),
                                width = unit(0.5*nrow(these_predictions$scores[,immune_names,drop = FALSE]), "mm"),
                                name = "Immune Scores",
                                border = TRUE,
                                column_split = split,
                                col = col_fun_immune,
                                cluster_rows = hm_cluster,
                                show_heatmap_legend = plot_hm_legend,
                                row_names_side = "left",
                                row_names_gp = gpar(fontsize = plot_font_size),
                                show_column_names = FALSE,
                                column_title = NULL,
                                row_title = NULL,
                                row_labels = immune_labels)
  }else{
    #immune scores heatmap
    hm_immune_scores <- Heatmap(t(scale(these_predictions$scores[,immune_names, drop = FALSE])),
                                top_annotation = hm,
                                column_order = sample_order,
                                height = unit(5*ncol(these_predictions$scores[,immune_names,drop = FALSE]), "mm"),
                                width = unit(0.5*nrow(these_predictions$scores[,immune_names,drop = FALSE]), "mm"),
                                name = "Immune Scores",
                                border = TRUE,
                                column_split = split,
                                col = col_fun_immune,
                                cluster_rows = hm_cluster,
                                show_heatmap_legend = plot_hm_legend,
                                row_names_side = "left",
                                row_names_gp = gpar(fontsize = plot_font_size),
                                show_column_names = FALSE,
                                column_title = NULL,
                                row_title = NULL,
                                row_labels = immune_labels)
  }
  
  #stroma scores heatmap
  hm_stroma_scores <- Heatmap(t(scale(these_predictions$scores[,stromal_names, drop = FALSE])),
                              column_order = sample_order,
                              height = unit(5*ncol(these_predictions$scores[,stromal_names, drop = FALSE]), "mm"),
                              width = unit(0.5*nrow(these_predictions$scores[,stromal_names,drop = FALSE]), "mm"),
                              name = "Stroma Scores",
                              border = TRUE,
                              column_split = split,
                              col = col_fun_immune,
                              cluster_rows = hm_cluster,
                              show_heatmap_legend = plot_hm_legend,
                              row_names_side = "left",
                              row_names_gp = gpar(fontsize = plot_font_size),
                              show_column_names = FALSE,
                              column_title = NULL,
                              row_title = NULL,
                              row_labels = stromal_labels)
  
  if(proportional_scores){
    #combine heatmaps
    hm_combined = draw(hm_immune_stromal %v% hm_immune_scores %v% hm_stroma_scores,
                       column_title = plot_title,     heatmap_legend_side = "right",
                       annotation_legend_side = "right",
                       column_title_gp = gpar(fontsize = (plot_font_size*2))); 
  }else{
    #combine heatmaps
    hm_combined = draw(hm_immune_scores %v% hm_stroma_scores,
                       column_title = plot_title,     heatmap_legend_side = "right",
                       annotation_legend_side = "right",
                       column_title_gp = gpar(fontsize = (plot_font_size*2)));
  }
    
    hm_sample_order <- column_order(hm_combined@ht_list$`Stroma Scores`)
    
    #close graphics device only if saving to file
    if(save_to_file){
      dev.off()
      if(verbose){
        message(paste0("Heatmap exported to ", out_path, plot_title, "_heatmap_scores.", out_format))
      }
    } else {
      if(verbose){
        message("No output path provided, displaying heatmap in R...")
      }
    }
  
  if(return_scores == TRUE){
    
    #subset scores data frame from list
    my_scores = these_predictions$scores
    my_scores = as.data.frame(my_scores)
    
    #add new column with cohort information
    my_scores$cohort = plot_title
    
    #convert correct columns to factors
    factor_cols <- c('molecular_grade_who_1999' ,'molecular_grade_who_2022',
                     'progression_risk', 'cohort')
    
    my_scores[,factor_cols] <- lapply(my_scores[,factor_cols] , factor)
    
    if(to_xlsx){
      my_scores <- tibble::rownames_to_column(my_scores, "SampleID")
      write.xlsx(my_scores, paste0(out_path, plot_title, "_scores.xlsx"))
      if(verbose){
        message(paste0("Prediction scores exported as xlsx to ", out_path, plot_title, "_scores.xlsx"))
      }
    }
    
    if(verbose){
      message("Prediction Scores for each sample and class sucessfully generated!")
    }
    
    #return data frame
    return(my_scores)
  }
}
