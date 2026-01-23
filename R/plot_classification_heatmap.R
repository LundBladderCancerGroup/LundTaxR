#' @title Visualize Output from the LundTaxR Classifier in a Heatmap
#'
#' @description Plot heatmap with classification results.
#'
#' @details This function plots a heatmap including genes and signatures of interest, with 
#' prediction results and scores on top. Optionally includes signature scores heatmap.
#' 
#' @param these_predictions Required parameter, should be the output from 
#' [LundTaxR::classify_samples()].
#' @param this_data Expression data used for predictions. 
#' Required if the output from [LundTaxR::classify_samples()] is run with 
#' include_data = FALSE (default).
#' @param gene_id Specify the type of gene identifier used in `this_data`. Accepted values are; 
#' hgnc_symbol (default) or ensembl_gene_id.
#' @param subtype_annotation Can be one of the following; "5_class" (default) or "7_class" 
#' annotation.
#' @param this_sample_order Optional, set sample order. Default, samples are split by subtype, and 
#' order within each subtype. By default, samples are order by late/early cell cycle ratio (low to 
#' high).
#' @param norm Boolean parameter. Set to TRUE (default) to normalize the data into Z-scaled values.
#' @param plot_scores Boolean parameter. Set to TRUE to plot prediction scores for each class.
#' Default is FALSE.
#' @param plot_signature_scores Boolean parameter. Set to TRUE to add signature scores heatmap 
#' after the classification heatmap. Default is FALSE.
#' @param custom_annotation Optional HeatmapAnnotation object to add below the existing 
#' subtype prediction annotations and score bars. Should be created with `get_custom_annotations()`.
#' Default is NULL.
#' @param include_sections Named list specifying which heatmap sections to include. Available sections:
#' early_late_cc, luminal_tfs, luminal_genes, fgfr3, circuit, tp63, basq, keratinization, erbb, 
#' adhesion, myc, neuronal. Each should be TRUE/FALSE. Default includes all sections (TRUE).
#' @param show_ann_legend Boolean parameter, set to TRUE to show annotation legend (Lund classes). 
#' Default is FALSE.
#' @param show_hm_legend Boolean parameter, set to TRUE to show heatmap legend, default is FALSE.
#' @param ann_height Plotting parameter, optional. Annotation height in cm. Default = 8.
#' @param plot_title Plotting parameter. The title for the generated heatmap. Default is "My Plot".
#' @param plot_width This parameter controls the width in inches. Default is 14 (4200 pixels at 300 
#' PPI).
#' @param plot_height This parameter controls the height in pixels. Default is 10 (3000 pixels at 
#' 300 PPI)
#' @param plot_font_size Optional parameter to control the size of the font in the generated 
#' heatmap. Note, the title of the plot will always be twice that of the set value here 
#' (default = 10).
#' @param plot_font_row_size Optional parameter to control the size of the font in the generated 
#' heatmap. Note, the title of the plot will always be twice that of the set value here 
#' (default = 8).
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can control the dimensions with `plot_width` and `plot_height`.
#' 
#' @return Draws heatmap and silently returns the sample order.
#' 
#' @import ComplexHeatmap circlize dplyr grid
#' @importFrom stats na.omit quantile
#' @importFrom grDevices dev.off pdf png
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
#' #plot
#' plot_classification_heatmap(these_predictions = sjodahl_classes, 
#'                             subtype_annotation = "5_class",
#'                             plot_scores = FALSE,
#'                             plot_title = "Classification Results (Sjodahl 2017)", 
#'                             show_ann_legend = TRUE,
#'                             ann_height = 0.5)
#'
plot_classification_heatmap = function(these_predictions = NULL,
                                       this_data = NULL,
                                       gene_id = "hgnc_symbols",
                                       subtype_annotation = "5_class",
                                       this_sample_order = NULL,
                                       norm = TRUE,
                                       plot_scores = FALSE,
                                       plot_signature_scores = FALSE,
                                       custom_annotation = NULL,
                                       include_sections = list(
                                         early_late_cc = TRUE,
                                         luminal_tfs = TRUE,
                                         luminal_genes = TRUE,
                                         fgfr3 = TRUE,
                                         circuit = TRUE,
                                         tp63 = TRUE,
                                         basq = TRUE,
                                         keratinization = TRUE,
                                         erbb = TRUE,
                                         adhesion = TRUE,
                                         myc = TRUE,
                                         neuronal = TRUE
                                       ),
                                       show_ann_legend = FALSE,
                                       show_hm_legend = FALSE,
                                       ann_height = 0.5,
                                       plot_title = "My Plot",
                                       plot_width = 14,
                                       plot_height = 11,
                                       plot_font_size = 10,
                                       plot_font_row_size = 7,
                                       out_path = NULL,
                                       out_format = "png"){
  
  #check incoming data and parameter combinations
  if(is.null(these_predictions) | !is.list(these_predictions)){
    stop("Input should be the result of applying landtax_predict_sub...")
  }else if(is.list(these_predictions) & "data" %in% names(these_predictions)){
    this_data = these_predictions$data
    score_matrix = these_predictions$subtype_scores
    pred_labels5 = these_predictions$predictions_5classes
    pred_labels7 = these_predictions$predictions_7classes
  }else if(is.list(these_predictions) & !"data" %in% names(these_predictions)){
    if(is.null(this_data)){
      stop("Data is missing. Include it in results object by running predict_lundtax with include_data = TRUE or provide it in the this_data argument...")
    }else if(!(class(this_data)[1] %in%  c("matrix","data.frame"))){
      stop("Data should be in matrix or data.frame format")
    }else if(class(this_data)[1] %in%  c("matrix","data.frame")){
      this_data = this_data[,names(these_predictions$predictions_7classes)]
      score_matrix = these_predictions$subtype_scores
      pred_labels5 = these_predictions$predictions_5classes
      pred_labels7 = these_predictions$predictions_7classes
    }
  }
  
  #scale and set this_data as matrix
  if(norm){
    this_data = as.matrix(this_data)
    this_data = scale(t(this_data))
    this_data = t(this_data)
  }else{
    this_data = as.matrix(this_data)
  }
  
  #convert ensembl gene id to hgnc symbols
  if(gene_id == "ensembl_gene_id"){
    
    #convert the rownames to first column
    mutated_data = dplyr::as_tibble(this_data, 
                                    rownames = "ensembl_gene_id")
    
    #left join with gene list to get ensembl IDs
    mutated_data = dplyr::left_join(gene_list, 
                                    mutated_data, 
                                    by = "ensembl_gene_id")
    
    #remove duplicated rows
    mutated_data = dplyr::filter(mutated_data, 
                                 duplicated(hgnc_symbol) == FALSE)
    
    #convert the first column back to rownames
    row.names(mutated_data) = unique(mutated_data$hgnc_symbol)
    mutated_data[1:2] = NULL
    
    #convert back to expected name
    this_data = mutated_data
  }
  
  #gene signatures for plotting
  genes_to_plot = list(Early_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "EarlyCellCycle"), 1]),
                       Late_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "LateCellCycle"), 1]),
                       Late_Early = NULL,
                       UroDiff = c("PPARG", "FOXA1", "GATA3", "ELF3"),
                       UPKs = c("UPK1A", "UPK1B", "UPK2", "UPK3A", "KRT20"),
                       Circuit = c("FGFR3", "CCND1", "E2F3", "RB1", "CDKN2A"),
                       Circuit_score = NULL,
                       FGFR3 = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "FGFR3"), 1]),
                       BaSq = c("KRT5", "KRT14", "FOXA1", "GATA3"),
                       BaSq_ratio = NULL,
                       Keratinization = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "Keratinization_QTC"), 1]),
                       Adhesion = c("EPCAM", "CDH1", "CDH3"),
                       MYC = c("MYCL", "MYCN", "MYC"),
                       ERBB = c("EGFR", "ERBB2", "ERBB3"),
                       ERBB_score = NULL,
                       ScNE = c("CHGA", "SYP", "ENO2"),
                       Immune141_UP = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "Immune141_UP"), 1]),
                       Stromal141_UP = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "Stromal141_UP"), 1]),
                       Immune141_UP_score = NULL,
                       Stromal141_UP_score = NULL)
  
  #create wrapper function for subtype score plots
  bar_anno = function(plot_scores = plot_scores,
                      subtype = NULL){
    
    if(plot_scores){
      this_bar = anno_barplot(as.numeric(score_matrix[,subtype]), 
                              axis = FALSE, 
                              ylim = c(0, 1), 
                              gp = gpar(fill = lund_colors$lund_colors[subtype], 
                                        border = NA, 
                                        col = NA), 
                              bar_width = 1, 
                              height = unit(5, "mm"))
    }else{
      this_bar = NULL
    }
    return(this_bar)
  }
  
  ##### DRAW HEATMAPS #####
  #################################### SUBTYPE SCORES ##############################################
  #get subtype colors - reorder for desired legend sequence
  #create ordered color mapping for 5-class system
  ordered_colors_5class = lund_colors$lund_colors[c("Uro", "GU", "BaSq", "Mes", "ScNE")]
  col = list(`Predicted Subtype` = ordered_colors_5class)

  #predictions
  pred_lab = pred_labels5

  split = factor(pred_lab, levels = c("Uro", "GU", "BaSq", "Mes", "ScNE"))
  bar1 = bar_anno(plot_scores = plot_scores, subtype = "Uro")
  bar2 = bar_anno(plot_scores = plot_scores, subtype = "GU")
  bar3 = bar_anno(plot_scores = plot_scores, subtype = "BaSq")
  bar4 = bar_anno(plot_scores = plot_scores, subtype = "Mes")
  bar5 = bar_anno(plot_scores = plot_scores, subtype = "ScNE")
  bar6 = bar_anno(plot_scores = FALSE, subtype = "UroA")
  bar7 = bar_anno(plot_scores = FALSE, subtype = "UroB")
  bar8 = bar_anno(plot_scores = FALSE, subtype = "UroC")
  
  #add Uro subtypes
  if(subtype_annotation == "7_class"){
    
    #predictions
    pred_lab = pred_labels7
    
    #overwrite predictions
    split = factor(pred_lab ,levels = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE"))
    
    #create ordered color mapping for 7-class system
    ordered_colors_7class = lund_colors$lund_colors[c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")]
    col = list(`Predicted Subtype` = ordered_colors_7class)
    
    #add score plots for Uro subtypes
    bar6 = bar_anno(plot_scores = plot_scores, subtype = "UroA")
    bar7 = bar_anno(plot_scores = plot_scores, subtype = "UroB")
    bar8 = bar_anno(plot_scores = plot_scores, subtype = "UroC")
    
  }else if(!subtype_annotation %in% c("5_class", "7_class")){
    stop("Please provide a valid subtype annotation (5 classes or 7 classes)...")
  }
  
  #Heatmap annotations, prediction scores
  hm_a_predictions = HeatmapAnnotation(`Predicted Subtype` = pred_lab,
                                       annotation_name_side = "left", 
                                       annotation_name_rot = 0,
                                       col = col,
                                       Uro = bar1,
                                       GU = bar2, 
                                       `Ba/Sq` = bar3, 
                                       Mes = bar4,
                                       ScNE = bar5,
                                       UroA = bar6,
                                       UroB = bar7,
                                       UroC = bar8,
                                       na_col = "gray83", 
                                       gap = unit(2, "mm"),
                                       simple_anno_size = unit(10, "mm"),
                                       simple_anno_size_adjust = TRUE,
                                       show_legend = show_ann_legend,
                                       border = TRUE,
                                       height = unit(ann_height, "cm"), 
                                       annotation_name_gp = gpar(fontsize = plot_font_size),
                                       annotation_legend_param = list(`Predicted Subtype` = list(
                                         at = if(subtype_annotation == "7_class") c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE") else c("Uro", "GU", "BaSq", "Mes", "ScNE"),
                                         labels = if(subtype_annotation == "7_class") c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE") else c("Uro", "GU", "BaSq", "Mes", "ScNE")
                                       )))
  
  ################################### EARLY / LATE CELL CYCLE ######################################
  #get genes in provided data for downstream filtering steps
  these_genes = row.names(this_data)
  genes_early = intersect(genes_to_plot$Early_CC, these_genes)
  genes_late = intersect(genes_to_plot$Late_CC, these_genes)
  
  #combine genes
  genes_cc = na.omit(c(genes_early, genes_late))
  
  #check if genes from both early and late are present
  if(length(genes_early) != 0 & !all(is.na(genes_late))){
    
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
    
    #create color palette for late/early
    col_fun_cc = circlize::colorRamp2(c(quantile(late_early_ratio, 0.05),
                                        median(late_early_ratio),
                                        quantile(late_early_ratio, 0.95)),
                                      c("#2166AC","white", "#B2182B"))
    
    if(is.null(this_sample_order)){
      #order samples by late_early cell cycle
      sample_order = order(late_early_ratio) 
    }else{
      message("Custom sample order will be used...")
      sample_order = this_sample_order
    }
    
    #heatmap annotation
    col = list(Predictions = lund_colors$lund_colors,
               `Late/Early Ratio` = col_fun_cc)
    
    #draw annotations track, late/early
    hm_a_lateearly = HeatmapAnnotation(`Late/Early Ratio` = genes_to_plot$Late_Early,
                                       annotation_name_side = "left",
                                       gap = unit(1, "mm"),
                                       simple_anno_size = unit(4, "mm"),
                                       simple_anno_size_adjust = TRUE,
                                       col = col,
                                       show_legend = FALSE,
                                       border = TRUE,
                                       annotation_name_gp = gpar(fontsize = plot_font_row_size))
    
  }else if(length(genes_early) == 0 & !all(is.na(genes_late))){
    message("All genes from the early cell cycle signature are missing.
            \nSamples will not be ordered by cell cycle...")
    sample_order = NULL
    ha1b = NULL
    row_split = NULL
    row_title_cc = "Late cell cycle"
  }else if(length(genes_early) != 0 & all(is.na(genes_late))){
    message("All genes from the late cell cycle signature are missing.
            \nSamples will not be ordered by cell cycle...")
    sample_order = NULL
    ha1b = NULL
    row_split = NULL
    row_title_cc = "Early cell cycle"
  }else{ #both are missing
    message("All genes from the late/early cell cycle signatures are missing.
            \nSamples will not be ordered by cell cycle...")
    sample_order = NULL
    ha1b = NULL
    row_split = NULL
    row_title_cc = NULL
  }
  
  #heatmap colors
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#4DF76F", "black", "#F74D4D"))
  
  #draw heatmap- late/early
  hm_pred_lateearly = Heatmap(this_data[genes_cc,, drop = FALSE],
                              top_annotation = if(!is.null(custom_annotation)) {
                                c(hm_a_predictions, custom_annotation)
                              } else {
                                hm_a_predictions
                              },
                              bottom_annotation = hm_a_lateearly,
                              name = "hm_early_late",
                              col = col_fun,
                              column_split = split, 
                              column_title = NULL,
                              row_split = row_split,
                              row_title = row_title_cc,
                              row_title_rot = 0,
                              cluster_row_slices = FALSE,
                              cluster_column_slices = FALSE,
                              cluster_columns = FALSE,
                              column_order = sample_order,
                              row_names_side = "left",
                              show_column_names = FALSE,
                              show_row_names = FALSE,
                              show_row_dend = FALSE,
                              border = TRUE,
                              row_names_gp = gpar(fontsize = plot_font_row_size),
                              row_title_gp = gpar(fontsize = plot_font_size),
                              clustering_distance_rows = "spearman",
                              clustering_method_rows = "ward.D2",
                              height = 8,
                              border_gp = gpar(lwd = 1),
                              show_heatmap_legend = FALSE)
  
  ################################## LUMINAL DIFFERENTIATION TFs ###################################
  genes_ud = genes_to_plot$UroDiff
  genes_ud = genes_ud[which(genes_ud %in% rownames(this_data))]
  
  hm_luminal_tfs = Heatmap(this_data[genes_ud,,drop = FALSE],
                           name = "hm_luminal_diff_tfs",
                           col = col_fun,
                           column_split = split,
                           column_title = NULL,
                           cluster_row_slices = FALSE,
                           cluster_column_slices = FALSE,
                           cluster_columns = FALSE,
                           column_order = sample_order,
                           row_names_side = "left",
                           row_title = "Luminal Differentiation TFs",
                           show_column_names = FALSE,
                           show_row_names = TRUE,
                           show_row_dend = FALSE,
                           border = TRUE,
                           row_names_gp = gpar(fontsize = plot_font_row_size),
                           row_title_gp = gpar(fontsize = plot_font_size),
                           cluster_rows = FALSE,
                           border_gp = gpar(lwd = 1),
                           show_heatmap_legend = FALSE,
                           row_title_rot = 0)
  
  ################################# LUMINAL DIFFERENTIATION GENES ##################################
  genes_upk = genes_to_plot$UPKs
  genes_upk = genes_upk[which(genes_upk %in% rownames(this_data))]
  
  hm_luminal_gen = Heatmap(this_data[genes_upk,,drop = FALSE],
                           name = "hm_luminal_diff_genes",
                           row_title = "Luminal Differentiation Genes",
                           column_title = NULL,
                           col = col_fun,
                           column_split = split,
                           cluster_row_slices = FALSE,
                           cluster_column_slices = FALSE,
                           cluster_columns = FALSE,
                           column_order = sample_order,
                           row_names_side = "left",
                           show_column_names = FALSE,
                           show_row_names = TRUE,
                           show_row_dend = TRUE,
                           border = TRUE,
                           row_names_gp = gpar(fontsize = plot_font_row_size),
                           row_title_gp = gpar(fontsize = plot_font_size),
                           cluster_rows = FALSE,
                           border_gp = gpar(lwd = 1),
                           show_heatmap_legend = FALSE,
                           row_title_rot = 0)
  
  ########################################## GU VS. URO ############################################
  genes_circ = genes_to_plot$Circuit
  genes_circ = genes_circ[which(genes_circ %in% rownames(this_data))]
  
  #check if all genes are present
  if(length(genes_circ) == 5){
    
    #calculate circuit sore
    circuit_score = apply(this_data, 2, function(col) sum(col[c("RB1", "FGFR3", "CCND1")]) - sum(col[c("E2F3", "CDKN2A")]))
    
    #add scores to genes_to_plot object
    genes_to_plot$Circuit_score = circuit_score
    
    #generate color palette
    col_fun_circ = circlize::colorRamp2(c(quantile(circuit_score, 0.10),
                                          median(circuit_score),
                                          quantile(circuit_score, 0.90)),
                                        c("#2166AC","white", "#B2182B"))
    
    col = list(`Circuit Score` = col_fun_circ)
    
    #create annotation track for heatmap 4
    hm_a_circuitscore = HeatmapAnnotation(`Circuit Score` = genes_to_plot$Circuit_score,
                                          simple_anno_size = unit(4, "mm"),
                                          simple_anno_size_adjust = TRUE,
                                          annotation_name_side = "left",
                                          col = col,
                                          show_legend = FALSE,
                                          border = TRUE,
                                          annotation_name_gp = gpar(fontsize = plot_font_row_size))
    
    #draw heatmap 4 - circuit score
    hm_circuitscore = Heatmap(this_data[genes_circ,, drop = FALSE],
                              name = "hm_gu_vs_uro",
                              row_title = "GU vs. Uro",
                              column_title = NULL,
                              bottom_annotation = hm_a_circuitscore,
                              col = col_fun,
                              column_split = split,
                              cluster_row_slices = FALSE,
                              cluster_column_slices = FALSE,
                              cluster_columns = FALSE,
                              column_order = sample_order,
                              row_names_side = "left",
                              show_column_names = FALSE, 
                              show_row_names = TRUE,
                              show_row_dend = FALSE,
                              border = TRUE,
                              row_names_gp = gpar(fontsize = plot_font_row_size),
                              row_title_gp = gpar(fontsize = plot_font_size),
                              cluster_rows = FALSE,
                              border_gp = gpar(lwd = 1),
                              show_heatmap_legend = FALSE,
                              row_title_rot = 0)
    
  }else if(length(genes_circ) != 0){
    hm_a_circuitscore = NULL
    message("Genes involved in the circuit score are missing.
            \nCircuit score will not be calculated...")
    
    #draw heatmap 4 - circuit score
    hm_circuitscore = Heatmap(this_data[genes_circ,,drop = FALSE],
                              name = "hm_gu_vs_uro",
                              row_title = "GU vs. Uro",
                              column_title = NULL,
                              bottom_annotation = hm_a_circuitscore,
                              col = col_fun,
                              column_split = split,
                              cluster_row_slices = FALSE,
                              cluster_column_slices = FALSE,
                              cluster_columns = FALSE,
                              column_order = sample_order,
                              row_names_side = "left",
                              show_column_names = FALSE,
                              show_row_names = TRUE,
                              show_row_dend = FALSE,
                              border = TRUE,
                              row_names_gp = gpar(fontsize = plot_font_row_size),
                              row_title_gp = gpar(fontsize = plot_font_size),
                              cluster_rows = FALSE,
                              border_gp = gpar(lwd = 1),
                              show_heatmap_legend = FALSE,
                              row_title_rot = 0)
  }else{
    message("No genes involved in the circuit score are found.
            \nCircuit score will not be calculated...")
    hm_circuitscore = NULL
  }
  
  ###################################### FGFR3 SIGNATURE ###########################################
  genes_fgfr3 = genes_to_plot$FGFR3
  genes_fgfr3 = genes_fgfr3[which(genes_fgfr3 %in% rownames(this_data))]
  
  if(length(genes_fgfr3) != 0){
    hm_fgfr3 = Heatmap(this_data[genes_fgfr3,, drop = FALSE],
                       name = "hm_fgfr3_signature",
                       col = col_fun,
                       column_split = split,
                       column_title = NULL,
                       cluster_row_slices = FALSE,
                       cluster_column_slices = FALSE,
                       cluster_columns = FALSE,
                       column_order = sample_order,
                       row_names_side = "left",
                       show_column_names = FALSE,
                       show_row_names = FALSE,
                       show_row_dend = FALSE,
                       row_title = "FGFR3 Signature",
                       row_title_rot = 0,
                       border = TRUE,
                       row_names_gp = gpar(fontsize = plot_font_row_size),
                       row_title_gp = gpar(fontsize = plot_font_size),
                       clustering_distance_rows= "pearson",
                       clustering_method_rows = "ward.D2",
                       height = 4,
                       border_gp = gpar(lwd = 1),
                       show_heatmap_legend = FALSE)
  }else{
    message("FGFR3 signature genes was not found in the expresison matrix.
             \nNo heatmap will be generated for this gene...")
    hm_fgfr3 = NULL
  }
  
  ########################################## TP63 ##################################################
  if("TP63" %in% rownames(this_data)){
    tp63 = this_data["TP63",, drop = FALSE]
    rownames(tp63) = c("TP63")
    
    #draw heatm 4_1 - TP63
    hm_tp63 = Heatmap(tp63,
                      name = "hm_tp63",
                      row_title = "TP63",
                      column_title = NULL,
                      col = col_fun,
                      column_split = split,
                      cluster_row_slices = FALSE,
                      cluster_column_slices = FALSE,
                      cluster_columns = FALSE,
                      column_order = sample_order,
                      row_names_side = "left",
                      show_column_names = FALSE,
                      show_row_names = TRUE,
                      show_row_dend = FALSE,
                      border = TRUE,
                      row_names_gp = gpar(fontsize = plot_font_row_size),
                      row_title_gp = gpar(fontsize = plot_font_size),
                      border_gp = gpar(lwd = 1),
                      show_heatmap_legend = FALSE,
                      row_title_rot = 0)
  }else{
    message("TP63 was not found in the expresison matrix.
            \nNo heatmap will be generated for this gene...")
    hm_tp63 = NULL
  }
  
  ##################################### DEFINITION BaSq ############################################
  genes_basq = genes_to_plot$BaSq
  genes_basq = genes_basq[which(genes_basq %in% rownames(this_data))]
  
  #check if all genes are present
  if(length(genes_basq) == 4) {
    
    #calculate BaSq ratio
    basq_ratio = apply(this_data, 2, function(col) sum(col[c("KRT5", "KRT14", "FOXA1")]) - sum(col[c("GATA3")]))
    
    genes_to_plot$BaSq_ratio = basq_ratio
    
    #generate color palette
    col_fun_basq = circlize::colorRamp2(c(quantile(basq_ratio, 0.10),
                                          median(basq_ratio),
                                          quantile(basq_ratio, 0.90)),
                                        c("#2166AC","white", "#B2182B"))
    
    col = list(`Ba/Sq Ratio` = col_fun_basq)
    
    #draw annotation track for heatmap 6 - BaSq ratio
    hm_a_basq = HeatmapAnnotation(`Ba/Sq Ratio` = genes_to_plot$BaSq_ratio,
                                  simple_anno_size = unit(4, "mm"),
                                  simple_anno_size_adjust = TRUE,
                                  annotation_name_side = "left",
                                  col = col,
                                  show_legend = FALSE,
                                  border = TRUE,
                                  annotation_name_gp = gpar(fontsize = plot_font_row_size))
    
    #draw heatmap 6 - BaSq ratio
    hm_basq = Heatmap(this_data[genes_basq,, drop = FALSE],
                      name = "hm_basq_definition",
                      bottom_annotation = hm_a_basq,
                      col = col_fun,
                      column_split = split,
                      column_title = NULL,
                      row_title = "Definition Ba/Sq",
                      cluster_row_slices = FALSE,
                      cluster_column_slices = FALSE,
                      cluster_columns = FALSE,
                      column_order = sample_order,
                      row_names_side = "left",
                      show_column_names = FALSE,
                      show_row_names = TRUE,
                      show_row_dend = FALSE,
                      border = TRUE,
                      row_names_gp = gpar(fontsize = plot_font_row_size),
                      row_title_gp = gpar(fontsize = plot_font_size),
                      cluster_rows = FALSE,
                      border_gp = gpar(lwd = 1),
                      show_heatmap_legend = FALSE,
                      row_title_rot = 0)
    
  }else if(length(genes_basq) != 0){
    #draw annotation track for heatmap 6 - BaSq ratio
    hm_a_basq = NULL
    message("Genes involved in the BaSq score are missing.
            \nBaSq ratio will not be calculated...")
    
    #draw heatmap 6 - BaSq ratio
    hm_basq = Heatmap(this_data[genes_basq,,drop = FALSE],
                      name = "hm_basq_definition",
                      bottom_annotation = ha6,
                      col = col_fun,
                      column_split = split,
                      column_title = NULL,
                      row_title = "Definition Ba/Sq",
                      cluster_row_slices = FALSE,
                      cluster_column_slices = FALSE,
                      cluster_columns = FALSE,
                      column_order = sample_order,
                      row_names_side = "left",
                      show_column_names = FALSE,
                      show_row_names = TRUE,
                      show_row_dend = FALSE,
                      border = TRUE,
                      row_names_gp = gpar(fontsize = plot_font_row_size),
                      row_title_gp = gpar(fontsize = plot_font_size),
                      cluster_rows = FALSE,
                      border_gp = gpar(lwd = 1),
                      show_heatmap_legend = FALSE,
                      row_title_rot = 0)
  }else{
    message("Genes involved in the BaSq score are missing.
            \nNo heatmap forBaSq ratio will be generated...")
    hm_basq = NULL
    col_fun_basq = NULL
  }
  
  ################################  KERATINIZATION SIGNATURE #######################################
  genes_krt = genes_to_plot$Keratinization
  genes_krt = genes_krt[which(genes_krt %in% rownames(this_data))]
  
  if(length(genes_krt) != 0){
    hm_keratinization = Heatmap(this_data[genes_krt,,drop = FALSE],
                                name = "hm_keratinization_signature",
                                col = col_fun,
                                column_split = split,
                                column_title = NULL,
                                row_title = "Keratinization Signature",
                                cluster_row_slices = FALSE,
                                cluster_column_slices = FALSE,
                                cluster_columns = FALSE,
                                column_order = sample_order,
                                row_names_side = "left",
                                show_column_names = FALSE,
                                show_row_names = FALSE,
                                show_row_dend = FALSE,
                                row_title_rot = 0,
                                border = TRUE,
                                row_names_gp = gpar(fontsize = plot_font_row_size),
                                row_title_gp = gpar(fontsize = plot_font_size),
                                clustering_distance_rows = "pearson",
                                clustering_method_rows = "ward.D2",
                                height = 4,
                                border_gp = gpar(lwd = 1),
                                show_heatmap_legend = FALSE)
  }else{
    message("Genes involved in the Keratinization signatures are missing.
            \nNo heatmap for keratinization signature ratio will be generated...")
    hm_keratinization = NULL
  }
  
  ##################################### ERBB RECEPTORS #############################################
  genes_erbb = genes_to_plot$ERBB
  genes_erbb = genes_erbb[which(genes_erbb %in% rownames(this_data))]
  
  #check if all genes are present
  if(length(genes_erbb) == 3){
    
    #calculate erbb score
    erbb_score = apply(this_data, 2, function(col) sum(col[c("EGFR")]) - sum(col[c("ERBB2", "ERBB3")]))
    
    #add scores back to the genes_to_plot object
    genes_to_plot$ERBB_score = erbb_score
    
    #generate color palette
    col_fun_erbb = circlize::colorRamp2(c(quantile(erbb_score, 0.10),
                                          median(erbb_score),
                                          quantile(erbb_score, 0.90)),
                                        c("#2166AC","white", "#B2182B"))
    
    col = list(`ERBB Score` = col_fun_erbb)
    
    #draw annotation track for ERBB score heatmap (heatmap 10)
    hm_a_erbscore = HeatmapAnnotation(`ERBB Score` = genes_to_plot$ERBB_score,
                                      simple_anno_size = unit(4, "mm"),
                                      simple_anno_size_adjust = TRUE,
                                      annotation_name_side = "left",
                                      col = col,
                                      annotation_legend_param = list(title = "Scores"),
                                      show_legend = FALSE,
                                      border = TRUE,
                                      annotation_name_gp = gpar(fontsize = plot_font_row_size))
    
    #draw heatmap 10, ERBB scores
    hm_erbscore = Heatmap(this_data[genes_erbb,,drop = FALSE],
                          name = "hm_erbb_receptors",
                          bottom_annotation = hm_a_erbscore,
                          col = col_fun,
                          column_split = split,
                          column_title = NULL,
                          row_title = "ERBB Receptors",
                          cluster_row_slices = FALSE,
                          cluster_rows = FALSE,
                          cluster_column_slices = FALSE,
                          cluster_columns = FALSE,
                          column_order = sample_order,
                          row_names_side = "left",
                          show_column_names = FALSE,
                          show_row_names = TRUE,
                          show_row_dend = FALSE,
                          border = TRUE,
                          row_names_gp = gpar(fontsize = plot_font_row_size),
                          row_title_gp = gpar(fontsize = plot_font_size),
                          clustering_distance_rows= "pearson",
                          clustering_method_rows = "ward.D2",
                          border_gp = gpar(lwd = 1),
                          show_heatmap_legend = FALSE,
                          row_title_rot = 0)
    
  }else if(length(genes_erbb) != 0){ 
    message("Some genes involved in the ERBB score are missing.
            \nERBB score will not be calculated.")
    
    #draw annotation track for ERBB score heatmap (heatmap 10)
    hm_a_erbscore = NULL
    
    #draw heatmap 10, ERBB scores
    hm_erbscore = Heatmap(this_data[genes_erbb,,drop = FALSE],
                          name = "hm_erbb_receptors",
                          bottom_annotation = ha10,
                          col = col_fun,
                          column_title = NULL,
                          row_title = "ERBB Receptors",
                          column_split = split,
                          cluster_row_slices = FALSE,
                          cluster_rows = FALSE,
                          cluster_column_slices = FALSE,
                          cluster_columns = FALSE,
                          column_order = sample_order,
                          row_names_side = "left",
                          show_column_names = FALSE,
                          show_row_names = TRUE,
                          show_row_dend = FALSE,
                          border = TRUE,
                          row_names_gp = gpar(fontsize = plot_font_row_size),
                          row_title_gp = gpar(fontsize = plot_font_size),
                          clustering_distance_rows= "pearson",
                          clustering_method_rows = "ward.D2",
                          border_gp = gpar(lwd = 1),
                          show_heatmap_legend = FALSE,
                          row_title_rot = 90)
  }else{
    message("Genes involved in the ERBB scores are missing.
            \nNo heatmap associated with ERBB will be generated...")
    hm_erbscore = NULL
  }
  
  ################################# CELL ADHESION ##################################################
  genes_ad = genes_to_plot$Adhesion
  genes_ad = genes_ad[which(genes_ad %in% rownames(this_data))]
  
  hm_adhesion = Heatmap(this_data[genes_ad,,drop = FALSE],
                        name = "hm_cell_adhesion",
                        col = col_fun,
                        column_split = split,
                        column_title = NULL,
                        row_title = "Cell Adhesion",
                        cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        cluster_columns = FALSE,
                        column_order = sample_order,
                        row_names_side = "left",
                        show_column_names = FALSE,
                        show_row_names = TRUE,
                        show_row_dend = FALSE,
                        border = TRUE,
                        row_names_gp = gpar(fontsize = plot_font_row_size),
                        row_title_gp = gpar(fontsize = plot_font_size),
                        cluster_rows = FALSE,
                        border_gp = gpar(lwd = 1),
                        show_heatmap_legend = FALSE,
                        row_title_rot = 0)
  
  ####################################### MYC TFs ##################################################
  genes_myc = genes_to_plot$MYC
  genes_myc = genes_myc[which(genes_myc %in% rownames(this_data))]
  
  hm_myc = Heatmap(this_data[genes_myc,,drop = FALSE],
                   name ="hm_myc_tfs",
                   col = col_fun,
                   column_split = split,
                   column_title = NULL,
                   row_title = "MYC TFs",
                   cluster_rows = FALSE,
                   cluster_row_slices = FALSE,
                   cluster_column_slices = FALSE,
                   cluster_columns = FALSE,
                   column_order = sample_order,
                   row_names_side = "left",
                   show_column_names = FALSE,
                   show_row_names = TRUE,
                   show_row_dend = FALSE,
                   border = TRUE,
                   row_names_gp = gpar(fontsize = plot_font_row_size),
                   row_title_gp = gpar(fontsize = plot_font_size),
                   border_gp = gpar(lwd = 1),
                   show_heatmap_legend = FALSE,
                   row_title_rot = 0)
  
  ################################ NEURONAL CELL MARKERS ###########################################
  genes_immune = intersect(rownames(this_data),genes_to_plot$Immune141_UP)
  genes_stroma = intersect(rownames(this_data),genes_to_plot$Stromal141_UP)
  
  #immune
  if(length(genes_immune) != 0){
    
    #calculate immune score
    immune_score = apply(this_data[genes_immune,],2,median)
    
    #add immune scores back to genes_to_plot object
    genes_to_plot$Immune141_UP_score = immune_score
    
    #generate color palette
    col_fun_immune = circlize::colorRamp2(c(quantile(immune_score, 0.10),
                                            median(immune_score),
                                            quantile(immune_score, 0.90)),
                                          c("#2166AC","white", "#B2182B"))
  }else{
    message("Genes involved in the Immune141_UP score are missing.
            \nImmune score will not be calculated.")
    immune_score = NULL
    col_fun_immune = NULL
  }
  
  #stroma
  if(length(genes_stroma) != 0){
    
    #calculate stromal score
    stromal_score <- apply(this_data[genes_stroma,],2,median)
    
    #add stromal scores back to genes_to_plot object
    genes_to_plot$Stromal141_UP_score = stromal_score
    
    #generate color palette
    col_fun_stromal = circlize::colorRamp2(c(quantile(stromal_score, 0.10),
                                             median(stromal_score),
                                             quantile(stromal_score, 0.90)),
                                           c("#2166AC","white", "#B2182B"))
  }else{
    message("Genes involved in the Stromal141_UP score are missing.
            \nStromal score will not be calculated.")
    stromal_score = NULL
    col_fun_stromal = NULL
  }
  
  col = list(`Immune 141_UP` = col_fun_immune,
             `Stromal 141_UP` = col_fun_stromal)
  
  #draw annotation track for heatmap 11 - SnNE and immune/stromal scores
  hm_a_immune_stroma = HeatmapAnnotation(`Immune 141_UP` = genes_to_plot$Immune141_UP_score,
                                         `Stromal 141_UP` = genes_to_plot$Stromal141_UP_score,
                                         simple_anno_size = unit(4, "mm"),
                                         simple_anno_size_adjust = TRUE,
                                         gap = unit(2, "mm"),
                                         annotation_name_side = "left",
                                         col = col,
                                         show_legend = FALSE,
                                         border = TRUE,
                                         annotation_name_gp = gpar(fontsize = plot_font_row_size))
  
  #ScNE genes
  genes_ne = genes_to_plot$ScNE
  genes_ne = genes_ne[which(genes_ne %in% rownames(this_data))]
  
  #draw heatmap 11 - ScNE, immune scores and stromal scores
  hm_neuronal = Heatmap(this_data[genes_ne,,drop = FALSE],
                        name = "hm_neuronal_cell_markers",
                        col = col_fun,
                        column_split = split,
                        column_title = NULL,
                        row_title = "Neuronal Cell Markers",
                        bottom_annotation = hm_a_immune_stroma,
                        cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        cluster_columns = FALSE,
                        column_order = sample_order, 
                        row_names_side = "left",
                        show_column_names = FALSE,
                        show_row_names = TRUE,
                        show_row_dend = FALSE,
                        border = TRUE,
                        gap = unit(10, "mm"),
                        row_names_gp = gpar(fontsize = plot_font_row_size),
                        row_title_gp = gpar(fontsize = plot_font_size),
                        cluster_rows = FALSE,
                        border_gp = gpar(lwd = 1),
                        show_heatmap_legend = show_hm_legend,
                        row_title_rot = 0)
  
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
  }else{
    message("No out_path provided, displaying heatmap in R session...")
  }
  
  ################################ SIGNATURE SCORES HEATMAP ####################################
  signature_heatmaps <- NULL
  
  if(plot_signature_scores) {
    message("Adding signature scores heatmap...")
    
    #get immune names
    immune_names <- c("b_cells", "t_cells",
                      "t_cells_cd8", "nk_cells",
                      "cytotoxicity_score", "neutrophils",
                      "monocytic_lineage", "macrophages",
                      "m2_macrophage", "myeloid_dendritic_cells")
      
    #get stromal names
    stromal_names <- c("endothelial_cells",
                       "fibroblasts", "smooth_muscle")  
    
    immune_labels = c("B Cells", "T Cells", "CD8+ T Cells", "NK Cells", 
                      "Cytotoxicity Score", "Neutrophils", "Monocytic Lineage", 
                      "Macrophages", "M2 Macrophages", "Myeloid DCs")
    
    stromal_labels = c("Endothelial Cells", "Fibroblasts", "Smooth Muscle")
      
    scores_data = these_predictions$scores
    
    #color function for immune/stromal signatures
    col_fun_immune_sig = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
    
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
    
    #check if columns exist in the scores data
    available_immune_names <- intersect(immune_names, colnames(scores_data))
    available_stromal_names <- intersect(stromal_names, colnames(scores_data))
    
    #build heatmap annotation (top)
    hm_a_sign <- HeatmapAnnotation(`Mol. grade (WHO1999)` = these_predictions$scores$molecular_grade_who_1999,
                                   `Mol. grade (WHO2022)` = these_predictions$scores$molecular_grade_who_2022,
                                   `Progression Risk` = these_predictions$scores$progression_risk,
                                   `Progression Score` = these_predictions$scores$progression_score,
                                   `Proliferation Score` = these_predictions$scores$proliferation_score,
                                   annotation_name_side = "left",
                                   show_legend = show_ann_legend,
                                   annotation_name_gp = gpar(fontsize = plot_font_size),
                                   border = TRUE,
                                   col = colour_obj,
                                   column_title = NULL,
                                   row_title = NULL)

    #immune scores heatmap
    if(length(available_immune_names) > 0) {
      hm_immune_scores_sig <- Heatmap(
        t(scale(scores_data[, available_immune_names, drop = FALSE])),
        top_annotation = hm_a_sign,
        column_order = sample_order,
        name = "Immune Signatures",
        border = TRUE,
        column_split = split,
        col = col_fun_immune_sig,
        cluster_rows = FALSE,
        show_heatmap_legend = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = plot_font_row_size),
        show_column_names = FALSE,
        column_title = NULL,
        row_title = NULL,
        row_labels = immune_labels[1:length(available_immune_names)]
      )
      
      if(is.null(signature_heatmaps)) {
        signature_heatmaps <- hm_immune_scores_sig
      } else {
        signature_heatmaps <- signature_heatmaps %v% hm_immune_scores_sig
      }
    }
    
    #stroma scores heatmap
    if(length(available_stromal_names) > 0) {
      hm_stroma_scores_sig <- Heatmap(
        t(scale(scores_data[, available_stromal_names, drop = FALSE])),
        column_order = sample_order,
        name = "Stromal Signatures",
        border = TRUE,
        column_split = split,
        col = col_fun_immune_sig,
        cluster_rows = FALSE,
        show_heatmap_legend = show_hm_legend,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = plot_font_row_size),
        show_column_names = FALSE,
        column_title = NULL,
        row_title = NULL,
        row_labels = stromal_labels[1:length(available_stromal_names)]
      )
      
      if(is.null(signature_heatmaps)) {
        signature_heatmaps <- hm_stroma_scores_sig
      } else {
        signature_heatmaps <- signature_heatmaps %v% hm_stroma_scores_sig
      }
    }
  }
  
  ################################ COMBINE SELECTED HEATMAP OBJECTS ############################
  #create list of all available heatmaps
  all_heatmaps <- list(
    early_late_cc = hm_pred_lateearly,
    luminal_tfs = hm_luminal_tfs,
    luminal_genes = hm_luminal_gen,
    fgfr3 = hm_fgfr3,
    circuit = hm_circuitscore,
    tp63 = hm_tp63,
    basq = hm_basq,
    keratinization = hm_keratinization,
    erbb = hm_erbscore,
    adhesion = hm_adhesion,
    myc = hm_myc,
    neuronal = hm_neuronal
  )
  
  #filter heatmaps based on user selection and availability
  selected_heatmaps <- list()
  
  for(section_name in names(include_sections)) {
    if(include_sections[[section_name]] && 
       section_name %in% names(all_heatmaps) && 
       !is.null(all_heatmaps[[section_name]])) {
      selected_heatmaps[[section_name]] <- all_heatmaps[[section_name]]
    }
  }
  
  #combine selected heatmaps
  if(length(selected_heatmaps) > 0) {
    main_heatmaps <- Reduce(`%v%`, selected_heatmaps)
  } else {
    stop("No valid heatmap sections selected or available...")
  }
  
  #add signature heatmaps if requested
  if(plot_signature_scores && !is.null(signature_heatmaps)) {
    final_hm = draw(main_heatmaps %v% signature_heatmaps,
                    column_title_gp = gpar("fontface", fontsize = 24),
                    column_title = plot_title)
  } else {
    final_hm = draw(main_heatmaps,
                    column_title_gp = gpar("fontface", fontsize = 24),
                    column_title = plot_title)
  }
  
  hm_sample_order <- column_order(final_hm@ht_list$hm_early_late)
  
  #only close graphics device if we opened one for file output
  if(save_to_file){
    dev.off()
    message("Heatmap saved to: ", paste0(out_path, plot_title, "_heatmap_scores.", out_format))
  }
  
  invisible(hm_sample_order)
}
