#' @title Plot Ranked Signature Scores in Subtypes
#'
#' @description Return a point plot with ranked scores for a set variable.
#'
#' @details This function takes predictions returned with [LundTaxR::classify_samples()],
#' together with a score variable, subtype class and returns a point plot (grob) with ranked scores 
#' along the x axis and scores along the y axis. The returned plot will also color fill the points 
#' based on subtype classification.
#'
#' @param these_predictions Required parameter, if no data provided with `this_data`.
#' Predictions returned with [LundTaxR::classify_samples()].
#' @param this_data Optional parameter, makes it possible for the user to give the plotting function 
#' their own data, preferably retrieved with this function with return_data = TRUE. If provided, 
#' the plotting function will disregard any other parameters and only draw the plot using the data 
#' provided.
#' @param this_score Required parameter, should be one of the numeric columns in the scores data 
#' frame from the list returned with [LundTaxR::classify_samples()].
#' @param this_subtype Optional parameter, lets the user subset the return to a specific subtype.
#' If nopt provided, all subtypes within the selected class will be returned.
#' @param subtype_class Required, one of the following; 5_class or 7_class. Needed for coloring the 
#' points based on subtype classification.
#' @param seg_plot Boolean parameter, set to TRUE to return a segment plot with ranked scores for 
#' the selected subtype class. Default is FALSE.
#' @param seg_width Controls the width of the segments if `seg_plot = TRUE`.
#' @param add_stat Boolean parameter, set to TRUE to add regression lines, p value for each 
#' regression. Default is FALSE.
#' @param plot_title Required parameter, if `out_path` is specified. plot title, will also be pasted
#' to the exported file.
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can further specify the dimensions of the returned plot with `plot_width` and 
#' `plot_height`.
#' @param plot_width This parameter controls the width in inches. Default is 4 (1200 pixels at 300 
#' PPI).
#' @param plot_height This parameter controls the height in inches. Default is 4 (1200 pixels at 300
#' PPI).
#' @param return_data Boolean parameter, set to TRUE and return the formatted data used for 
#' plotting. Default is FALSE
#'
#' @return Nothing.
#' 
#' @import ggplot2 dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom ggpubr stat_cor
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
#' #Seg plot, all 5 classes                                   
#' plot_ranked_score(these_predictions = sjodahl_classes, 
#'                   this_score = "proliferation_score", 
#'                   subtype_class = "5_class", 
#'                   seg_plot = TRUE, 
#'                   seg_width = 1)
#'
#' #Seg plot, all 7 classes
#' plot_ranked_score(these_predictions = sjodahl_classes, 
#'                   this_score = "proliferation_score", 
#'                   subtype_class = "7_class", 
#'                   seg_plot = TRUE, 
#'                   seg_width = 1)
#'
#' #rank plot, proliferation score vs 5 class subtypes
#' plot_ranked_score(these_predictions = sjodahl_classes, 
#'                   this_score = "proliferation_score", 
#'                   subtype_class = "5_class")
#'
#' #rank plot with statis added
#' plot_ranked_score(these_predictions = sjodahl_classes, 
#'                   this_score = "proliferation_score", 
#'                   subtype_class = "5_class",
#'                   add_stat = TRUE,
#'                   this_subtype = "GU")
#'                   
#' #return data
#' ranked_proliferation = plot_ranked_score(these_predictions = sjodahl_classes, 
#'                                          this_score = "proliferation_score", 
#'                                          return_data = TRUE,
#'                                          subtype_class = "7_class")
#' #view data
#' head(ranked_proliferation)
#'
plot_ranked_score = function(these_predictions = NULL, 
                             this_data = NULL,
                             this_score = NULL, 
                             this_subtype = NULL,
                             subtype_class = NULL,
                             seg_plot = FALSE, 
                             seg_width = 0.1,
                             add_stat = FALSE,
                             out_path = NULL,
                             out_format = "png",
                             plot_width = 8,
                             plot_height = 5,
                             plot_title = NULL,
                             return_data = FALSE){
  
  if(is.null(this_data)){
    #checks
    if(is.null(these_predictions)){
      stop("these_predictions is missing, should be the return from lund_predict_sub...")
    }
    
    #impute plot title
    if(is.null(plot_title)){
      plot_title = this_score
    }
    
    if(is.null(this_score)){
      stop("this_score is missing, should be a numeric variable in the scores data frame within these_predictions...")
    }else if(!this_score %in% colnames(these_predictions$scores)){
      stop(paste0(this_score, " is non-existing in the scores data frame within these_predictions..."))
    }else if(!is.numeric(these_predictions$scores[,this_score])){
      stop(paste0(this_score, " is not of type numeric..."))
    }
    
    #subset on subtype
    if(subtype_class == "5_class"){
      subtypes = data.frame(these_predictions$predictions_5classes) %>% 
        tibble::rownames_to_column("sample_id")
      if(is.null(this_subtype)){
        this_subtype = c("Uro", "GU", "ScNE", "BaSq", "Mes")
      }else{
        if(!this_subtype %in% c("Uro", "GU", "ScNE", "BaSq", "Mes")){
          stop("The requested subtype is not available in the selected class..")
        }
      }
    }else if(subtype_class == "7_class"){
      subtypes = data.frame(these_predictions$predictions_7classes) %>% 
        tibble::rownames_to_column("sample_id") 
      if(is.null(this_subtype)){
        this_subtype = c("UroA","UroB","UroC", "GU", "ScNE", "BaSq", "Mes")
      }else{
        if(!this_subtype %in% c("UroA","UroB","UroC", "GU", "ScNE", "BaSq", "Mes")){
          stop("The requested subtype is not available in the selected class..")
        }
      }
    }else{
      stop("A valid subtype class must be provided, one of the following; 5_class or 7_class...")
    }
    
    #subset score, left join with subtypes, add rank
    this_data = these_predictions$scores %>% 
      dplyr::select(all_of(this_score)) %>% 
      tibble::rownames_to_column("sample_id") %>% 
      dplyr::left_join(subtypes, by = "sample_id") %>% 
      dplyr::rename(score = 2, subtype = 3) %>%
      dplyr::filter(subtype %in% this_subtype) %>% 
      dplyr::mutate_at(vars(subtype), factor) %>%
      mutate(rank = as.numeric(as.factor(score)))
    
    if(subtype_class == "5_class"){
      this_data = this_data %>% 
        mutate(subtype = factor(subtype, levels = c("Uro", "GU", "BaSq", "Mes", "ScNE")))
    }else if(subtype_class == "7_class"){
      this_data = this_data %>% 
        mutate(subtype = factor(subtype, levels = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")))
      
    }
    
    if(return_data){
      message("No plot generated, returning plot data instead...")
      return(this_data)
    }
  }else{
    
    if(!is.null(this_subtype)){
      message("You have provided your own data, the function will not subset plotting data on the requested subtype...")
    }
    
    if(is.null(this_score)){
      stop("You must specify the type of score with `this_score` in your provided dataset...")
    }
  }
  
  #get subtitle
  formatted_subtitle <- tools::toTitleCase(gsub("_", " ", this_score))
  
  if(seg_plot){
    #get segment data for each subtype
    gu_lines = this_data %>% filter(subtype == "GU") %>% 
      select(rank) %>% 
      arrange(rank) %>% 
      rename(x_start = rank) %>% 
      mutate(x_end = x_start) %>% 
      mutate(subtype = factor("GU")) 
    
    basq_lines = this_data %>% filter(subtype == "BaSq") %>% 
      select(rank) %>% 
      arrange(rank) %>% 
      rename(x_start = rank) %>% 
      mutate(x_end = x_start) %>% 
      mutate(subtype = factor("BaSq")) 
    
    mes_lines = this_data %>% filter(subtype == "Mes") %>% 
      select(rank) %>% 
      arrange(rank) %>% 
      rename(x_start = rank) %>% 
      mutate(x_end = x_start) %>% 
      mutate(subtype = factor("Mes")) 
    
    scne_lines = this_data %>% filter(subtype == "ScNE") %>% 
      select(rank) %>% 
      arrange(rank) %>% 
      rename(x_start = rank) %>% 
      mutate(x_end = x_start) %>% 
      mutate(subtype = factor("ScNE"))
    
    if(subtype_class == "5_class"){
      uro_lines = this_data %>% filter(subtype == "Uro") %>% 
        select(rank) %>% 
        arrange(rank) %>% 
        rename(x_start = rank) %>% 
        mutate(x_end = x_start) %>% 
        mutate(subtype = factor("Uro")) 
      
      #set factor levels
      this_data$subtype = factor(this_data$subtype, levels = c("Uro", "GU", "BaSq", "Mes", "ScNE"))
      
    }else if(subtype_class == "7_class"){
      uroa_lines = this_data %>% filter(subtype == "UroA") %>% 
        select(rank) %>% 
        arrange(rank) %>% 
        rename(x_start = rank) %>% 
        mutate(x_end = x_start) %>% 
        mutate(subtype = factor("UroA")) 
      
      urob_lines = this_data %>% filter(subtype == "UroB") %>% 
        select(rank) %>% 
        arrange(rank) %>% 
        rename(x_start = rank) %>% 
        mutate(x_end = x_start) %>% 
        mutate(subtype = factor("UroB")) 
      
      uroc_lines = this_data %>% filter(subtype == "UroC") %>% 
        select(rank) %>% 
        arrange(rank) %>% 
        rename(x_start = rank) %>% 
        mutate(x_end = x_start) %>% 
        mutate(subtype = factor("UroC")) 
      
      #set factor levels
      this_data$subtype = factor(this_data$subtype, levels = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE"))
    }
    
    #add subtype labels
    if(subtype_class == "5_class"){
      subtype_labels <- data.frame(
        subtype = c("Uro", "GU", "BaSq", "Mes", "ScNE"),
        y_pos = c(5.5, 4.5, 3.5, 2.5, 1.5),
        x_pos = rep(-30, 5)
      ) 
    }else if(subtype_class == "7_class"){
      subtype_labels <- data.frame(
        subtype = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE"),
        y_pos = c(7.5, 6.5, 5.5, 4.5, 3.5, 2.5, 1.5),
        x_pos = rep(-30, 7)
      ) 
    }
    
    #set the minor grid line color based on subtype_class
    if(subtype_class == "5_class") {
      minor_grid_color <- "white"
    } else if(subtype_class == "7_class") {
      minor_grid_color <- "black"
    }
    
    
    #segment plot
    #segment plot
    my_plot = ggplot() +
      {if(subtype_class == "5_class") geom_segment(data = uro_lines, aes(x = x_start, xend = x_end, y = 5, yend = 5.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "5_class") geom_segment(data = gu_lines, aes(x = x_start, xend = x_end, y = 4, yend = 4.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "5_class") geom_segment(data = basq_lines, aes(x = x_start, xend = x_end, y = 3, yend = 3.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "5_class") geom_segment(data = mes_lines, aes(x = x_start, xend = x_end, y = 2, yend = 2.98, color = subtype), linewidth = seg_width, lineend = "butt")} + 
      {if(subtype_class == "5_class") geom_segment(data = scne_lines, aes(x = x_start, xend = x_end, y = 1, yend = 1.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "7_class") geom_segment(data = uroa_lines, aes(x = x_start, xend = x_end, y = 7, yend = 7.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "7_class") geom_segment(data = urob_lines, aes(x = x_start, xend = x_end, y = 6, yend = 6.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "7_class") geom_segment(data = uroc_lines, aes(x = x_start, xend = x_end, y = 5, yend = 5.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "7_class") geom_segment(data = gu_lines, aes(x = x_start, xend = x_end, y = 4, yend = 4.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "7_class") geom_segment(data = basq_lines, aes(x = x_start, xend = x_end, y = 3, yend = 3.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "7_class") geom_segment(data = mes_lines, aes(x = x_start, xend = x_end, y = 2, yend = 2.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      {if(subtype_class == "7_class") geom_segment(data = scne_lines, aes(x = x_start, xend = x_end, y = 1, yend = 1.98, color = subtype), linewidth = seg_width, lineend = "butt")} +
      scale_color_manual(values = lund_colors$lund_colors) + 
      xlab("Index") + 
      geom_text(data = subtype_labels, 
                aes(x = x_pos, y = y_pos, label = subtype),
                hjust = 0, vjust = 0.5, size = 4) +
      ggtitle(label = "Segments Plot", subtitle = formatted_subtitle) +
      guides(colour = guide_legend(nrow = 1)) +
      coord_cartesian(xlim = c(-30, max(this_data$rank)), expand = FALSE) +
      scale_x_continuous(limits = c(-30, max(this_data$rank))) +
      {if(subtype_class == "5_class") scale_y_continuous(limits = c(1, 6))} +
      {if(subtype_class == "7_class") scale_y_continuous(limits = c(1, 8))} +
      theme(legend.position = "none", 
            legend.title = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line( size = 0.2, color = "black"),
            panel.grid.minor.y = element_line(size = 0.1, color = minor_grid_color),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            panel.border = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 10), 
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.ticks.y = element_blank(), 
            axis.line.x = element_line(),
            axis.text.y = element_blank())
  }else{
    #draw default plot
    my_plot = ggplot(data = this_data, aes(x = rank, y = score, color = subtype)) +
      geom_point(size = 2, shape = 16) +
      {if(add_stat) geom_smooth(method = lm, se = FALSE, linewidth = 0.5)} +
      {if(add_stat) stat_cor(method = "pearson", label.x = 10)} +
      scale_color_manual(values = lund_colors$lund_colors) + 
      xlab("Index") + 
      labs(color = "Subtype") +
      ylab(formatted_subtitle) +
      ggtitle(label = paste0("Ranked Signature Score")) +
      theme_bw() +
      theme(legend.position = "right", 
            axis.line.x = element_blank())
  }
  
  if(is.null(this_subtype)){
    file_name = paste0(subtype_class, "_", this_score)
  }else{
    file_name = paste0(this_subtype, "_", this_score)
  }
  
  if(seg_plot){
    file_name = paste0(file_name, "_seg_plot")
  }else{
    file_name = paste0(file_name, "_ranked_score")
  }
  
  if(!is.null(out_path)){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, file_name, ".pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, file_name, ".png"),
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
  }else{
    return(my_plot) 
  }
}
