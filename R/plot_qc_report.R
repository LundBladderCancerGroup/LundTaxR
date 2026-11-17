#' @title Plot a LundTaxR QC Report
#'
#' @description Visualize the output from [LundTaxR::get_qc_report()] in a multi-panel figure.
#'
#' @details The figure has up to four panels:
#' \itemize{
#'   \item **Prediction confidence**: the collapsed prediction delta (difference between the
#'   highest and second highest prediction score) for each sample, grouped and colored by subtype.
#'   The dashed line marks the `delta_threshold` used to flag low confidence calls.
#'   \item **Signature genes present**: the percentage of the expected genes found in the data for
#'   each signature score.
#'   \item **Classifier rules complete** (only if expression data was available for the report):
#'   the percentage of rules (gene pairs) with both genes present for the subtype and grade
#'   classifiers.
#'   \item **Expression distribution per sample** (only if expression data was available): one
#'   point per sample, showing its median expression against its interquartile range. Samples that
#'   were transformed or normalized differently, or have many zero values, separate from the main
#'   cloud. The dashed lines mark the outlier limits for the median used by
#'   [LundTaxR::get_qc_report()] (`outlier_mads`), and flagged samples are highlighted (and named, 
#'   if there are 15 or fewer).
#'   \item **Missing classifier rules per sample** (only if expression data was available): the
#'   fraction of 5-class classifier rules that could not be evaluated in each sample. Samples above
#'   the dashed `impute_reject` line are not classified.
#' }
#' The number of samples, unclassified samples and low confidence calls is given in the subtitle.
#'
#' @param these_qc Output from [LundTaxR::get_qc_report()]. Required if `these_predictions` is
#' not provided.
#' @param these_predictions Optional. Output from [LundTaxR::classify_samples()]. If provided
#' instead of `these_qc`, the QC report is generated internally.
#' @param this_data Optional. Expression data, passed to [LundTaxR::get_qc_report()] when
#' `these_predictions` is used.
#' @param subtype_class Can be one of the following; 5_class (default) or 7_class. Decides how
#' samples are grouped and colored in the prediction confidence panel.
#' @param delta_threshold Optional. Low confidence threshold, used when the QC report is generated
#' from `these_predictions`. Default is 0.2.
#' @param impute_reject The `impute_reject` value used for the classification, drawn as a reference
#' line in the missing rules panel. Default is 0.67 (the default in [LundTaxR::classify_samples()]).
#' @param plot_title Title of the figure. Default is "LundTax QC Report".
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' @param file_name File name (without extension) for the exported plot. Default is "qc_report".
#' @param plot_width This parameter controls the width in inches. Default is 10 (3000 pixels at 300
#' PPI).
#' @param plot_height This parameter controls the height in inches. Default is 11 (3300 pixels at 
#' 300 PPI).
#'
#' @return A patchwork object (if no `out_path` is provided), otherwise the plot is written to disk.
#'
#' @import ggplot2 patchwork
#' @importFrom grDevices dev.off pdf png
#'
#' @export
#'
#' @examples
#' #run classifier
#' sjodahl_classes = classify_samples(this_data = sjodahl_2017,
#'                                    include_data = TRUE,
#'                                    verbose = FALSE)
#'
#' #get QC report
#' sjodahl_qc = get_qc_report(these_predictions = sjodahl_classes)
#'
#' #plot
#' plot_qc_report(these_qc = sjodahl_qc)
#'
plot_qc_report = function(these_qc = NULL,
                          these_predictions = NULL,
                          this_data = NULL,
                          subtype_class = "5_class",
                          delta_threshold = 0.2,
                          impute_reject = 0.67,
                          plot_title = "LundTax QC Report",
                          out_path = NULL,
                          out_format = "png",
                          file_name = "qc_report",
                          plot_width = 10,
                          plot_height = 11){

  #checks
  if(is.null(these_qc)){
    if(is.null(these_predictions)){
      stop("Provide a QC report with `these_qc`, or the output from classify_samples() with `these_predictions`...")
    }
    these_qc = get_qc_report(these_predictions = these_predictions,
                             this_data = this_data,
                             delta_threshold = delta_threshold)
  }

  if(!inherits(these_qc, "lundtax_qc")){
    stop("these_qc must be the output from get_qc_report()...")
  }

  if(!subtype_class %in% c("5_class", "7_class")){
    stop("subtype_class must be one of the following; 5_class or 7_class...")
  }

  #shared styling
  bar_col = "#5B6B7F"
  ref_col = "grey30"
  qc_theme = theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 0.3, colour = "grey40"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey40"),
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8, colour = "grey20"))

  #################################### PREDICTION CONFIDENCE ######################################
  if(subtype_class == "5_class"){
    subtype_levels = c("Uro", "GU", "BaSq", "Mes", "ScNE")
    subtype_col = "subtype_5_class"
  }else{
    subtype_levels = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")
    subtype_col = "subtype_7_class"
  }

  conf_data = these_qc$samples[these_qc$samples$classified & !is.na(these_qc$samples$prediction_delta_collapsed), ]
  conf_data$subtype = factor(conf_data[[subtype_col]], levels = subtype_levels)
  conf_data = conf_data[!is.na(conf_data$subtype), ]
  conf_data$confidence = factor(ifelse(conf_data$low_confidence, "Low confidence", "Confident"), levels = c("Confident", "Low confidence"))

  #number of low confidence calls per subtype, shown above each group
  n_low = tapply(conf_data$low_confidence, conf_data$subtype, sum)
  n_all = tapply(conf_data$low_confidence, conf_data$subtype, length)
  group_labels = data.frame(subtype = factor(subtype_levels, levels = subtype_levels),
                            label = ifelse(is.na(n_all), "n = 0", paste0("n = ", n_all, "\n", n_low, " low")),
                            stringsAsFactors = FALSE)

  threshold = these_qc$summary$delta_threshold

  p_conf = ggplot(conf_data, aes(x = subtype, y = prediction_delta_collapsed)) +
    geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 0.4, colour = ref_col) +
    geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.3, colour = "grey45", fill = NA) +
    geom_jitter(aes(colour = subtype, shape = confidence), width = 0.18, height = 0, size = 1.6, alpha = 0.8) +
    geom_text(data = group_labels, aes(x = subtype, y = 1.13, label = label), inherit.aes = FALSE, size = 2.6, colour = "grey25", lineheight = 0.9) +
    scale_colour_manual(values = lund_colors$lund_colors, guide = "none", drop = FALSE) +
    scale_shape_manual(values = c(Confident = 16, `Low confidence` = 4), name = NULL, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, by = 0.25), expand = c(0, 0.01)) +
    labs(title = "Prediction confidence",
         subtitle = paste0("Difference between the top two prediction scores, dashed line at ", threshold),
         x = NULL,
         y = "Prediction delta") +
    qc_theme +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.margin = margin(0, 0, 0, 0))

  #################################### COVERAGE PANELS ############################################
  coverage_bar = function(df, label_col, value_col, n_present, n_total, title, subtitle){
    df$value = df[[value_col]]
    df$label = factor(df[[label_col]], levels = rev(df[[label_col]]))
    df$text = paste0(format(round(df$value, 1), nsmall = 1), "% (", n_present, "/", n_total, ")")
    ggplot(df, aes(x = value, y = label)) +
      geom_col(width = 0.6, fill = bar_col) +
      geom_text(aes(label = text), hjust = 0, nudge_x = 2, size = 2.6, colour = "grey20") +
      scale_x_continuous(breaks = seq(0, 100, by = 25), expand = c(0, 0)) +
      coord_cartesian(xlim = c(0, 145), clip = "off") +
      labs(title = title, subtitle = subtitle, x = "% present", y = NULL) +
      qc_theme +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
            axis.ticks.y = element_blank())
  }

  p_sig = NULL
  if(!is.null(these_qc$signature_coverage)){
    sc = these_qc$signature_coverage
    sc$signature = c(proliferation = "Proliferation", progression = "Progression",
                     immune = "Immune", score141up = "141_UP (stromal/immune)")[sc$signature]
    p_sig = coverage_bar(sc, "signature", "pct_present",
                         sc$genes_expected - sc$genes_missing, sc$genes_expected,
                         "Signature genes present", "Genes found in the data")
  }

  p_clf = NULL
  p_rules = NULL
  if(!is.null(these_qc$classifier_coverage)){
    cc = these_qc$classifier_coverage
    cc$classifier = c(lundtax_5_class = "LundTax 5-class", lundtax_7_class = "LundTax 7-class",
                      grade_who_1999 = "Grade WHO1999", grade_who_2022 = "Grade WHO2022")[cc$classifier]
    p_clf = coverage_bar(cc, "classifier", "pct_rules_complete",
                         cc$rules_complete, cc$rules_total,
                         "Classifier rules complete", "Gene pairs with both genes in the data")

    #missing rules per sample, samples without missing rules are left out (and counted in the subtitle)
    n_complete = sum(these_qc$samples$missing_rules_5_class == 0, na.rm = TRUE)
    rules_data = these_qc$samples[!is.na(these_qc$samples$missing_rules_5_class) & these_qc$samples$missing_rules_5_class > 0, ]
    rules_subtitle = paste0("5-class rules that could not be evaluated\n", n_complete, " of ",
                            nrow(these_qc$samples), " samples have no missing rules (not shown)")

    p_rules = ggplot(rules_data, aes(x = missing_rules_5_class)) +
      {if(nrow(rules_data) > 0) geom_histogram(binwidth = 0.025, boundary = 0, fill = bar_col, colour = "white", linewidth = 0.2)} +
      {if(nrow(rules_data) == 0) annotate("text", x = if(impute_reject > 0.5) impute_reject / 2 else (1 + impute_reject) / 2, y = 0.5, label = "No missing rules in any sample", size = 3.2, colour = "grey30")} +
      {if(nrow(rules_data) == 0) list(scale_y_continuous(limits = c(0, 1), breaks = NULL),
                                      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()))} +
      geom_vline(xintercept = impute_reject, linetype = "dashed", linewidth = 0.4, colour = ref_col) +
      annotate("text", x = impute_reject, y = Inf, label = paste0(" impute_reject = ", impute_reject),
               hjust = 0, vjust = 1.5, size = 2.6, colour = "grey25") +
      scale_x_continuous(breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
      {if(nrow(rules_data) > 0) scale_y_continuous(expand = expansion(mult = c(0, 0.1)))} +
      coord_cartesian(xlim = c(-0.025, 1.025)) +
      labs(title = "Missing classifier rules per sample",
           subtitle = rules_subtitle,
           x = "Fraction of rules missing",
           y = "Samples") +
      qc_theme
  }

  #################################### EXPRESSION DISTRIBUTION ####################################
  #one point per sample (median vs IQR), which stays readable for large cohorts
  p_dist = NULL
  if(!is.null(these_qc$samples$expr_median)){
    dist_data = these_qc$samples[!is.na(these_qc$samples$expr_median), ]
    dist_data$status = factor(ifelse(dist_data$distribution_outlier, "Outlier", "Within limits"),
                              levels = c("Within limits", "Outlier"))
    
    #draw outliers last so they end up on top
    dist_data = dist_data[order(dist_data$distribution_outlier), ]
    limits = attr(these_qc$samples, "outlier_limits")
    n_samples = nrow(dist_data)

    p_dist = ggplot(dist_data, aes(x = expr_median, y = expr_iqr, colour = status)) +
      geom_point(data = dist_data[!dist_data$distribution_outlier, ],
                 size = if(n_samples > 500) 0.9 else 1.6,
                 alpha = if(n_samples > 500) 0.5 else 0.8) +
      geom_point(data = dist_data[dist_data$distribution_outlier, ], size = 1.8) +
      scale_colour_manual(values = c(`Within limits` = "#5B6B7F", Outlier = "#B2182B"), name = NULL, drop = FALSE) +
      labs(title = "Expression distribution per sample",
           subtitle = "One point per sample, dashed lines mark the outlier limits for the median",
           x = "Median expression",
           y = "Interquartile range") +
      qc_theme +
      theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
            legend.position = "bottom",
            legend.text = element_text(size = 8),
            legend.margin = margin(0, 0, 0, 0)) +
      guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

    if(!is.null(limits) && all(is.finite(limits)) && limits[["upper"]] > limits[["lower"]]){
      p_dist = p_dist +
        geom_vline(xintercept = limits[c("lower", "upper")], linetype = "dashed", linewidth = 0.4, colour = ref_col)
    }

    #name the outliers, if there are not too many
    outliers = dist_data[dist_data$distribution_outlier, ]
    if(nrow(outliers) > 0 && nrow(outliers) <= 15){
      p_dist = p_dist +
        geom_text(data = outliers, aes(label = sample_id), hjust = -0.15, vjust = 0.5,
                  size = 2.5, colour = "grey20", check_overlap = TRUE, show.legend = FALSE)
    }else if(nrow(outliers) > 15){
      p_dist = p_dist +
        labs(caption = paste0(nrow(outliers), " outliers, see `$samples` for sample IDs")) +
        theme(plot.caption = element_text(size = 8, colour = "grey30"))
    }
  }

  #################################### COMBINE ####################################################
  s = these_qc$summary
  summary_text = paste0(s$n_samples, " samples  |  ",
                        s$n_classified, " classified  |  ",
                        s$n_not_classified, " not classified  |  ",
                        s$n_low_confidence, " low confidence")
  if(!is.null(these_qc$data_check)){
    summary_text = paste0(summary_text, "  |  ",
                          ifelse(these_qc$data_check$log2_transformed, "data looks log2-transformed", "data does NOT look log2-transformed"))
  }

  bottom_panels = Filter(Negate(is.null), list(p_sig, p_clf, p_rules))

  if(length(bottom_panels) == 0){
    combined = p_conf
  }else if(length(bottom_panels) == 1){
    combined = p_conf / bottom_panels[[1]] + plot_layout(heights = c(2.5, 1))
  }else{
    right_col = wrap_plots(bottom_panels[c(1, 2)], ncol = 1)
    bottom = if(length(bottom_panels) == 3) (right_col | bottom_panels[[3]]) else right_col
    if(is.null(p_dist)){
      combined = p_conf / bottom + plot_layout(heights = c(1.3, 1))
    }else{
      combined = p_conf / p_dist / bottom + plot_layout(heights = c(1.2, 0.9, 1))
    }
  }

  combined = combined +
    plot_annotation(title = plot_title,
                    subtitle = summary_text,
                    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                                  plot.subtitle = element_text(size = 10, colour = "grey30")))

  #export or return
  if(!is.null(out_path)){
    if(out_format == "pdf"){
      pdf(file.path(out_path, paste0(file_name, ".pdf")),
          width = plot_width,
          height = plot_height)
    }else if(out_format == "png"){
      png(file.path(out_path, paste0(file_name, ".png")),
          width = plot_width,
          height = plot_height,
          units = "in",
          res = 300,
          pointsize = 10,
          bg = "white")
    }else{
      stop("Enter a valid output format (pdf or png)...")
    }
    print(combined)
    dev.off()
    message(paste0("Plot exported to ", file.path(out_path, paste0(file_name, ".", out_format))))
  }else{
    return(combined)
  }
}
