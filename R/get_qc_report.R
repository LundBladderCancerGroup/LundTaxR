#' @title Quality Control Report for Classification Results
#'
#' @description Summarise how reliable the classification results are for a cohort: which samples
#' could not be classified, which calls are low confidence, and how many of the genes needed by
#' the classifiers and signatures were found in the data.
#'
#' @details This function takes the output from [LundTaxR::classify_samples()] and collects quality
#' control information that is otherwise spread over the returned object, or not returned at all.
#' The report has the following parts:
#' \itemize{
#'   \item `samples`: one row per sample with the 5-class and 7-class calls, the highest 5-class
#'   prediction score, the prediction deltas (difference between the highest and second highest
#'   prediction score), and flags for samples that were not classified, low confidence calls
#'   (`prediction_delta_collapsed` below `delta_threshold`) and ties. If expression data is
#'   available, the fraction of classifier rules that could not be evaluated in each sample is also
#'   reported. Samples where this fraction is above `impute_reject` (see
#'   [LundTaxR::classify_samples()]) are not classified.
#'   \item `signature_coverage`: for each signature score (proliferation, progression, immune and
#'   the 141_UP stromal/immune scores), the number of expected genes and how many of them were
#'   missing from the data.
#'   \item `classifier_coverage` (only if expression data is available): for the 5-class and
#'   7-class subtype classifiers and the two grade classifiers, how many of their genes and rules
#'   are present in the data.
#'   \item `data_check` (only if expression data is available): the result of the log2
#'   transformation check that is also run by [LundTaxR::classify_samples()].
#' }
#' If expression data is available, `samples` also describes each sample's expression
#' distribution (`expr_median`, `expr_q25`, `expr_q75`, `expr_iqr` and `pct_zero`, the percentage of
#' values that are exactly zero). Samples whose median is more than `outlier_mads` median absolute
#' deviations away from the cohort's median of sample medians are flagged as
#' `distribution_outlier`. Since sample medians are almost identical in a well normalized cohort,
#' the allowed difference is never smaller than a quarter of the cohort's median within-sample IQR. This catches samples that were not normalized or transformed like the
#' rest of the cohort, or that are of low quality. The subtype calls are based on within-sample
#' gene pairs and are robust to such shifts, but the immune and stromal scores are not.
#' Expression data is taken from `this_data` if provided, otherwise from the `data` element of
#' `these_predictions` (i.e. [LundTaxR::classify_samples()] run with `include_data = TRUE`). The
#' data must use HGNC symbols as row names, which is the case for the `data` element even when
#' the classifier was run on Ensembl IDs.
#'
#' The default `delta_threshold` is a heuristic for flagging calls worth a second look, not a
#' validated cut-off.
#'
#' @param these_predictions Required parameter. Output from [LundTaxR::classify_samples()], run
#' with `include_pred_scores = TRUE` (default).
#' @param this_data Optional parameter. Expression data (genes as rows, HGNC symbols as row names,
#' samples as columns) used for the classification. Only needed if `these_predictions` does not
#' include the data.
#' @param delta_threshold Numeric value between 0 and 1. Samples with a collapsed prediction delta
#' below this value are flagged as low confidence. Default is 0.2.
#' @param outlier_mads Number of median absolute deviations a sample's median expression can differ
#' from the cohort before it is flagged as a distribution outlier. Default is 3.
#'
#' @return A list of class `lundtax_qc` with the elements `summary`, `samples`,
#' `signature_coverage`, `classifier_coverage` and `data_check`. Elements that cannot be computed
#' with the available input are NULL. Printing the object shows a short summary.
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
#' #print summary
#' sjodahl_qc
#'
#' #samples flagged as low confidence
#' subset(sjodahl_qc$samples, low_confidence)
#'
get_qc_report = function(these_predictions = NULL,
                         this_data = NULL,
                         delta_threshold = 0.2,
                         outlier_mads = 3){

  #checks
  if(is.null(these_predictions) || !is.list(these_predictions)){
    stop("these_predictions must be the output from classify_samples()...")
  }

  if(is.null(these_predictions$subtype_scores)){
    stop("these_predictions has no prediction scores. Run classify_samples() with include_pred_scores = TRUE...")
  }

  if(!is.numeric(delta_threshold) || length(delta_threshold) != 1 || delta_threshold < 0 || delta_threshold > 1){
    stop("delta_threshold must be a single number between 0 and 1...")
  }

  #get expression data, if available
  if(is.null(this_data)){
    this_data = these_predictions$data
  }
  if(!is.null(this_data)){
    this_data = as.matrix(this_data)
  }

  scores = these_predictions$subtype_scores
  sample_ids = names(these_predictions$predictions_5classes)

  #sample level information
  class_scores = scores[, c("Uro", "GU", "BaSq", "Mes", "ScNE"), drop = FALSE]
  classified = !is.na(these_predictions$predictions_5classes)
  top_score = apply(class_scores, 1, function(x) if(all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  delta = scores[, "prediction_delta_collapsed"]

  samples = data.frame(sample_id = sample_ids,
                       subtype_5_class = unname(these_predictions$predictions_5classes),
                       subtype_7_class = unname(these_predictions$predictions_7classes[sample_ids]),
                       top_score_5_class = unname(top_score),
                       prediction_delta_5_class = unname(scores[, "prediction_delta_5_class"]),
                       prediction_delta_7_class = unname(scores[, "prediction_delta_7_class"]),
                       prediction_delta_collapsed = unname(delta),
                       classified = unname(classified),
                       low_confidence = unname(classified & !is.na(delta) & delta < delta_threshold),
                       tie = unname(classified & !is.na(delta) & delta == 0),
                       stringsAsFactors = FALSE)

  #classifier coverage, only possible with expression data
  classifier_coverage = NULL
  data_check = NULL

  if(!is.null(this_data)){

    #gene availability per sample (gene present in the data and not NA in that sample)
    available = function(genes){
      out = matrix(FALSE,
                   nrow = length(genes),
                   ncol = length(sample_ids),
                   dimnames = list(genes, sample_ids))
      present = genes[genes %in% rownames(this_data)]
      out[present, ] = !is.na(this_data[present, sample_ids, drop = FALSE])
      out
    }

    #fraction of rules (gene pairs) that cannot be evaluated in each sample
    missing_rules = function(gene1, gene2){
      genes = unique(c(gene1, gene2))
      avail = available(genes)
      complete = avail[gene1, , drop = FALSE] & avail[gene2, , drop = FALSE]
      1 - colMeans(complete)
    }

    rules_5c = classifier_lundtax_5c$RF_scheme$rules
    rules_7c = classifier_lundtax_7c$RF_scheme$rules

    samples$missing_rules_5_class = unname(missing_rules(rules_5c$gene1, rules_5c$gene2))
    samples$missing_rules_7_class = unname(missing_rules(rules_7c$gene1, rules_7c$gene2))

    #coverage per classifier (based on genes present in the data, regardless of NA values)
    coverage_row = function(name, genes, gene1, gene2){
      data.frame(classifier = name,
                 genes_total = length(unique(genes)),
                 genes_present = sum(unique(genes) %in% rownames(this_data)),
                 rules_total = length(gene1),
                 rules_complete = sum(gene1 %in% rownames(this_data) & gene2 %in% rownames(this_data)),
                 stringsAsFactors = FALSE)
    }

    classifier_coverage = rbind(
      coverage_row("lundtax_5_class", classifier_lundtax_5c$RF_scheme$genes, rules_5c$gene1, rules_5c$gene2),
      coverage_row("lundtax_7_class", classifier_lundtax_7c$RF_scheme$genes, rules_7c$gene1, rules_7c$gene2),
      coverage_row("grade_who_1999", as.vector(classifier_grade_who_1999$TSPs),
                   classifier_grade_who_1999$TSPs[, "gene1"], classifier_grade_who_1999$TSPs[, "gene2"]),
      coverage_row("grade_who_2022", as.vector(classifier_grade_who_2022$TSPs),
                   classifier_grade_who_2022$TSPs[, "gene1"], classifier_grade_who_2022$TSPs[, "gene2"]))

    classifier_coverage$pct_rules_complete = round(100 * classifier_coverage$rules_complete / classifier_coverage$rules_total, 1)

    #expression distribution per sample
    expr = this_data[, sample_ids, drop = FALSE]
    quants = apply(expr, 2, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    samples$expr_median = unname(quants[2, ])
    samples$expr_q25 = unname(quants[1, ])
    samples$expr_q75 = unname(quants[3, ])
    samples$expr_iqr = samples$expr_q75 - samples$expr_q25
    samples$pct_zero = unname(round(100 * colMeans(expr == 0, na.rm = TRUE), 2))

    #flag samples with a median far from the rest of the cohort. In well normalized cohorts the 
    #sample medians are almost identical, so the allowed shift is never smaller than a quarter of 
    #the typical within-sample IQR, to avoid flagging trivial differences
    cohort_median = median(samples$expr_median, na.rm = TRUE)
    cohort_mad = mad(samples$expr_median, na.rm = TRUE)
    min_shift = 0.25 * median(samples$expr_iqr, na.rm = TRUE)
    allowed_shift = max(outlier_mads * cohort_mad, min_shift, na.rm = TRUE)
    if(!is.finite(allowed_shift) || allowed_shift == 0){
      samples$distribution_outlier = FALSE
    }else{
      samples$distribution_outlier = abs(samples$expr_median - cohort_median) > allowed_shift
    }
    attr(samples, "outlier_limits") = c(lower = cohort_median - allowed_shift,
                                        center = cohort_median,
                                        upper = cohort_median + allowed_shift)

    #log2 check
    data_check = int_check_log2_transformation(expression_df = this_data)
  }

  #signature coverage
  expected_genes = list(
    proliferation = signatures$proliferation$hgnc_symbol,
    progression = signatures$progression$hgnc_symbol,
    immune = signatures$immune$hgnc_symbol,
    score141up = signatures$signatures_plot$hgnc_symbol[signatures$signatures_plot$signature %in% c("Stromal141_UP", "Immune141_UP")])

  signature_coverage = NULL

  if(!is.null(this_data) || !is.null(these_predictions$na_genes)){
    signature_coverage = do.call(rbind, lapply(names(expected_genes), function(sig){
      genes = unique(na.omit(expected_genes[[sig]]))
      if(!is.null(this_data)){
        n_missing = sum(!genes %in% rownames(this_data))
      }else{
        na_genes = these_predictions$na_genes
        n_missing = length(intersect(genes, na_genes$genes[na_genes$signature == sig]))
      }
      data.frame(signature = sig,
                 genes_expected = length(genes),
                 genes_missing = n_missing,
                 pct_present = round(100 * (1 - n_missing / length(genes)), 1),
                 stringsAsFactors = FALSE)
    }))
  }

  #summary
  summary = data.frame(n_samples = nrow(samples),
                       n_classified = sum(samples$classified),
                       n_not_classified = sum(!samples$classified),
                       n_low_confidence = sum(samples$low_confidence),
                       n_ties = sum(samples$tie),
                       delta_threshold = delta_threshold,
                       stringsAsFactors = FALSE)

  if(!is.null(data_check)){
    summary$n_distribution_outliers = sum(samples$distribution_outlier, na.rm = TRUE)
    summary$outlier_mads = outlier_mads
    summary$log2_transformed = data_check$log2_transformed
  }

  report = list(summary = summary,
                samples = samples,
                signature_coverage = signature_coverage,
                classifier_coverage = classifier_coverage,
                data_check = data_check)

  class(report) = c("lundtax_qc", "list")

  return(report)
}

#' @title Print a LundTaxR QC Report
#'
#' @description Print method for the object returned by [LundTaxR::get_qc_report()].
#'
#' @param x Object returned by [LundTaxR::get_qc_report()].
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#'
#' @export
#'
print.lundtax_qc = function(x, ...){

  s = x$summary
  samples = x$samples

  list_samples = function(ids, n_max = 10){
    if(length(ids) > n_max){
      paste0(paste(ids[seq_len(n_max)], collapse = ", "), ", ... (", length(ids) - n_max, " more)")
    }else{
      paste(ids, collapse = ", ")
    }
  }

  cat(rep("=", 60), "\n", sep = "")
  cat("  LUNDTAX QC REPORT\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("  Samples:          ", s$n_samples, "\n", sep = "")
  cat("  Classified:       ", s$n_classified, "\n", sep = "")

  if(s$n_not_classified > 0){
    cat("  Not classified:   ", s$n_not_classified, " (", list_samples(samples$sample_id[!samples$classified]), ")\n", sep = "")
  }

  cat("  Low confidence:   ", s$n_low_confidence, " (prediction delta < ", s$delta_threshold, ")\n", sep = "")

  if(s$n_ties > 0){
    cat("  Ties:             ", s$n_ties, " (", list_samples(samples$sample_id[samples$tie]), ")\n", sep = "")
  }

  if(!is.null(s$n_distribution_outliers)){
    cat("  Expr. outliers:   ", s$n_distribution_outliers, sep = "")
    if(s$n_distribution_outliers > 0){
      cat(" (", list_samples(samples$sample_id[samples$distribution_outlier]), ")", sep = "")
    }
    cat("\n")
  }

  if(!is.null(x$data_check)){
    cat("  Data check:       ", gsub("\\s*\n\\s*", " ", x$data_check$message), "\n", sep = "")
  }

  if(!is.null(x$classifier_coverage)){
    cat(rep("-", 60), "\n", sep = "")
    cat("  Classifier rules with both genes present:\n")
    cc = x$classifier_coverage
    for(i in seq_len(nrow(cc))){
      cat(sprintf("    %-18s %5.1f%% (%d/%d)\n", cc$classifier[i], cc$pct_rules_complete[i], cc$rules_complete[i], cc$rules_total[i]))
    }
  }

  if(!is.null(x$signature_coverage)){
    cat(rep("-", 60), "\n", sep = "")
    cat("  Signature genes present:\n")
    sc = x$signature_coverage
    for(i in seq_len(nrow(sc))){
      cat(sprintf("    %-18s %5.1f%% (%d/%d)\n", sc$signature[i], sc$pct_present[i], sc$genes_expected[i] - sc$genes_missing[i], sc$genes_expected[i]))
    }
  }

  cat(rep("=", 60), "\n", sep = "")
  cat("  Per-sample details in $samples\n")

  invisible(x)
}
