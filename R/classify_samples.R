#' @title Lund Taxonomy Classifier
#'
#' @description Predict Lund Taxonomy subtypes based on rule-based Random Forest classifiers.
#'
#' @details This function uses 2 classifiers to classify the samples: 5-class
#' classifier first  classifies samples into Uro, GU, BaSq, Mes or ScNE.
#' Samples classified as Uro receive a second classification as UroA, B or C by
#' the second classifier. This function internally calls
#' [LundTaxR::int_calc_signatures()] for retrieving signature scores.
#'
#' @param this_data Required parameter. Data frame or matrix with expression values.
#' @param gene_id Specify the type of gene identifier used in `this_data`. Accepted values are;
#' hgnc_symbol (default) or ensembl_gene_id.
#' @param threshold_progression Threshold to flag a sample as high risk of progression, default is
#' 0.58.
#' @param threshold_grade Threshold to flag a sample as high grade, default is 0.5.
#' @param log_transform Boolean parameter. If TRUE, the function log transforms the incoming
#' expression values. Default is FALSE.
#' @param adjust Boolean parameter. If TRUE (default), the function will proceed with adjusting the
#' scores based on stable genes. If FALSE, no adjustment will be made and the original score values
#' will be retained.
#' @param adj_factor Only applicable if adjust is set to TRUE. Allows users to apply a proportional
#' adjustment to the normalized scores, enabling finer control over the final output values. After
#' dividing each score by the mean expression of stable genes, the result is multiplied by this
#' factor.
#' @param impute From [multiclassPairs::predict_RF()]. Boolean. To determine if missed genes and NA
#' values should be imputed or not. The non missed rules will be used to determine the closest
#' samples in the training binary matrix (i.e. which is stored in the classifier object). For each
#' sample, the mode value for nearest samples in the training data will be assigned to the missed
#' rules. Default is TRUE
#' @param impute_kNN From [multiclassPairs::predict_RF()]. Integer determines the number of the
#' nearest samples in the training data to be used in the imputation. Default is 5. It is not
#' recommended to use large number (i.e. >10).
#' @param impute_reject From [multiclassPairs::predict_RF()]. A number between 0 and 1 indicating
#' the threshold of the missed rules in the sample. Based on this threshold the sample will be
#' rejected (i.e. skipped if higher than the impute_reject threshold) and the missed rules will not
#' be imputed in this sample. Default is 0.67. NOTE, The results object will not have any results
#' for this sample.
#' @param verbose A logical value indicating whether processing messages will be printed or not.
#' Default is TRUE.
#' @param subtype_only Boolean parameter. Set to TRUE to return subtypes and nothing else. Default
#' is FALSE.
#' @param include_data Boolean parameter. Set to TRUE to include data in output, FALSE is default.
#' @param include_pred_scores Boolean parameter. Set to TRUE (default) to include prediction scores
#' for each sample and class in output.
#'
#' @return Returns a list object including: Data (optional, not included by default), Prediction
#' scores for all classes (optional, included by default), Predicted LundTax class for 7-class
#' system, Predicted LundTax class for 5-class system, as well a data frame with missing genes
#' information.
#'
#' @import dplyr
#' @importFrom stats setNames median na.omit quantile
#' @importFrom multiclassPairs predict_RF
#'
#' @export
#'
#' @examples
#' #load your gene expression data, genes as rows, samples as columns
#' data("sjodahl_2017")
#'
#' #classify samples
#' sjodahl_classes = classify_samples(this_data = sjodahl_2017,
#'                                    log_transform = FALSE,
#'                                    adjust = TRUE,
#'                                    impute = TRUE,
#'                                    include_data = TRUE,
#'                                    verbose = FALSE)
#'
#' #view data
#' summary(sjodahl_classes)
#'
classify_samples = function(this_data = NULL,
                            gene_id = "hgnc_symbol",
                            threshold_progression = 0.58,
                            threshold_grade = 0.5,
                            log_transform = FALSE,
                            adjust = TRUE,
                            adj_factor = 5.1431,
                            impute = TRUE,
                            impute_reject = 0.67,
                            impute_kNN = 5,
                            subtype_only = FALSE,
                            include_data = FALSE,
                            include_pred_scores = TRUE,
                            verbose = TRUE){

  #checks
  if(is.null(this_data)){
    stop("Input data is missing...")
  }

  if(verbose) {
    cat("\n", rep("=", 60), sep = "")
    cat("\n  LUND TAXONOMY CLASSIFIER")
    cat("\n", rep("=", 60), sep = "", "\n")
  }

  #check if data is log2 transformed
  log_check = int_check_log2_transformation(expression_df = this_data)

  if(verbose){
    cat("  Data Check:", log_check$message, "\n")
    cat("  Input samples:", ncol(this_data), "\n")
    cat("  Input genes:", nrow(this_data), "\n")
    cat("  Gene ID format:", gene_id)
  }

  if(log_transform){
    if(log_check$log2_transformed){
      message("WARNING: Data appears to already be log2 transformed.\nPlease ensure data type and evaluate `log_transform` argument accordingly...")
    }
  }

  #check the incoming data
  if(!class(this_data)[1] %in% c("data.frame","matrix")){
    stop("Data must be in dataframe or matrix format...")
  }

  #if input is a matrix, convert to data.frame
  if(is.matrix(this_data)){
    this_data <- as.data.frame(this_data)
  }

  #check for duplicated samples
  if(ncol(this_data) != length(unique(colnames(this_data)))){
    stop("Sample names (column names) should not be duplicated")
  }

  #check gene format
  if(!gene_id %in% c("hgnc_symbol", "ensembl_gene_id")){
    stop("gene_id must be one of the following: 'hgnc_symbol' or 'ensembl_gene_id'")
  }

  #stop if ensembl_gene_id is selected
  if(gene_id == "ensembl_gene_id"){
    if(verbose) cat("  Converting: Ensembl IDs to HGNC symbols\n")

    if(!exists("tx2gene")) {
      stop("Mapping table 'tx2gene' not found. Please provide a data.frame with columns 'gene_id' and 'gene_name'.")
    }

    hgnc_mapped <- tx2gene$gene_name[match(rownames(this_data), tx2gene$gene_id)]
    valid <- !is.na(hgnc_mapped) & hgnc_mapped != ""
    this_data <- this_data[valid, , drop = FALSE]
    rownames(this_data) <- hgnc_mapped[valid]
    gene_id <- "hgnc_symbol"

    if(verbose) {
      cat("  Genes mapped:", sum(valid), "\n")
      cat("  Genes unmapped:", sum(!valid), "\n")
    }
  }

  #log transform
  if(log_transform){
    if(verbose) cat("  Applying: log2 transformation\n")
    this_data = log2(this_data + 1)
  }

  if(verbose) {
    cat("\n", rep("=", 60), sep = "")
    cat("\n  SUBTYPE CLASSIFICATION")
    cat("\n", rep("=", 60), sep = "", "\n")
  }

  #initiate results object
  results_suburo <- list(data = this_data,
                         subtype_scores = NULL,
                         predictions_7classes = NULL,
                         predictions_5classes = NULL,
                         scores = NULL,
                         na_genes = NULL)

  # Step 1: 5-class prediction
  if(verbose) {
    cat("  1. Predicting 5-class subtypes")
  }

  prediction = predict_RF(classifier = classifier_lundtax_5c,
                          Data = this_data,
                          impute = impute,
                          impute_reject = impute_reject,
                          impute_kNN = impute_kNN,
                          verbose = FALSE)

  #reorder scores
  pred_order = prediction$predictions[,c("Uro","GU","BaSq","Mes","ScNE"), drop = FALSE]
  prediction$predictions = pred_order
  prediction$predictions_classes = colnames(pred_order)[max.col(replace(pred_order, is.na(pred_order), -Inf), ties.method = "first")]

  subtype_counts_5c <- table(factor(prediction$predictions_classes, levels = c("Uro", "GU", "BaSq", "Mes", "ScNE")))

  # Step 2: Uro subclassification
  if("Uro" %in% prediction$predictions_classes){

    n_uro <- sum(prediction$predictions_classes == "Uro")

    if(verbose) {
      cat("\n", rep("-", 60), "\n", sep = "")
      cat("  2. Subclassifying Uro samples")
    }

    pred_uro = this_data[,which(prediction$predictions_classes == "Uro"), drop = FALSE]
    pred_nouro = this_data[,which(prediction$predictions_classes != "Uro"), drop = FALSE]

    prediction_suburo = predict_RF(classifier = classifier_lundtax_7c,
                                   Data = pred_uro,
                                   impute = impute,
                                   impute_reject = impute_reject,
                                   impute_kNN = impute_kNN,
                                   verbose = FALSE)

    suburo_counts <- table(factor(colnames(prediction_suburo$predictions)[max.col(prediction_suburo$predictions)],
                                  levels = c("UroA", "UroB", "UroC")))

    names_uro = colnames(pred_uro)
    names_all = colnames(this_data)

    score_matrix = int_merge_suburo_matrix(score_matrix1 = prediction_suburo$predictions,
                                           score_matrix2 = prediction$predictions,
                                           row.names = list(names_uro,names_all))
  } else {
    score_matrix <- cbind("Uro" = prediction$predictions[,"Uro"],
                          "UroA" = NA,
                          "UroB" = NA,
                          "UroC" = NA,
                          "GU" = prediction$predictions[,"GU"],
                          "BaSq" = prediction$predictions[,"BaSq"],
                          "Mes" = prediction$predictions[,"Mes"],
                          "ScNE" = prediction$predictions[,"ScNE"])
  }

  # Step 3: Additional signatures
  if(verbose){
    cat("\n", rep("=", 60), sep = "")
    cat("\n  SIGNATURE CALCULATION")
    cat("\n", rep("=", 60), sep = "", "\n")
  }

  if(!subtype_only){
    all_scores = int_calc_signatures(this_data = this_data,
                                     gene_id = gene_id,
                                     threshold_progression = threshold_progression,
                                     threshold_grade = threshold_grade,
                                     adjust = adjust,
                                     adj_factor = adj_factor,
                                     verbose = verbose)
  } else {
    all_scores = NULL
  }

  #gather return
  results_suburo$scores = all_scores$merged_scores
  results_suburo$na_genes = all_scores$na_genes

  #calculate prediction_delta
  prediction_delta = apply(score_matrix, 1, function(row) {
    valid_scores = sort(row[!is.na(row)], decreasing = TRUE)
    if(length(valid_scores) >= 2){
      return(valid_scores[1] - valid_scores[2])
    } else if(length(valid_scores) == 1){
      return(valid_scores[1])
    } else {
      return(NA)
    }
  })

  score_matrix = cbind(score_matrix, prediction_delta = prediction_delta)
  results_suburo$subtype_scores = score_matrix

  score_matrix_suburo = score_matrix[,2:4, drop = FALSE]
  score_matrix_5c = score_matrix[,c(1,5:8), drop = FALSE]

  #5 class level
  results_suburo$predictions_5classes = colnames(score_matrix_5c)[max.col(replace(score_matrix_5c,is.na(score_matrix_5c),-Inf), ties.method = "first")]

  #7 class level
  results_suburo$predictions_7classes = results_suburo$predictions_5classes
  max_suburo = colnames(score_matrix_suburo)[max.col(replace(score_matrix_suburo, is.na(score_matrix_suburo), -Inf), ties.method = "first")]

  for(i in 1:length(results_suburo$predictions_7classes)){
    p = results_suburo$predictions_7classes[i]
    if(p == "Uro"){
      suburo = max_suburo[i]
      results_suburo$predictions_7classes[i] = suburo
    }
  }

  #score ties
  first5 = setNames(colnames(score_matrix_5c)[max.col(replace(score_matrix_5c, is.na(score_matrix_5c), -Inf), ties.method = "first")], rownames(score_matrix_5c))
  last5  = setNames(colnames(score_matrix_5c)[max.col(replace(score_matrix_5c, is.na(score_matrix_5c), -Inf), ties.method = "last")], rownames(score_matrix_5c))

  if(sum(first5 != last5) > 0){
    int_check_ties(first5,last5)
  }

  #7 class level
  score_matrix_suburo_ties = score_matrix_suburo[!is.na(score_matrix_suburo[,1]),]
  first7 = setNames(colnames(score_matrix_suburo_ties)[max.col(score_matrix_suburo_ties, ties.method = "first")], rownames(score_matrix_suburo_ties))
  last7 = setNames(colnames(score_matrix_suburo_ties)[max.col(score_matrix_suburo_ties, ties.method = "last")], rownames(score_matrix_suburo_ties))

  if(sum(first7 != last7) > 0){
    int_check_ties(first7,last7)
  }

  #final results
  names(results_suburo$predictions_7classes) = colnames(this_data)
  names(results_suburo$predictions_5classes) = colnames(this_data)

  if(subtype_only){
    predictions_suburo = list(predictions_7classes = results_suburo$predictions_7classes,
                              predictions_5classes = results_suburo$predictions_5classes)

    results_suburo_nodata = list(subtype_scores = results_suburo$subtype_scores,
                                 predictions_7classes = results_suburo$predictions_7classes,
                                 predictions_5classes = results_suburo$predictions_5classes)
  } else {
    predictions_suburo = list(predictions_7classes = results_suburo$predictions_7classes,
                              predictions_5classes = results_suburo$predictions_5classes,
                              scores = results_suburo$scores,
                              na_genes = results_suburo$na_genes)

    results_suburo_nodata = list(subtype_scores = results_suburo$subtype_scores,
                                 predictions_7classes = results_suburo$predictions_7classes,
                                 predictions_5classes = results_suburo$predictions_5classes,
                                 scores = results_suburo$scores,
                                 na_genes = all_scores$na_genes)
  }

  if(include_data & include_pred_scores){
    result = results_suburo
  } else if(include_data == FALSE & include_pred_scores){
    result = results_suburo_nodata
  } else if(include_data == FALSE & include_pred_scores == FALSE){
    result = predictions_suburo
  }

  if(verbose){
    cat(rep("=", 60), "\n", sep = "")
    cat("  CLASSIFICATION COMPLETE\n")
    cat(rep("=", 60), "\n", sep = "")

    # Print final summary in desired order: UroA, UroB, UroC, then 5-class
    for(st in c("UroA", "UroB", "UroC")) {
      if(suburo_counts[st] > 0) {
        cat("    ", st, ": ", suburo_counts[st], " samples\n", sep = "")
      }
    }

    for(st in c("GU", "BaSq", "Mes", "ScNE")) {
      if(subtype_counts_5c[st] > 0) {
        cat("    ", st, ": ", subtype_counts_5c[st], " samples\n", sep = "")
      }
    }
    cat(rep("-", 60), "\n", sep = "")
  }
  return(result)
}
