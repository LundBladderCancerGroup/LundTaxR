#' @title Predict Tumor Grade
#'
#' @description Apply SwitchBox Classifier to Gene Expression Data
#'
#' @details Internal function called by [LundTaxR::int_calc_signatures()].
#' Not meant for out of package use. This function applies a pre-trained SwitchBox classifier to
#' gene expression data to predict molecular tumor grades. The classifier uses Top Scoring Pairs
#' (TSPs) methodology, comparing expression levels between gene pairs to generate predictions.
#' Each rule evaluates whether gene1 > gene2, and satisfied rules contribute their
#' weights to the final score.
#'
#' @param data Gene expression matrix (genes × samples). Row names must be gene identifiers.
#' @param classifier SwitchBox classifier object with $TSPs, $score components
#' @param grade_threshold Numeric cutoff for binary classification (default 0.5)
#' @param grade_labels Character vector of length 2: c("low_label", "high_label")
#' @param verbose Print diagnostic information
#'
#' @return Data frame with sample IDs as row names and two columns:
#'   \describe{
#'     \item{prediction_score}{Numeric scores (0-1 scale)}
#'     \item{predicted_class}{Character classification based on threshold}
#'   }
#'
#' @examples
#' \dontrun{
#' # No examples provided
#' }
#'
int_predict_grade <- function(data,
                              classifier,
                              grade_threshold = 0.5,
                              grade_labels = NULL,
                              verbose = TRUE){

  #validate inputs
  if(!all(c("TSPs", "score") %in% names(classifier))){
    stop("Classifier must have $TSPs and $score components")
  }

  if(is.null(grade_labels) || length(grade_labels) != 2){
    stop("grade_labels must be a character vector of length 2: c('low_label', 'high_label')")
  }

  #extract gene pairs from classifier
  gene1_ids <- classifier$TSPs[, "gene1"]
  gene2_ids <- classifier$TSPs[, "gene2"]

  #find which rules have both genes present in data
  gene1_present <- gene1_ids %in% rownames(data)
  gene2_present <- gene2_ids %in% rownames(data)
  both_present <- gene1_present & gene2_present

  #filter to valid rules only
  valid_gene1 <- gene1_ids[both_present]
  valid_gene2 <- gene2_ids[both_present]
  valid_scores <- classifier$score[both_present]

  #print diagnostic info
  if(verbose){
    cat("\n         Total rules in classifier:", length(gene1_ids))
    cat("\n         Rules with both genes present:", sum(both_present))
    cat("\n         Rules removed:", sum(!both_present))
    cat("\n         Samples to score:", ncol(data))
    cat("\n         Grade threshold:", grade_threshold)
    cat("\n         Grade labels:", paste(grade_labels, collapse = " / "))
    cat("\n", rep("-", 60), sep = "")
  }

  #check if we have enough rules
  if(sum(both_present) == 0){
    stop("No valid rules found! Check that gene IDs match between data and classifier.")
  }

  if(sum(both_present) < 10){
    warning("Very few rules remaining (", sum(both_present), "). Results may be unreliable.")
  }

  #extract expression values for valid gene pairs
  gene1_expr <- data[valid_gene1, , drop = FALSE]
  gene2_expr <- data[valid_gene2, , drop = FALSE]

  #check which rules are satisfied (gene1 > gene2)
  rules_satisfied <- gene1_expr > gene2_expr

  #calculate weighted score for each sample
  #score = (sum of weights for satisfied rules) / (total weight)
  weighted_rules <- (rules_satisfied * 1) * valid_scores
  scores <- colSums(weighted_rules) / sum(valid_scores)

  #assign class based on threshold
  predicted_class <- ifelse(scores >= grade_threshold,
                            grade_labels[2],
                            grade_labels[1])

  #create result data frame
  result <- data.frame(
    prediction_score = scores,
    predicted_class = predicted_class,
    row.names = names(scores)
  )

  return(result)
}
