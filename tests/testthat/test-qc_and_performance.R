library(testthat)
library(LundTaxR)

####################################################################################################
#fast gene ranking gives the same result as rank()
test_that("int_rank_genes matches apply(rank) / nrow, with and without ties and NAs", {

  data("sjodahl_2017")

  #rounding creates many ties
  expr = round(as.matrix(sjodahl_2017[1:2000, 1:10]), 1)
  genes = rownames(expr)[seq(1, 2000, by = 7)]
  expected = (apply(expr, 2, rank) / nrow(expr))[genes, ]

  expect_equal(LundTaxR:::int_rank_genes(this_data = expr, genes = genes), expected, tolerance = 0)

  #NA values use the rank() fallback
  expr[5, 3] = NA
  expected_na = apply(expr, 2, rank)
  rownames(expected_na) = rownames(expr)
  expected_na = (expected_na / nrow(expr))[genes, ]

  expect_equal(LundTaxR:::int_rank_genes(this_data = expr, genes = genes), expected_na, tolerance = 0)
})

####################################################################################################
#data frame and matrix input give the same result
test_that("classify_samples gives the same result for data frame and matrix input", {

  data("sjodahl_2017")

  expr = sjodahl_2017[, 1:20]
  res_df = classify_samples(this_data = expr, verbose = FALSE)
  res_mat = classify_samples(this_data = as.matrix(expr), verbose = FALSE)

  expect_equal(res_df$predictions_7classes, res_mat$predictions_7classes)
  expect_equal(res_df$subtype_scores, res_mat$subtype_scores)
  expect_equal(res_df$scores, res_mat$scores)
})

####################################################################################################
#samples rejected by predict_RF are kept as NA
test_that("samples with too many missing genes are kept with NA predictions", {

  data("sjodahl_2017")

  expr = sjodahl_2017[, 1:8]
  classifier_genes = LundTaxR::classifier_lundtax_5c$RF_scheme$genes
  expr[rownames(expr) %in% classifier_genes[1:650], 1] = NA

  res = suppressMessages(classify_samples(this_data = expr, verbose = FALSE))

  expect_length(res$predictions_5classes, 8)
  expect_equal(names(res$predictions_5classes), colnames(expr))
  expect_true(is.na(res$predictions_5classes[1]))
  expect_true(is.na(res$predictions_7classes[1]))
  expect_true(all(is.na(res$subtype_scores[1, c("Uro", "GU", "BaSq", "Mes", "ScNE")])))
  expect_false(anyNA(res$predictions_5classes[-1]))
  expect_message(classify_samples(this_data = expr, verbose = FALSE), "could not be classified")
})

####################################################################################################
#qc report
test_that("get_qc_report returns the expected components", {

  data("sjodahl_2017")

  res = classify_samples(this_data = sjodahl_2017[, 1:30], include_data = TRUE, verbose = FALSE)
  qc = get_qc_report(these_predictions = res)

  expect_s3_class(qc, "lundtax_qc")
  expect_named(qc, c("summary", "samples", "signature_coverage", "classifier_coverage", "data_check"))
  expect_equal(nrow(qc$samples), 30)
  expect_equal(qc$samples$sample_id, colnames(sjodahl_2017)[1:30])
  expect_equal(qc$summary$n_samples, 30)
  expect_equal(nrow(qc$signature_coverage), 4)
  expect_equal(nrow(qc$classifier_coverage), 4)
  expect_true(all(qc$samples$missing_rules_5_class == 0))
  expect_output(print(qc), "LUNDTAX QC REPORT")

  #low confidence flag follows the threshold
  qc_strict = get_qc_report(these_predictions = res, delta_threshold = 1)
  expect_equal(qc_strict$summary$n_low_confidence, sum(!is.na(qc_strict$samples$prediction_delta_collapsed) &
                                                         qc_strict$samples$prediction_delta_collapsed < 1))
})

test_that("get_qc_report works without expression data and flags unclassified samples", {

  data("sjodahl_2017")

  expr = sjodahl_2017[, 1:8]
  classifier_genes = LundTaxR::classifier_lundtax_5c$RF_scheme$genes
  expr[rownames(expr) %in% classifier_genes[1:650], 1] = NA

  res = suppressMessages(classify_samples(this_data = expr, verbose = FALSE))
  qc = get_qc_report(these_predictions = res)

  expect_null(qc$classifier_coverage)
  expect_null(qc$data_check)
  expect_equal(qc$summary$n_not_classified, 1)
  expect_false(qc$samples$classified[1])
  expect_false(qc$samples$low_confidence[1])

  #with the data supplied, the rejected sample has all its classifier rules missing
  qc_data = get_qc_report(these_predictions = res, this_data = expr)
  expect_equal(qc_data$samples$missing_rules_5_class[1], 1)
})

test_that("get_qc_report validates its input", {
  expect_error(get_qc_report(), "classify_samples")
  expect_error(get_qc_report(these_predictions = list(predictions_5classes = c(a = "Uro"))), "include_pred_scores")
})

####################################################################################################
#expression distribution
test_that("get_qc_report flags samples with a shifted expression distribution", {

  data("sjodahl_2017")

  expr = sjodahl_2017[, 1:40]
  res_clean = classify_samples(this_data = expr, include_data = TRUE, verbose = FALSE)
  qc_clean = get_qc_report(these_predictions = res_clean)

  expect_true(all(c("expr_median", "expr_q25", "expr_q75", "expr_iqr", "pct_zero", "distribution_outlier") %in% colnames(qc_clean$samples)))
  expect_equal(qc_clean$summary$n_distribution_outliers, 0)
  expect_equal(qc_clean$samples$expr_median[1], median(expr[, 1]))

  #shift one sample and set part of another to zero
  expr[, 5] = expr[, 5] - 2.5
  expr[1:4000, 6] = 0
  res_shift = classify_samples(this_data = expr, include_data = TRUE, verbose = FALSE)
  qc_shift = get_qc_report(these_predictions = res_shift)

  expect_equal(qc_shift$samples$sample_id[qc_shift$samples$distribution_outlier], colnames(expr)[c(5, 6)])
  expect_gt(qc_shift$samples$pct_zero[6], 20)
})

####################################################################################################
#qc plot
test_that("plot_qc_report returns a plot with and without expression data", {

  data("sjodahl_2017")

  res = classify_samples(this_data = sjodahl_2017[, 1:30], include_data = TRUE, verbose = FALSE)

  p_full = plot_qc_report(these_qc = get_qc_report(these_predictions = res))
  expect_s3_class(p_full, "patchwork")

  res$data = NULL
  p_nodata = plot_qc_report(these_predictions = res, subtype_class = "7_class")
  expect_s3_class(p_nodata, "patchwork")

  #large cohorts
  qc_big = get_qc_report(these_predictions = classify_samples(this_data = sjodahl_2017, include_data = TRUE, verbose = FALSE))
  qc_big$samples = qc_big$samples[rep(seq_len(nrow(qc_big$samples)), 3), ]
  expect_s3_class(plot_qc_report(these_qc = qc_big), "patchwork")

  expect_error(plot_qc_report(), "get_qc_report|classify_samples")
  expect_error(plot_qc_report(these_qc = list()), "get_qc_report")
})
