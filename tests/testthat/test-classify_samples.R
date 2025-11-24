library(testthat)
library(LundTaxR)
library(tibble)
library(dplyr)

####################################################################################################
#classifier overview
test_that("classify_samples returns a list with expected components", {

  #load test data
  data("sjodahl_2017")

  #run the classifier
  result <- classify_samples(this_data = sjodahl_2017, log_transform = FALSE)
  expect_type(result, "list")
  expect_true("predictions_5classes" %in% names(result))
  expect_true("predictions_7classes" %in% names(result))
})

####################################################################################################
#scores object
test_that("classify_samples returns correct scores for first row", {

  #load test data
  data("sjodahl_2017")

  #run the classifier
  result <- classify_samples(
    this_data = sjodahl_2017,
    log_transform = FALSE,
    adjust = TRUE,
    impute = TRUE,
    include_data = TRUE,
    verbose = FALSE
  )

  #subset the scores to the first sample and add sample ID as first column
  scores_results <- as.data.frame(result$scores[1,]) %>%
    rownames_to_column("sample_id")

  #get the expected scores for the first sample
  expected <- data.frame(sample_id = "1.CEL",
    proliferation_score = 0.7310869,
    progression_score = 0.7675475,
    progression_risk = "HR",
    molecular_grade_who_2022_score = 0.7736667,
    molecular_grade_who_2022 = "HG",
    molecular_grade_who_1999_score = 0.8074286,
    molecular_grade_who_1999 = "G3",
    stromal141_up = 4.047838,
    immune141_up = 4.27489,
    b_cells = 3.211749,
    b_cells_proportion = 0.06186246,
    t_cells = 3.569395,
    t_cells_proportion = 0.06875119,
    t_cells_cd8 = 1.966563,
    t_cells_cd8_proportion = 0.03787857,
    nk_cells = 3.162512,
    nk_cells_proportion = 0.06091409,
    cytotoxicity_score = 3.530327,
    cytotoxicity_score_proportion = 0.06799869,
    neutrophils = 3.405826,
    neutrophils_proportion = 0.06560065,
    monocytic_lineage = 4.281524,
    monocytic_lineage_proportion = 0.08246773,
    macrophages = 4.844223,
    macrophages_proportion = 0.09330604,
    m2_macrophage = 3.806195,
    m2_macrophage_proportion = 0.07331228,
    myeloid_dendritic_cells = 2.997839,
    myeloid_dendritic_cells_proportion = 0.05774229,
    endothelial_cells = 3.781402,
    endothelial_cells_proportion = 0.07283472,
    fibroblasts = 6.999977,
    fibroblasts_proportion = 0.1348287,
    smooth_muscle = 6.360039,
    smooth_muscle_proportion = 0.1225026
  )

  #make expected columns factors with the same levels as scores_results
  expected$progression_risk <- factor(expected$progression_risk, levels = levels(scores_results$progression_risk))
  expected$molecular_grade_who_2022 <- factor(expected$molecular_grade_who_2022, levels = levels(scores_results$molecular_grade_who_2022))
  expected$molecular_grade_who_1999 <- factor(expected$molecular_grade_who_1999, levels = levels(scores_results$molecular_grade_who_1999))

  #do we get back the exact same values for the first sample
  expect_equal(expected, scores_results, tolerance = 1e-6, check.attributes = FALSE)

  #check if column names are the expected
  expect_equal(colnames(expected), colnames(scores_results))

})

####################################################################################################
#5class subtypes
test_that("classify_samples returns correct 5-class subtypes", {
  #load test data
  data("sjodahl_2017")
  #run the classifier
  result <- classify_samples(
    this_data = sjodahl_2017,
    log_transform = FALSE,
    adjust = TRUE,
    impute = TRUE,
    include_data = TRUE,
    verbose = FALSE
  )
  #expected subtypes
  expected_subtypes <- data.frame(
    sample_id = c(
      "1.CEL", "2.CEL", "3.CEL", "4.CEL", "5.CEL", "6.CEL", "7.CEL", "8.CEL", "9.CEL", "10.CEL", "11.CEL", "12.CEL", "13.CEL",
      "14.CEL", "15.CEL", "16.CEL", "17.CEL", "18.CEL", "19.CEL", "20.CEL", "21.CEL", "22.CEL", "23.CEL", "24.CEL", "25.CEL", "26.CEL",
      "27.CEL", "28.CEL", "29.CEL", "30.CEL", "32.CEL", "33.CEL", "34.CEL", "35.CEL", "36.CEL", "37.CEL", "38.CEL", "39.CEL", "40.CEL",
      "41.CEL", "43.CEL", "44.CEL", "45.CEL", "46.CEL", "47.CEL", "48.CEL", "49.CEL", "50.CEL", "51.CEL", "52.CEL", "53.CEL", "54.CEL",
      "56.CEL", "57.CEL", "58.CEL", "59.CEL", "60.CEL", "61.CEL", "62.CEL", "63.CEL", "64.CEL", "65.CEL", "66.CEL", "67.CEL", "68.CEL",
      "69.CEL", "70.CEL", "71.CEL", "72.CEL", "73.CEL", "74.CEL", "75.CEL", "76.CEL", "77.CEL", "78.CEL", "79.CEL", "80.CEL", "81.CEL",
      "82.CEL", "83.CEL", "84.CEL", "85.CEL", "87.CEL", "88.CEL", "89.CEL", "90.CEL", "91.CEL", "92.CEL", "93.CEL", "94.CEL", "95.CEL",
      "96.CEL", "97.CEL", "98.CEL", "99.CEL", "100.CEL", "101.CEL", "102.CEL", "104.CEL", "105.CEL", "106.CEL", "107.CEL", "108.CEL", "109.CEL",
      "110.CEL", "111.CEL", "112.CEL", "113.CEL", "114.CEL", "115.CEL", "116.CEL", "117.CEL", "118.CEL", "119.CEL", "120.CEL", "121.CEL", "122.CEL",
      "123.CEL", "124.CEL", "125.CEL", "126.CEL", "127.CEL", "128.CEL", "129.CEL", "130.CEL", "132.CEL", "133.CEL", "135.CEL", "136.CEL", "137.CEL",
      "138.CEL", "139.CEL", "140.CEL", "141.CEL", "142.CEL", "143.CEL", "144.CEL", "145.CEL", "147.CEL", "148.CEL", "149.CEL", "150.CEL", "151.CEL",
      "153.CEL", "154.CEL", "155.CEL", "156.CEL", "157.CEL", "158.CEL", "159.CEL", "160.CEL", "161.CEL", "162.CEL", "163.CEL", "164.CEL", "165.CEL",
      "166.CEL", "167.CEL", "169.CEL", "170.CEL", "171.CEL", "172.CEL", "173.CEL", "174.CEL", "175.CEL", "176.CEL", "178.CEL", "179.CEL", "180.CEL",
      "181.CEL", "182.CEL", "183.CEL", "184.CEL", "185.CEL", "186.CEL", "187.CEL", "188.CEL", "189.CEL", "190.CEL", "191.CEL", "192.CEL", "193.CEL",
      "194.CEL", "195.CEL", "196.CEL", "197.CEL", "198.CEL", "199.CEL", "200.CEL", "201.CEL", "202.CEL", "205.CEL", "206.CEL", "207.CEL", "208.CEL",
      "211.CEL", "212.CEL", "213.CEL", "214.CEL", "215.CEL", "216.CEL", "217.CEL", "218.CEL", "219.CEL", "220.CEL", "221.CEL", "222.CEL", "223.CEL",
      "224.CEL", "225.CEL", "226.CEL", "227.CEL", "228.CEL", "229.CEL", "230.CEL", "231.CEL", "233.CEL", "235.CEL", "236.CEL", "237.CEL", "238.CEL",
      "239.CEL", "240.CEL", "241.CEL", "242.CEL", "243.CEL", "244.CEL", "245.CEL", "246.CEL", "247.CEL", "248.CEL", "249.CEL", "250.CEL", "252.CEL",
      "253.CEL", "254.CEL", "255.CEL", "257.CEL", "258.CEL", "259.CEL", "260.CEL", "261.CEL", "263.CEL", "264.CEL", "266.CEL", "268.CEL", "270.CEL",
      "271.CEL", "272.CEL", "273.CEL", "274.CEL", "278.CEL", "279.CEL", "281.CEL", "284.CEL", "288.CEL", "289.CEL", "290.CEL", "291.CEL", "292.CEL",
      "294.CEL", "295.CEL", "296.CEL", "299.CEL", "301.CEL", "304.CEL", "305.CEL"
    ),
    subtype = c(
      "Uro", "Mes", "GU", "Uro", "Uro", "GU", "Uro", "ScNE", "GU", "Uro", "GU", "Uro", "Uro",
      "GU", "Uro", "GU", "GU", "BaSq", "Mes", "Uro", "Uro", "Uro", "BaSq", "BaSq", "BaSq", "Uro",
      "BaSq", "Uro", "Uro", "GU", "Uro", "GU", "ScNE", "ScNE", "BaSq", "Uro", "Uro", "Mes", "ScNE",
      "BaSq", "Uro", "BaSq", "Uro", "GU", "GU", "BaSq", "Uro", "Uro", "Uro", "Uro", "Uro", "BaSq",
      "GU", "Uro", "Uro", "Uro", "Uro", "ScNE", "BaSq", "GU", "Uro", "Uro", "Mes", "Uro", "Uro",
      "GU", "Mes", "Uro", "Uro", "GU", "ScNE", "BaSq", "Uro", "GU", "GU", "Uro", "BaSq", "BaSq",
      "GU", "GU", "Uro", "ScNE", "Uro", "BaSq", "Uro", "Uro", "BaSq", "ScNE", "Uro", "BaSq", "Mes",
      "Uro", "GU", "Uro", "Uro", "GU", "Uro", "Uro", "Uro", "BaSq", "GU", "Uro", "Uro", "BaSq",
      "GU", "Uro", "Uro", "BaSq", "BaSq", "Mes", "GU", "Uro", "Uro", "Uro", "Uro", "Uro", "Uro",
      "BaSq", "Uro", "ScNE", "BaSq", "BaSq", "Uro", "GU", "GU", "Uro", "GU", "Mes", "Mes", "Uro",
      "Uro", "Uro", "Uro", "Uro", "BaSq", "BaSq", "Uro", "BaSq", "GU", "BaSq", "BaSq", "Uro", "GU",
      "BaSq", "GU", "ScNE", "BaSq", "BaSq", "Mes", "ScNE", "Uro", "GU", "Uro", "BaSq", "Uro", "BaSq",
      "ScNE", "Uro", "Uro", "Uro", "Uro", "GU", "Uro", "BaSq", "Uro", "Uro", "Uro", "Uro", "Uro",
      "BaSq", "Uro", "Uro", "BaSq", "BaSq", "Uro", "ScNE", "GU", "Uro", "Uro", "Mes", "Uro", "BaSq",
      "BaSq", "Uro", "GU", "ScNE", "Uro", "Uro", "ScNE", "Uro", "Uro", "GU", "GU", "Uro", "BaSq",
      "ScNE", "BaSq", "ScNE", "Mes", "Uro", "BaSq", "Uro", "BaSq", "Mes", "Uro", "GU", "BaSq", "ScNE",
      "BaSq", "Mes", "BaSq", "Mes", "BaSq", "Uro", "BaSq", "GU", "GU", "GU", "BaSq", "GU", "Uro",
      "Uro", "Uro", "Uro", "GU", "ScNE", "Uro", "BaSq", "Uro", "GU", "GU", "Uro", "Uro", "Uro",
      "Uro", "Uro", "Uro", "BaSq", "Uro", "Uro", "Mes", "GU", "GU", "GU", "Uro", "Uro", "GU",
      "BaSq", "Uro", "Uro", "GU", "Uro", "Uro", "Uro", "BaSq", "GU", "GU", "ScNE", "Uro", "GU",
      "GU", "BaSq", "Uro", "Uro", "BaSq", "Uro", "GU"
    ),
    stringsAsFactors = FALSE
  )

  #actual subtypes
  actual_subtypes <- as.data.frame(result$predictions_5classes) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::rename(subtype = `result$predictions_5classes`)

  #check that names and values match
  expected_subtypes$subtype <- as.character(expected_subtypes$subtype)
  actual_subtypes$subtype <- as.character(actual_subtypes$subtype)
  expected_subtypes <- expected_subtypes[order(expected_subtypes$sample_id), ]
  actual_subtypes <- actual_subtypes[order(actual_subtypes$sample_id), ]
  expect_equal(actual_subtypes, expected_subtypes)
})

####################################################################################################
#7class subtypes
test_that("classify_samples returns correct 7-class subtypes", {

  #load test data
  data("sjodahl_2017")

  #run the classifier
  result <- classify_samples(
    this_data = sjodahl_2017,
    log_transform = FALSE,
    adjust = TRUE,
    impute = TRUE,
    include_data = TRUE,
    verbose = FALSE
  )

  #expected subtypes
  expected_subtypes <- data.frame(
    sample_id = c(
      "1.CEL", "2.CEL", "3.CEL", "4.CEL", "5.CEL", "6.CEL", "7.CEL", "8.CEL", "9.CEL", "10.CEL", "11.CEL", "12.CEL", "13.CEL",
      "14.CEL", "15.CEL", "16.CEL", "17.CEL", "18.CEL", "19.CEL", "20.CEL", "21.CEL", "22.CEL", "23.CEL", "24.CEL", "25.CEL", "26.CEL",
      "27.CEL", "28.CEL", "29.CEL", "30.CEL", "32.CEL", "33.CEL", "34.CEL", "35.CEL", "36.CEL", "37.CEL", "38.CEL", "39.CEL", "40.CEL",
      "41.CEL", "43.CEL", "44.CEL", "45.CEL", "46.CEL", "47.CEL", "48.CEL", "49.CEL", "50.CEL", "51.CEL", "52.CEL", "53.CEL", "54.CEL",
      "56.CEL", "57.CEL", "58.CEL", "59.CEL", "60.CEL", "61.CEL", "62.CEL", "63.CEL", "64.CEL", "65.CEL", "66.CEL", "67.CEL", "68.CEL",
      "69.CEL", "70.CEL", "71.CEL", "72.CEL", "73.CEL", "74.CEL", "75.CEL", "76.CEL", "77.CEL", "78.CEL", "79.CEL", "80.CEL", "81.CEL",
      "82.CEL", "83.CEL", "84.CEL", "85.CEL", "87.CEL", "88.CEL", "89.CEL", "90.CEL", "91.CEL", "92.CEL", "93.CEL", "94.CEL", "95.CEL",
      "96.CEL", "97.CEL", "98.CEL", "99.CEL", "100.CEL", "101.CEL", "102.CEL", "104.CEL", "105.CEL", "106.CEL", "107.CEL", "108.CEL", "109.CEL",
      "110.CEL", "111.CEL", "112.CEL", "113.CEL", "114.CEL", "115.CEL", "116.CEL", "117.CEL", "118.CEL", "119.CEL", "120.CEL", "121.CEL", "122.CEL",
      "123.CEL", "124.CEL", "125.CEL", "126.CEL", "127.CEL", "128.CEL", "129.CEL", "130.CEL", "132.CEL", "133.CEL", "135.CEL", "136.CEL", "137.CEL",
      "138.CEL", "139.CEL", "140.CEL", "141.CEL", "142.CEL", "143.CEL", "144.CEL", "145.CEL", "147.CEL", "148.CEL", "149.CEL", "150.CEL", "151.CEL",
      "153.CEL", "154.CEL", "155.CEL", "156.CEL", "157.CEL", "158.CEL", "159.CEL", "160.CEL", "161.CEL", "162.CEL", "163.CEL", "164.CEL", "165.CEL",
      "166.CEL", "167.CEL", "169.CEL", "170.CEL", "171.CEL", "172.CEL", "173.CEL", "174.CEL", "175.CEL", "176.CEL", "178.CEL", "179.CEL", "180.CEL",
      "181.CEL", "182.CEL", "183.CEL", "184.CEL", "185.CEL", "186.CEL", "187.CEL", "188.CEL", "189.CEL", "190.CEL", "191.CEL", "192.CEL", "193.CEL",
      "194.CEL", "195.CEL", "196.CEL", "197.CEL", "198.CEL", "199.CEL", "200.CEL", "201.CEL", "202.CEL", "205.CEL", "206.CEL", "207.CEL", "208.CEL",
      "211.CEL", "212.CEL", "213.CEL", "214.CEL", "215.CEL", "216.CEL", "217.CEL", "218.CEL", "219.CEL", "220.CEL", "221.CEL", "222.CEL", "223.CEL",
      "224.CEL", "225.CEL", "226.CEL", "227.CEL", "228.CEL", "229.CEL", "230.CEL", "231.CEL", "233.CEL", "235.CEL", "236.CEL", "237.CEL", "238.CEL",
      "239.CEL", "240.CEL", "241.CEL", "242.CEL", "243.CEL", "244.CEL", "245.CEL", "246.CEL", "247.CEL", "248.CEL", "249.CEL", "250.CEL", "252.CEL",
      "253.CEL", "254.CEL", "255.CEL", "257.CEL", "258.CEL", "259.CEL", "260.CEL", "261.CEL", "263.CEL", "264.CEL", "266.CEL", "268.CEL", "270.CEL",
      "271.CEL", "272.CEL", "273.CEL", "274.CEL", "278.CEL", "279.CEL", "281.CEL", "284.CEL", "288.CEL", "289.CEL", "290.CEL", "291.CEL", "292.CEL",
      "294.CEL", "295.CEL", "296.CEL", "299.CEL", "301.CEL", "304.CEL", "305.CEL"
    ),
    subtype = c(
      "UroA", "Mes", "GU", "UroA", "UroB", "GU", "UroB", "ScNE", "GU", "UroC", "GU", "UroC", "UroB",
      "GU", "UroA", "GU", "GU", "BaSq", "Mes", "UroC", "UroB", "UroC", "BaSq", "BaSq", "BaSq", "UroC",
      "BaSq", "UroA", "UroB", "GU", "UroB", "GU", "ScNE", "ScNE", "BaSq", "UroA", "UroC", "Mes", "ScNE",
      "BaSq", "UroB", "BaSq", "UroA", "GU", "GU", "BaSq", "UroA", "UroC", "UroA", "UroA", "UroB", "BaSq",
      "GU", "UroC", "UroA", "UroC", "UroC", "ScNE", "BaSq", "GU", "UroC", "UroC", "Mes", "UroC", "UroC",
      "GU", "Mes", "UroA", "UroC", "GU", "ScNE", "BaSq", "UroB", "GU", "GU", "UroA", "BaSq", "BaSq",
      "GU", "GU", "UroC", "ScNE", "UroA", "BaSq", "UroA", "UroC", "BaSq", "ScNE", "UroC", "BaSq", "Mes",
      "UroC", "GU", "UroB", "UroC", "GU", "UroA", "UroC", "UroB", "BaSq", "GU", "UroC", "UroC", "BaSq",
      "GU", "UroC", "UroB", "BaSq", "BaSq", "Mes", "GU", "UroA", "UroB", "UroB", "UroC", "UroB", "UroA",
      "BaSq", "UroC", "ScNE", "BaSq", "BaSq", "UroC", "GU", "GU", "UroA", "GU", "Mes", "Mes", "UroA",
      "UroB", "UroA", "UroA", "UroC", "BaSq", "BaSq", "UroA", "BaSq", "GU", "BaSq", "BaSq", "UroA", "GU",
      "BaSq", "GU", "ScNE", "BaSq", "BaSq", "Mes", "ScNE", "UroC", "GU", "UroC", "BaSq", "UroB", "BaSq",
      "ScNE", "UroA", "UroB", "UroA", "UroA", "GU", "UroB", "BaSq", "UroC", "UroC", "UroB", "UroC", "UroC",
      "BaSq", "UroA", "UroC", "BaSq", "BaSq", "UroC", "ScNE", "GU", "UroC", "UroC", "Mes", "UroA", "BaSq",
      "BaSq", "UroC", "GU", "ScNE", "UroB", "UroA", "ScNE", "UroC", "UroB", "GU", "GU", "UroA", "BaSq",
      "ScNE", "BaSq", "ScNE", "Mes", "UroA", "BaSq", "UroA", "BaSq", "Mes", "UroC", "GU", "BaSq", "ScNE",
      "BaSq", "Mes", "BaSq", "Mes", "BaSq", "UroB", "BaSq", "GU", "GU", "GU", "BaSq", "GU", "UroC",
      "UroC", "UroC", "UroA", "GU", "ScNE", "UroB", "BaSq", "UroB", "GU", "GU", "UroA", "UroC", "UroC",
      "UroA", "UroA", "UroC", "BaSq", "UroA", "UroA", "Mes", "GU", "GU", "GU", "UroA", "UroC", "GU",
      "BaSq", "UroC", "UroC", "GU", "UroB", "UroA", "UroA", "BaSq", "GU", "GU", "ScNE", "UroA", "GU",
      "GU", "BaSq", "UroC", "UroA", "BaSq", "UroB", "GU"
    ),
    stringsAsFactors = FALSE
  )

  #actual subtypes
  actual_subtypes <- as.data.frame(result$predictions_7classes) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::rename(subtype = `result$predictions_7classes`)

  #check that names and values match
  expected_subtypes$subtype <- as.character(expected_subtypes$subtype)
  actual_subtypes$subtype <- as.character(actual_subtypes$subtype)
  expected_subtypes <- expected_subtypes[order(expected_subtypes$sample_id), ]
  actual_subtypes <- actual_subtypes[order(actual_subtypes$sample_id), ]
  expect_equal(actual_subtypes, expected_subtypes)
})
