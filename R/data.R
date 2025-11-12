#' LundTax2023
#'
#' This packages implements a Random Forest rule-based single-sample predictor that classifies
#' transcriptomic samples into the 5 (or 7, including subclasses) Lund Taxonomy molecular subtypes.
#' The final classifier is composed of two separate predictors applied sequentially: first a sample
#' is classified as one of the 5 main classes (Uro, GU, BaSq, Mes or ScNE), and then samples classified as
#' Uro are subclassified into UroA, UroB or UroC by a second predictor
"_PACKAGE"

#' SwitchBox Classifier for WHO 1999 Tumor Grading
#'
#' A pre-trained SwitchBox classifier using Top Scoring Pairs (TSPs) methodology
#' to predict bladder tumor grades according to the WHO 1999 classification system.
#' The classifier distinguishes between low-grade (G1/G2) and high-grade (G3) tumors
#' based on gene expression patterns.
#'
#' @format A list with 5 components:
#' \describe{
#'   \item{name}{Character string; classifier identifier (e.g., "283TSPs")}
#'   \item{TSPs}{Character matrix with 2 columns ("gene1", "gene2") containing
#'     gene pair identifiers for each comparison rule. Gene identifiers are symbols
#'     (e.g., "TP53", "BRCA1")}
#'   \item{score}{Named numeric vector of rule weights/importance scores, one per
#'     gene pair. Higher values indicate more discriminative rules}
#'   \item{tieVote}{Factor indicating voting behavior for tied comparisons}
#'   \item{labels}{Character vector of length 2: class labels c("G1G2", "G3")}
#' }
#'
#' @details
#' This classifier was trained on bladder cancer gene expression data to predict
#' WHO 1999 histological grades. It uses rank-based gene pair comparisons that
#' are robust to technical batch effects and normalization differences.
#' 
#' **Usage:**
#' - Scores < 0.5 indicate low-grade tumors (G1/G2)
#' - Scores ≥ 0.5 indicate high-grade tumors (G3)
#' 
#' **Gene identifiers:** Gene symbols (HUGO nomenclature)
#' 
#' **Number of rules:** Varies (check \code{nrow(classifier_grade_who_1999$TSPs)})
#'
#' @seealso
#' \code{\link{classifier_grade_who_2022}} for WHO 2022 grading,
#' \code{\link{int_predict_grade}} to apply the classifier
#'
"classifier_grade_who_1999"

#' SwitchBox Classifier for WHO 2022 Tumor Grading
#'
#' A pre-trained SwitchBox classifier using Top Scoring Pairs (TSPs) methodology
#' to predict bladder tumor grades according to the WHO 2022 classification system.
#' The classifier distinguishes between low-grade (LG) and high-grade (HG) tumors
#' incorporating molecular and morphological features.
#'
#' @format A list with 5 components:
#' \describe{
#'   \item{name}{Character string; classifier identifier (e.g., "140TSPs")}
#'   \item{TSPs}{Character matrix with 2 columns ("gene1", "gene2") containing
#'     gene pair identifiers for each comparison rule. Gene identifiers are symbols
#'     (e.g., "TP53", "BRCA1")}
#'   \item{score}{Named numeric vector of rule weights/importance scores, one per
#'     gene pair. Higher values indicate more discriminative rules}
#'   \item{tieVote}{Factor indicating voting behavior for tied comparisons}
#'   \item{labels}{Character vector of length 2: class labels c("LG", "HG")}
#' }
#'
#' @details
#' This classifier was trained on bladder cancer gene expression data to predict
#' WHO 2022 molecular grades, which integrate both molecular alterations and
#' morphological features. The rank-based approach is robust to technical
#' variation between datasets.
#' 
#' **Usage:**
#' - Scores < 0.5 indicate low-grade tumors (LG)
#' - Scores ≥ 0.5 indicate high-grade tumors (HG)
#' 
#' **Gene identifiers:** Gene symbols (HUGO nomenclature)
#' 
#' **Number of rules:** Varies (check \code{nrow(classifier_grade_who_2022$TSPs)})
#'
#' @seealso
#' \code{\link{classifier_grade_who_1999}} for WHO 1999 grading,
#' \code{\link{int_predict_grade}} to apply the classifier
#'
"classifier_grade_who_2022"

#' Classifier LundTax 5c.
#'
#' Classifier as a 'rule_based_RandomForest' object. Predicts samples as one of the 5 main Lund 
#' Taxonomy molecular subtypes, Uro, GU, BaSq, Mes, or ScNE. Object includes the final RF classifier,
#' the used genes and rules in the final model, the Boruta results, and the training matrix. The 
#' training matrix is a binary matrix containing the rule values for the training data and it is
#' used for imputation purposes during the prediction if values are missing in the sample. This object
#' was generated using the [multiclassPairs::predict_RF()] function.
#'
#' A Large rule_based_RandomForest object. A list of 6.
#'
#' \itemize{
#'  \item genes. Genes in hgnc format.
#'  \item rules. A set of rules for the classifier in hgnc format.
#'  \item TrainingMatrix. Binary matrix for the rules in the training data.
#'  \item boruta. Boruta results for the classifier.
#'  \item RF_classifier. Random forest classifier details.
#'  \item calls. Information on how the model was generated.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name classifier_lundtax_5c
#' @usage data(classifier_lundtax_5c)
#' @format A list of 6.
NULL

#' Classifier LundTax 7c.
#'
#' Classifier as a rule_based_RandomForest object. Predicts samples as one of the 3 Uro subclasses, 
#' UroA, UroB, or UroC. Object includes the final RF classifier, the used genes and rules in the 
#' final model, the Boruta results, and the training matrix. The training matrix is a binary matrix
#' containing the rule values for the training data and it is used for imputation purposes during 
#' the prediction if values are missing in the sample. This object was generated using the
#' [multiclassPairs::predict_RF()] function.
#'
#' A Large rule_based_RandomForest object. A list of 6.
#'
#' \itemize{
#'  \item genes. Genes in hgnc format.
#'  \item rules. A set of rules for the classifier in hgnc format.
#'  \item TrainingMatrix. Binary matrix for the rules in the training data.
#'  \item boruta. Boruta results for the classifier.
#'  \item RF_classifier. Random forest classifier details.
#'  \item calls. Information on how the model was generated.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name classifier_lundtax_7c
#' @usage data(classifier_lundtax_7c)
#' @format A list of 6.
NULL

#' Gene List.
#'
#' Gene annotations in hgnc and ensembl format. Features, for the classifier relevant genes. 
#' For convenience, both formats are available for seamless conversion of IDs, regardless of the 
#' format of the incoming data.
#'
#' A data frame with gene information in different formats.
#'
#' \itemize{
#'  \item ensembl_gene_id. Gene annotation in ensembl format.
#'  \item hgnc_symbol. Gene annotation in HUGO format.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name gene_list
#' @usage data(gene_list)
#' @format A data frame with 1900 rows (genes) and 2 columns (hgnc_symbol and ensembl_gene_id).
NULL

#' Lund Colors.
#'
#' Standardized colors used for Lund Taxonomy subtypes and datasets.
#'
#' A list of 4 with color palettes frequently used in this package.
#'
#' \itemize{
#'  \item lund_colors. Color palette for the LundTax subtypes.
#'  \item lund_colors_transp. Same as `lund_colors` but transparent.
#'  \item stage_colors. Color palette for stages.
#'  \item dataset_colors. Color palette for the LundTax datasets.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name lund_colors
#' @usage data(lund_colors)
#' @format A list of 4.
NULL

#' Signatures
#'
#' Bundles a set of signatures in different molecular processes and the genes associated.
#' This dataset is needed for calculating the signature scores.
#'
#' A list of 6 with different signatures and the associated genes.
#'
#' \itemize{
#'  \item signatures_plot. Genes in hgnc and ensembl format associated with specific signatures.
#'  \item proliferation. Genes in hgnc and ensembl format associated with proliferation.
#'  \item progression. Genes in hgnc and ensembl format associated with progression.
#'  \item prostate. Genes in hgnc and ensembl format associated with prostate cancer.
#'  \item immune. Genes in hgnc and ensembl format associated with immuno process.
#'  \item stable_genes. Genes in hgnc and ensembl format associated with stable genes.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name signatures
#' @usage data(signatures)
#' @format A list of 6.
NULL

#' Sjodahl 2017.
#'
#' Gene expression data derived from the Sjödahl et. al. (2017) cohort. A matrix of RMA normalized
#' and ComBat adjusted gene expression values for 15697 genes (hgnc format) for 267 samples
#'
#' A data frame with gene information in different formats.
#'
#' @docType data
#' @keywords datasets
#' @name sjodahl_2017
#' @usage data(sjodahl_2017)
#' @format A data frame with 15697 rows (gene expressions) and 267 columns (samples).
NULL

#' Sjodahl 2017 Meta.
#'
#' metadata associated with the bundled expresison data.
#'
#' A data frame with meta data for the bundled samples.
#' 
#' \itemize{
#'  \item sample_id. Unique sample ID.
#'  \item cluster_order. Cluster order.
#'  \item age. Age of patient, years.
#'  \item gender. Patient gender.
#'  \item region_cx. Geographical region of CX.
#'  \item turb_grade. Pathological grading of TURB-specimen. 0 = WHO1999 Grade 2, 1 = WHO1999 Grade 3.
#'  \item turb_stage. Pathological staging of TURB-specimen. 0 = Non muscle-invasive, 1 = Muscle invasive.
#'  \item turb_skiv. Presence of keratinization or squamous differentiation in TURB-specimen. 0 = No presence, 1 = presence.
#'  \item turb_cis. Presence of cis in the TURB-specimen. 0 = No presence, 1 = presence.
#'  \item turb_lvi. Presence of lymphovascular invasion in the TURB-specimen. 0 = No presence, 1 = presence.
#'  \item turb_necros. Presence of necrosis in TURB-specimen. 0 = no, 1 = present, 2 = extensive.
#'  \item turb_apoptos. Presence of apoptosis in the TURB-specimen. 0 = No presence, 1 = presence.
#'  \item variant_histology. Variant histology present 1 = yes, 0 = no.
#'  \item cT. Clinical tumor state
#'  \item cN. Clinical node state
#'  \item hydroneph. Presence of hydronephrose. 0 = No presence, 1 = presence.
#'  \item pTCx. Pathological tumor state at cystectomy.
#'  \item pN. Patological node state at cystectomy.
#'  \item adj_chemo. Adjeuvant chemotherapy. 0 = No presence, 1 = presence.
#'  \item surv_os_event. Survival data, 1 = event, 0 = non event. 
#'  \item surv_os_time. Survival data, in months.
#'  \item surv_css_event. Cancer specific survival, 1 = event, 0 = non event. 
#'  \item surv_css_time. Cancer specific survival, in months.
#'  \item surv_pfs_event. Progression-free survival, 1 = event, 0 = non event. 
#'  \item surv_pfs_time. Progression-free survival, in months.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sjodahl_2017_meta
#' @usage data(sjodahl_2017_meta)
#' @format A data frame with 15697 rows (gene expressions) and 267 columns (samples).
NULL

#' GENCODE v44 Transcript-to-Gene Mapping Table
#'
#' A data frame mapping Ensembl transcript IDs and gene IDs to HGNC gene symbols,
#' based on GENCODE v44.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{transcript_id}{Ensembl transcript ID (e.g., "ENST00000426406")}
#'   \item{gene_id}{Ensembl gene ID (e.g., "ENSG00000284733")}
#'   \item{gene_name}{HGNC gene symbol (e.g., "OR4F29")}
#' }
#' @source \url{https://www.gencodegenes.org/human/release_44.html}
"tx2gene"
