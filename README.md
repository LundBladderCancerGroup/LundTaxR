# LundTaxR

<img src="man/figures/logo.png" align="right" height="280" alt="LundTaxR logo" />

<!-- GitHub Actions Badges -->

[![R-CMD-check](https://github.com/mattssca/LundTaxR/actions/workflows/r-cmd-check.yml/badge.svg)](https://github.com/mattssca/LundTaxR/actions/workflows/r-cmd-check.yml)
[![testthat](https://github.com/mattssca/LundTaxR/actions/workflows/testthat.yml/badge.svg)](https://github.com/mattssca/LundTaxR/actions/workflows/testthat.yml)
[![pkgdown](https://github.com/mattssca/LundTaxR/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/mattssca/LundTaxR/actions/workflows/pkgdown.yaml)

<!-- R Package Badges -->

![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)
![License](https://img.shields.io/badge/license-GPL%20(%E2%89%A5%202)-blue)
![Version](https://img.shields.io/badge/version-2.0.0-blue)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen)

<!-- Domain Badges -->

![Bioinformatics](https://img.shields.io/badge/domain-bioinformatics-purple)
![Cancer Research](https://img.shields.io/badge/application-cancer%20research-red)
![Classification](https://img.shields.io/badge/method-machine%20learning-orange)

> **Robust molecular subtyping of bladder cancer using the Lund Taxonomy classification system**

LundTaxR implements a Random Forest-based classifier for molecular subtyping of bladder cancer samples using gene expression data. The package provides a two-stage classification system based on the established Lund Taxonomy, enabling precise molecular characterization.

## LundTax Subtypes

The LundTax subtypes are defined at the cancer phenotype level and does not include features of the microenvironment. Thus, cancer cell phenotypes are treated independent from features from the microenvironment.

## Overviews

The Lund Taxonomy is a single-sample molecular classification system that divides bladder cancer into distinct subtypes with different clinical behaviors and treatment responses. LundTaxR automates this classification process through:

- **Two-stage classification**: First classifies samples into 5 main subtypes (Uro, GU, BaSq, Mes, ScNE), then sub-classifies Uro samples into UroA, UroB, or UroC
- **Comprehensive molecular profiling**: Calculates signature scores for molecular grade, cell proliferation, progression risk, as well as signature scores for various immune and stromal cell types
- **Advanced visualizations**: Generates publication-ready heatmaps and plots
- **Robust data handling**: Independent of preprocessing, supports multiple gene ID formats, and includes missing data imputation

## Key Features

### 🎯 **Accurate Classification**

- Single sample classification
- Random Forest-based predictors trained on extensive bladder cancer datasets
- Confidence scores and prediction reliability metrics

### 📊 **Molecular Signatures**

- Molecular grade, proliferation and progression scores
- Scores for individual immune cell types (e.g. CD8+ T cells, NK cells, macrophages)
- Scores for individual stromal cell types

### 📈 **Visualization Tools**

- Classification result heatmaps
- Signature score distributions
- Subtype comparison plots
- Survival analysis visualizations

### 🔧 **User-Friendly**

- Simple single-function classification
- Comprehensive documentation and tutorials
- Flexible input formats
- Integration with standard R workflows

## Installation

### From GitHub (Recommended)

```r
# Install devtools if you haven't already
if (!require(devtools)) install.packages("devtools")

# Install LundTaxR
devtools::install_github("mattssca/LundTaxR")
```

### System Requirements

- R ≥ 4.0.0
- Required packages will be installed automatically

## Quick Start

```r
library(LundTaxR)

# 1. Load your gene expression data (genes as rows, samples as columns)
data("sjodahl_2017") # Bundled dataset

# 2. Run Classifier
sjodahl_classes = classify_samples(this_data = sjodahl_2017, 
                                   log_transform = FALSE, 
                                   adjust = TRUE, 
                                   impute = TRUE, 
                                   include_data = TRUE)

```

## Advanced Usage

For a more comprehensive tutorial and usage examples, please refer to the vignettes.

## Classification System

### Main Subtypes (5-class)

- **Uro** (Urothelial-like): Luminal, resembles urothelium
- **GU** (Genomically Unstable): Luminal, complex genome, lacking basal cells
- **BaSq** (Basal/Squamous): Non-luminal, basal- or squamous-like characteristics
- **Mes** (Mesenchymal): Non-luminal, mesenchymal/sarcomatoid cancer cells or extremely stroma-rich biopsy
- **ScNE** (Small cell/Neuroendocrine): Non-luminal, very aggressive subtype often with neuroendocrine histology

### Uro Subclasses (7-class)

- **UroA**: Intact basal cell stratification
- **UroB**: Extended basal cell stratification
- **UroC**: Absent or discontinuous stratification

## Documentation

- **Package Website**: [https://mattssca.github.io/LundTaxR/](https://mattssca.github.io/LundTaxR/) - Complete documentation, tutorials, and function reference
- **Package Documentation**: Access help files with `?function_name`
- **Vignettes**: Comprehensive tutorials and examples

```r
# View main tutorial
vignette("getting_started_with_LundTaxR_fig_1", package = "LundTaxR")

# Browse all documentation
help(package = "LundTaxR")
```

## Example Workflows

### Basic Classification

```r
# Classify samples
my_predicted = classify_samples(this_data = expression_data)
```

## Citation

If you use LundTaxR in your research, please cite:

```
Mattsson A, Aramendía E, Eriksson P, Sjödahl G, Bernardo C, Höglund M (2024). 
LundTaxR: Molecular subtyping of bladder cancer using the Lund Taxonomy. 
R package version 2.0.0.
```

## References

Sjödahl, G., Lauss, M., Lövgren, K., Chebil, G., Gudjonsson, S., Veerla, S., Patschan, O., Aine, M., Fernö, M., Ringnér, M., Månsson, W., Liedberg, F., Lindgren, D., & Höglund, M. (2012). A molecular taxonomy for urothelial carcinoma. *Clinical Cancer Research*, 18(12), 3377–3386. https://doi.org/10.1158/1078-0432.CCR-12-0077-T

Sjödahl, G., Eriksson, P., Liedberg, F., & Höglund, M. (2017). Molecular classification of urothelial carcinoma: global mRNA classification versus tumour-cell phenotype classification. *The Journal of Pathology*, 242(1), 113–125. https://doi.org/10.1002/path.4886

Marzouka, N. A., Eriksson, P., Rovira, C., Liedberg, F., Sjödahl, G., & Höglund, M. (2018). A validation and extended description of the Lund taxonomy for urothelial carcinoma using the TCGA cohort. *Scientific Reports*, 8, 3737. https://doi.org/10.1038/s41598-018-22126-x

Bernardo, C., Eriksson, P., Marzouka, N. A., Liedberg, F., Sjödahl, G., & Höglund, M. (2019). Molecular pathology of the luminal class of urothelial tumors. *The Journal of Pathology*, 249(3), 308–318. https://doi.org/10.1002/path.5318

Marzouka, N. A., & Eriksson, P. (2021). multiclassPairs: an R package to train multiclass pair-based classifier. *Bioinformatics*, 37(18), 3043–3044. https://doi.org/10.1093/bioinformatics/btab088

Bernardo, C., Eriksson, P., Marzouka, N. A., Liedberg, F., Sjödahl, G., & Höglund, M. (2022). Molecular pathology of the non-luminal Ba/Sq-like and Sc/NE-like classes of urothelial tumours: An integrated immunohistochemical analysis. *Human Pathology*, 122, 11–24. https://doi.org/10.1016/j.humpath.2022.01.006

Eriksson, P., Marzouka, N. A., Sjödahl, G., Bernardo, C., Liedberg, F., & Höglund, M. (2022). A comparison of rule-based and centroid single-sample multiclass predictors for transcriptomic classification. *Bioinformatics*, 38(4), 1022–1029. https://doi.org/10.1093/bioinformatics/btab725

Marzouka, N. A., Eriksson, P., Bernardo, C., Hurst, C. D., Knowles, M. A., Sjödahl, G., Liedberg, F., & Höglund, M. (2022). The Lund Molecular Taxonomy Applied to Non-Muscle-Invasive Urothelial Carcinoma. *The Journal of Molecular Diagnostics*, 24(9), 992–1008. https://doi.org/10.1016/j.jmoldx.2022.05.006

Höglund, M., Bernardo, C., Sjödahl, G., Eriksson, P., Axelson, H., & Liedberg, F. (2023). The Lund taxonomy for bladder cancer classification - from gene expression clustering to cancer cell molecular phenotypes, and back again. *The Journal of Pathology*, 259(4), 369–375. https://doi.org/10.1002/path.6062

Cotillas, E. A., Bernardo, C., Veerla, S., Liedberg, F., Sjödahl, G., & Eriksson, P. (2024). A Versatile and Upgraded Version of the LundTax Classification Algorithm Applied to Independent Cohorts. *The Journal of Molecular Diagnostics*, 26(12), 1081–1101. https://doi.org/10.1016/j.jmoldx.2024.08.005

## Contributing

We welcome contributions!

- Reporting bugs
- Suggesting enhancements
- Submitting pull requests

## Support

- **Issues**: [GitHub Issues](https://github.com/mattssca/LundTaxR/issues)
- **Email**: adam.mattsson@med.lu.se
- **Documentation**: Package help files and vignettes

## License

This project is licensed under the GPL (≥ 2) License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Lund Bladder Cancer Group
- Contributors to the original Lund Taxonomy classification system
- The R community for excellent bioinformatics tools

---

**Developed by the Lund Bladder Cancer Group** 🇸🇪
