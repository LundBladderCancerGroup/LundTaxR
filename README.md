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

LundTaxR implements a Random Forest-based classifier for molecular subtyping of bladder cancer samples using gene expression data. The package provides a two-stage classification system based on the established Lund Taxonomy. This classification can be applied to both non–muscle-invasive and muscle-invasive disease and provides consistent and biologically meaningful groups.

## Overviews

The Lund Taxonomy is a single-sample molecular classification system that divides bladder cancer into distinct subtypes with different clinical behaviors and treatment responses. LundTaxR automates this classification process through:

- **Two-stage classification**: First classifies samples into 5 main subtypes (Uro, GU, BaSq, Mes, ScNE), then sub-classifies Uro samples into UroA, UroB, or UroC
- **Comprehensive molecular profiling**: Calculates signature scores for immune infiltration, stromal content, cell proliferation, and progression risk
- **Advanced visualizations**: Generates publication-ready heatmaps and plots
- **Robust data handling**: Supports multiple gene ID formats and includes missing data imputation

## Key Features

### 🎯 **Accurate Classification**
- Random Forest-based predictors trained on extensive bladder cancer datasets
- Confidence scores and prediction reliability metrics
- Single sample classification

### 📊 **Molecular Signatures**
- Immune infiltration scores (CD8+ T cells, NK cells, macrophages)
- Stromal content assessment
- Cell proliferation and progression risk markers

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
Mattsson A, Aramendía E, Eriksson P, et al. (2024). 
LundTaxR: Molecular subtyping of bladder cancer using the Lund Taxonomy. 
R package version 2.0.0.
```

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
