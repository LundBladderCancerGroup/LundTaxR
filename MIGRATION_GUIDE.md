# LundTaxR Repository Migration Guide

> **Migration from LundTax2023Classifier to LundTaxR and transfer to research group GitHub account**

## Overview

This guide outlines the complete process for migrating from the old `LundTax2023Classifier` package to the new `LundTaxR` package and transferring it to the research group's GitHub account.

## 📋 Pre-Migration Checklist

- [ ] Ensure all functionality from old package is implemented in LundTaxR
- [ ] All tests are passing
- [ ] Documentation is complete and accurate
- [ ] GitHub Actions workflow is working
- [ ] Website builds successfully
- [ ] Version number reflects major change (2.0.0+)

---

## Phase 1: Prepare LundTaxR Package

### 1.1 Package Metadata Updates

**File: `DESCRIPTION`**
```r
Package: LundTaxR
Title: Molecular Subtype Classification for Bladder Cancer
Version: 2.0.0
Description: Molecular subtype classification for bladder cancer using the Lund 
    Taxonomy system. This package replaces LundTax2023Classifier with improved 
    functionality, cleaner design, and modern R package standards.
URL: https://lundbladdercancergroup.github.io/LundTaxR, https://github.com/LundBladderCancerGroup/LundTaxR
BugReports: https://github.com/LundBladderCancerGroup/LundTaxR/issues
```

### 1.2 Package-Level Documentation

**File: `R/LundTaxR-package.R`** (create if doesn't exist)
```r
#' @details 
#' LundTaxR replaces the previous LundTax2023Classifier package with improved 
#' functionality, cleaner design, and enhanced documentation. Users migrating 
#' from the previous package will find similar functionality with updated 
#' function names and improved performance.
#' 
#' For migration guidance, see the "Getting Started" vignette.
#' 
#' @keywords internal
"_PACKAGE"
```

### 1.3 Main Function Documentation

**File: `R/classify_samples.R`** (add to existing documentation)
```r
#' @details 
#' This is the main classification function in LundTaxR, replacing the 
#' functionality previously available in LundTax2023Classifier with 
#' enhanced features and improved usability.
```

### 1.4 Homepage Updates

**File: `pkgdown/index.md`** (update existing file)
```markdown
> **Robust molecular subtyping of bladder cancer using the Lund Taxonomy classification system**

*LundTaxR is the successor to LundTax2023Classifier with enhanced functionality and improved design.*

## Migration from LundTax2023Classifier

If you're upgrading from the previous package:

### Installation
```r
# Remove old package (optional)
remove.packages("LundTax2023Classifier")

# Install new package
devtools::install_github("LundBladderCancerGroup/LundTaxR")
```

### Key Changes
- **Package name**: `LundTax2023Classifier` → `LundTaxR`
- **Main function**: Enhanced `classify_samples()` with improved parameters
- **Visualization**: New plotting functions with better customization
- **Documentation**: Comprehensive vignettes and examples
- **Performance**: Optimized algorithms and reduced dependencies
```

---

## Phase 2: Deprecate Old Package

### 2.1 Add Deprecation Warning

**File: `R/zzz.R` in LundTax2023Classifier repo**
```r
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    crayon::red("WARNING: This package has been replaced by 'LundTaxR'."), "\n",
    "Please install the new version: ", crayon::green("devtools::install_github('LundBladderCancerGroup/LundTaxR')"), "\n",
    "Documentation: ", crayon::blue("https://lundbladdercancergroup.github.io/LundTaxR"), "\n",
    "GitHub: ", crayon::blue("https://github.com/LundBladderCancerGroup/LundTaxR"), "\n",
    "\n",
    "LundTaxR offers improved functionality, better documentation, and ongoing support.\n"
  )
}
```

### 2.2 Update Old Package DESCRIPTION

```r
Package: LundTax2023Classifier
Title: [DEPRECATED] Use LundTaxR instead
Version: 1.9.9  # Increment to show it's the final version
Description: This package has been replaced by LundTaxR. Please install 
    the new package: devtools::install_github('LundBladderCancerGroup/LundTaxR').
    LundTaxR provides the same functionality with improvements and ongoing support.
URL: https://github.com/LundBladderCancerGroup/LundTaxR
BugReports: https://github.com/LundBladderCancerGroup/LundTaxR/issues
Imports: crayon
```

### 2.3 Update Old Package README

```markdown
# ⚠️ DEPRECATED PACKAGE ⚠️

**This package has been replaced by [LundTaxR](https://github.com/LundBladderCancerGroup/LundTaxR).**

## 🚀 Install the new package:
```r
devtools::install_github('LundBladderCancerGroup/LundTaxR')
```

## ✨ What's new in LundTaxR:
- **Improved performance** and cleaner code architecture
- **Enhanced visualizations** with better customization options
- **Comprehensive documentation** with detailed vignettes
- **Active maintenance** and ongoing support
- **Modern R package standards** and best practices

## 📚 Resources:
- **Documentation**: https://lundbladdercancergroup.github.io/LundTaxR
- **GitHub Repository**: https://github.com/LundBladderCancerGroup/LundTaxR
- **Getting Started Guide**: https://lundbladdercancergroup.github.io/LundTaxR/articles/getting_started_with_LundTaxR.html

---

*This repository is maintained for archival purposes only. All new development happens in LundTaxR.*
```

---

## Phase 3: Repository Transfer

### 3.1 Pre-Transfer Preparation

1. **Verify GitHub Actions are working**
   ```bash
   git push origin main
   # Check Actions tab for successful build
   ```

2. **Document current URLs** (for later updates)
   - Repository: `https://github.com/mattssca/LundTaxR`
   - Website: `https://mattssca.github.io/LundTaxR`
   - Clone URL: `https://github.com/mattssca/LundTaxR.git`

### 3.2 GitHub Repository Transfer

1. **Go to repository settings**: `https://github.com/mattssca/LundTaxR/settings`
2. **Scroll to "Danger Zone"**
3. **Click "Transfer ownership"**
4. **Enter target organization name**: `[RESEARCH-GROUP-NAME]`
5. **Type repository name to confirm**: `LundTaxR`
6. **Click "I understand, transfer this repository"**

### 3.3 Post-Transfer Updates

**Immediately after transfer:**

1. **Update _pkgdown.yml**
   ```yaml
   url: https://lundbladdercancergroup.github.io/LundTaxR
   
   navbar:
     components:
       github:
         icon: fab-github-alt
         href: https://github.com/LundBladderCancerGroup/LundTaxR
   ```

2. **Update DESCRIPTION**
   ```r
   URL: https://lundbladdercancergroup.github.io/LundTaxR, https://github.com/LundBladderCancerGroup/LundTaxR
   BugReports: https://github.com/LundBladderCancerGroup/LundTaxR/issues
   ```

3. **Update README.md installation instructions**
   ```r
   devtools::install_github("LundBladderCancerGroup/LundTaxR")
   ```

4. **Update pkgdown/index.md installation instructions**
   ```r
   devtools::install_github("LundBladderCancerGroup/LundTaxR")
   ```

### 3.4 Verify GitHub Pages

1. **Check Pages settings**: `Settings → Pages`
2. **Ensure Source is "GitHub Actions"**
3. **Verify site builds at new URL**: `https://lundbladdercancergroup.github.io/LundTaxR`

### 3.5 Update Local Development Environment

```bash
# Update remote URL for local repo
git remote set-url origin https://github.com/LundBladderCancerGroup/LundTaxR.git

# Verify the change
git remote -v

# Pull any changes
git pull origin main
```

---

## Phase 4: Communication & Rollout

### 4.1 Update Documentation References

**Files to update with new URLs:**
- [ ] All vignettes (`vignettes/`)
- [ ] Function documentation referencing GitHub
- [ ] Any citation files
- [ ] Conference presentations/papers

### 4.2 Communication Timeline

**Week 1-2: Soft Launch**
- [ ] Update old package with deprecation warnings
- [ ] Notify close collaborators
- [ ] Test new installation process

**Week 3-4: Public Announcement**
- [ ] Social media announcements
- [ ] Email to user mailing lists
- [ ] Update any published papers/preprints

**Month 2-3: Full Transition**
- [ ] Archive old repository (add archive banner)
- [ ] Remove old package from active development

### 4.3 User Migration Support

**Create migration checklist for users:**

```r
# 1. Remove old package (optional)
remove.packages("LundTax2023Classifier")

# 2. Install new package  
devtools::install_github("LundBladderCancerGroup/LundTaxR")

# 3. Update library calls
library(LundTaxR)  # instead of library(LundTax2023Classifier)

# 4. Update function calls (check documentation for any changes)
results <- classify_samples(data)  # Main function likely unchanged

# 5. Update any saved scripts with new package name
```

---

## Phase 5: Cleanup & Maintenance

### 5.1 Archive Old Repository

**Add to old repository README:**
```markdown
# 📦 ARCHIVED REPOSITORY

This repository has been archived and is no longer maintained.

**➡️ Active development has moved to: https://github.com/LundBladderCancerGroup/LundTaxR**
```

**GitHub settings:**
1. Go to old repo settings
2. Scroll to "Danger Zone"  
3. Click "Archive this repository"
4. Confirm archival

### 5.2 Monitor Transition

**Track for first 3 months:**
- [ ] GitHub stars/forks on new repo
- [ ] Issues opened (migration problems?)
- [ ] Download statistics
- [ ] User feedback

### 5.3 Update Citations

**Template for papers/presentations:**
```
We used the LundTaxR R package (version 2.0.0, https://github.com/LundBladderCancerGroup/LundTaxR) 
for molecular subtype classification of bladder cancer samples.
```

---

## 🚨 Troubleshooting

### Common Issues After Transfer

**GitHub Actions failing:**
- Check that secrets/tokens transferred correctly
- Verify organization permissions for Actions
- Update any hardcoded repository references in workflows

**GitHub Pages not working:**
- Confirm Pages is enabled for the organization
- Check that the `gh-pages` branch exists
- Verify GitHub Actions has Pages deployment permissions

**Installation issues:**
- Ensure repository is public (or users have access)
- Check that release tags transferred correctly
- Verify DESCRIPTION file URLs are updated

### Rollback Plan

If major issues arise:
1. Keep old repository accessible (don't archive immediately)
2. Document issues and timeline for fixes
3. Communicate clearly with users about temporary solutions
4. Consider gradual migration if needed

---

## 📊 Success Metrics

**Migration considered successful when:**
- [ ] New package installs cleanly from new location
- [ ] Website builds and displays correctly at new URL
- [ ] All GitHub Actions workflows pass
- [ ] No broken links in documentation
- [ ] Users can follow migration guide successfully
- [ ] Old package shows appropriate deprecation warnings

---

## 📝 Notes

**Important decisions to make:**
- [ ] **Research group GitHub organization name**
- [ ] **Exact timeline for deprecation phases**
- [ ] **Communication channels for announcements**
- [ ] **Whether to maintain old package with bug fixes during transition**

**Post-migration TODO:**
- [ ] Add team members to new repository with appropriate permissions
- [ ] Set up branch protection rules if needed
- [ ] Configure organization-level settings (Actions permissions, etc.)
- [ ] Update any CI/CD integrations or external services

---

*Last updated: October 16, 2025*
