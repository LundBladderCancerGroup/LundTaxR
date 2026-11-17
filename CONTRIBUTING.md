# Git Workflow Guide for LundTaxR Contributors

This guide will help you contribute to the LundTaxR project, even if you have little or no experience with git or GitHub. Follow these steps to make changes and submit them for review.

## Review Objects
These are the following objects that needs review from you:

- [Advanced Classification Workflows](https://github.com/LundBladderCancerGroup/LundTaxR/blob/main/vignettes/advanced_classification_workflows.Rmd)

- [Data Retrieval and Utility Workflows](https://github.com/LundBladderCancerGroup/LundTaxR/blob/main/vignettes/data_retrieval_and_utility_workflows.Rmd)

- [Getting Started with LundTaxR](https://github.com/LundBladderCancerGroup/LundTaxR/blob/main/vignettes/getting_started_with_LundTaxR.Rmd)

- [Quality Control of Classification Results](https://github.com/LundBladderCancerGroup/LundTaxR/blob/main/vignettes/quality_control.Rmd)

- [Survival and Statistical Analysis Workflows](https://github.com/LundBladderCancerGroup/LundTaxR/blob/main/vignettes/survival_and_statistical_analysis_workflows.Rmd)

- [Visualization Workflows](https://github.com/LundBladderCancerGroup/LundTaxR/blob/main/vignettes/visualization_workflows.Rmd)

- [README](https://github.com/LundBladderCancerGroup/LundTaxR/blob/main/README.md)


## 1. Install Git
- Download and install git from https://git-scm.com/downloads
- On Windows, you may also want to install [GitHub Desktop](https://desktop.github.com/) for a graphical interface.

## 2. Clone the Repository
Open a terminal (or Git Bash on Windows) and run:

```
git clone https://github.com/LundBladderCancerGroup/LundTaxR.git
cd LundTaxR
```

## 3. Switch to Your Own Branch
In the same temrinal, run the following (replacing `name` with your name):

```
git checkout -b name-dev

#e.g
git checkout -b pontus-dev
git checkout -b carina-dev
git checkout -b gottfrid-dev

```

## 4. Make Your Changes
- Edit any code, vignettes, or the README as needed.
- Save your changes.

## 5. Stage and Commit Your Changes
```
git add .
git commit -m "Describe your changes here"
```

## 6. Push Your Branch to GitHub
```
git push origin your-branch-name
```

## 7. Create a Pull Request (PR)
- Go to https://github.com/LundBladderCancerGroup/LundTaxR in your web browser.
- You should see a prompt to "Compare & pull request" for your branch. Click it.
- Add a description of your changes and submit the pull request.

## 8. Wait for Review
- The project maintainer will review your PR and may request changes or approve it.

---

### Tips
- If you get stuck, ask for help in the project issues or by email.
- Always work on your own branch, not `main`.
- Pull the latest changes from `main` before starting new work:
  ```
  git checkout main
  git pull origin main
  git checkout your-branch-name
  git merge main
  ```

Thank you for contributing to LundTaxR!
