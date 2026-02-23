# Survival Prediction in Advanced Breast Cancer Trials

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the data and code for the study entitled **"Interpretable Machine Learning to Predict Progression-free and Overall Survival in Advanced Breast Cancer Clinical Trials"**

We analyzed 2,622 efficacy endpoints from 892 phase II trials in advanced breast cancer to:
1. Examine correlations among surrogate endpoints (ORR, mPFS, mOS)
2. Develop interpretable machine learning models to predict median progression-free survival (mPFS) and median overall survival (mOS)

## Key Findings

- **Sequential endpoint pattern**: ORR correlates strongly with mPFS (ρ = 0.72), and mPFS with mOS (ρ = 0.68), while direct ORR-mOS correlation is weaker (ρ = 0.45)
- **Predictive performance**: Random forest achieved test MAE of 1.58 months for mPFS and 3.96 months for mOS
- **Interpretable hierarchy**: Decision trees identified ORR as the dominant predictor for mPFS, while mPFS — not ORR — was the principal driver of mOS

## Repository Structure

```
├── data/
│   ├── mBC_Dataset.xlsx       # Main dataset (1,088 trial arms, 25 variables)
│   └── mBC_Combo.xlsx         # Treatment combination flags (22 drug class indicators)
├── code/
│   └── Survival_Modeling_mBC.Rmd   # Complete analysis pipeline
├── figures/                   # Generated figures (after running code)
├── README.md
├── LICENSE
└── .gitignore
```

## Data Description

### mBC_Dataset.xlsx

Primary dataset containing trial-level efficacy endpoints and characteristics:

| Variable | Description |
|----------|-------------|
| PMID | PubMed identifier |
| Year | Publication year |
| Subtype | Tumor subtype (HR+, HER2+, Mixed, TNBC) |
| Pretreated | Line of therapy (First-line, Pretreated) |
| Age | Median patient age |
| N | Sample size per trial/arm|
| wECOG | Weighted ECOG performance status |
| ORR | Objective response rate (%) |
| mPFS | Median progression-free survival (months) |
| mOS | Median overall survival (months) |
| Therapy | Treatment class (Chemotherapy, Targeted, Chemo-Targeted) |
| Size | Number of agents (1-Agent, 2-Agent, 3-Agent) |
| T1–T4 | Individual treatment names |

### mBC_Combo.xlsx

Binary indicators for 22 drug classes:

- **Endocrine**: Antiestrogen, Aromatase inhibitor, GnRH agonist
- **Targeted**: CDK4/6 inhibitor, mTOR inhibitor, PARP inhibitor, HER2 antibody, HER2 kinase inhibitor, EGFR kinase inhibitor, VEGF inhibitor, Checkpoint inhibitor
- **Chemotherapy**: Taxane, Platinum, Anthracycline, Antimetabolite, Alkylating, Vinca alkaloid, Topoisomerase inhibitor

## Requirements

### R Version
- R ≥ 4.0

### Required Packages

```r
# Data manipulation
install.packages(c("tidyverse", "readxl"))

# Missing data handling
install.packages(c("mice", "VIM"))

# Machine learning
install.packages(c("tidymodels", "ranger", "rpart", "rpart.plot", "vip", "finetune"))

# Visualization
install.packages(c("ggplot2", "ggstatsplot", "viridis", "cowplot", "patchwork"))

# Parallel processing
install.packages("doParallel")
```

## Usage

1. Clone the repository:
```bash
git clone https://github.com/[username]/survival-prediction-breast-cancer.git
cd survival-prediction-breast-cancer
```

2. Open the R Markdown file in RStudio:
```r
file.edit("code/Survival_Modeling_mBC.Rmd")
```

3. Set working directory and run the analysis:
```r
setwd("path/to/repository")
rmarkdown::render("code/Survival_Modeling_mBC.Rmd")
```

## Analysis Pipeline

The R Markdown file contains the following sections:

1. **Data Import & Preprocessing**: Load datasets, format variables
2. **Missing Data Imputation**: MICE with CART algorithm (m=20, maxit=30)
3. **Feature Engineering**: Create derived variables, combine datasets
4. **Imputation Quality Assessment**: Density plots comparing distributions
5. **Correlation Analysis**: Spearman correlations among ORR, mPFS, mOS
6. **mOS Modeling**:
   - mOS random forest model
   - mOS decision tree model 
7. **mPFS Modeling**:
   - mPFS random forest model
   - mPFS decision tree model

## Reproducibility

- Random seed: `set.seed(786)` for train/test splits; `set.seed(329)` for imputation
- Train/test split: 80/20 with stratification
- Bootstrap resamples: 30 iterations for hyperparameter tuning

## Citation

If you use this code or data, please cite:

```
[Author names]. Interpretable Machine Learning to Predict Progression-free and Overall Survival in Advanced Breast Cancer Clinical Trials. [Journal]. [Year]. DOI: [pending]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaboration inquiries, please contact:
Faruque Azam at walid2372@gmail.com

## Acknowledgments

- Data extracted from published phase II clinical trials (1992–2022)
- Analysis conducted using R and the tidymodels framework
