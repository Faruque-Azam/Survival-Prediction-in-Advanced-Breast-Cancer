# R Package Requirements for Survival Prediction in Advanced Breast Cancer
# Install using: source("install_packages.R")

packages <- c(
  # Data manipulation
  "tidyverse",
  "readxl",
  
  # Missing data handling
  "mice",
  "VIM",
  
 # Machine learning (tidymodels ecosystem)
  "tidymodels",
  "ranger",
  "rpart",
  "rpart.plot",
  "vip",
  "finetune",
  "xgboost",
  
  # Visualization
  "ggplot2",
  "ggstatsplot",
  "viridis",
  "cowplot",
  "patchwork",
  "ggtext",
  
  # Parallel processing
  "doParallel"
)

# Install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(sapply(packages, install_if_missing))

cat("All required packages installed successfully.\n")
