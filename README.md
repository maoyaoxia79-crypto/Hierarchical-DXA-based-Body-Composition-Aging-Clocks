# Hierarchical DXA-based Body Composition Aging Clocks

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)

This repository provides the core R codebase for building and validating **Hierarchical DXA-based Body Composition Aging Clocks** and modeling biological age gaps.

---

## 🔬 Scientific Framework

Dual-energy X-ray Absorptiometry (DXA) measures regional fat, muscle, and bone mineral content. This project implements a hierarchical framework to analyze biological aging across multiple systems:

* **Level 0 (L0 - Biomarkers):** Regional DXA traits (axial/peripheral fat, muscle, and bone mineral content).
* **Level 1 (L1 - Domain Clocks):** Domain-specific biological age clocks constructed using Gompertz survival models and LASSO-Cox selection (Central Fat, Peripheral Fat, Total Fat, Axial Muscle, Peripheral Muscle, Total Muscle, Axial Skeletal, Peripheral Skeletal, Total Skeletal, and overall Body Composition).
* **Level 2 (L2 - Age Gaps):** Biological Age Gaps (BAG = Biological Age - Chronological Age) representing accelerated or decelerated aging in specific systems.
* **Level 3 (L3 - Causal Networks):** Causal network structure learning using PC/GES algorithms to map pathways between system BAGs, biomarkers, and clinical outcomes (mortality and multi-system diseases).

---

## 📂 Repository Structure

```
Hierarchical-DXA-based-Body-Composition-Aging-Clocks/
├── README.md                 # Project documentation
├── LICENSE                   # MIT License
└── code/
    ├── clock_construction_functions.R 
    │                         # Core R functions for clock calibration, selection, and metrics
    └── demo_simulation.R     # Self-contained script simulating cohort data and running the pipeline
```

* **`clock_construction_functions.R`**: Contains modular functions for collinearity filtering (`caret::findCorrelation`), LASSO variable selection (`glmnet::cv.glmnet`), Gompertz-Cox regression (`flexsurv::flexsurvreg`), biological age scoring, and performance evaluation.
* **`demo_simulation.R`**: Generates synthetic (mock) data representing standardized DXA traits, age, sex, and survival outcomes. It calls the core functions to build domain-specific clocks, outputs evaluation metrics, and saves performance summaries.

---

## 🛠️ Prerequisites & Installation

The code runs in **R** (version $\ge 4.0$). Install the required libraries from CRAN:

```R
# Install dependencies from CRAN
install.packages(c("survival", "flexsurv", "glmnet", "caret", "ggplot2", 
                   "Metrics", "ggExtra", "tidyr", "dplyr"))
```

---

## 🚀 Quick Start (Demo)

To test the clock construction pipeline on synthetic data, run the `demo_simulation.R` script:

1. Open R or RStudio.
2. Set the working directory to `code/`.
3. Source `demo_simulation.R`:

```R
# Set working directory to the code folder
setwd("path/to/repository/code")

# Run the simulation and biological clock building pipeline
source("demo_simulation.R")
```

### Outputs Generated:
* Model performance metrics (Correlation, MAE, R²) are summarized and saved to `output/All_Models_Performance_Summary.csv`.
* Visual scatter plots are saved as PDFs in the corresponding model output folders (e.g., `output/Central_Fat_Age_Results_staged_LASSO/`).

---

## 🔒 Data Availability

The clinical and DXA data used in this study were obtained from the **UK Biobank** under Application Number [Insert Application Number]. Due to licensing restrictions, raw individual-level data cannot be directly shared in this repository. External researchers can apply for access directly via the [UK Biobank Resource](https://www.ukbiobank.ac.uk/).

---

## 📝 License & Citation

This repository is licensed under the **MIT License**. If you use this code in your research, please cite our paper:

```
[Maoyao Xia, et al. 2026(IN SUBMISSION)]
```
