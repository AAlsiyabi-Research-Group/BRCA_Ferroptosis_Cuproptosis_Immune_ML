# BRCA_Ferroptosis_Cuproptosis_Immune_ML
Title: An Integrative machine learning approach identifies the centrality of ferroptosis, cuproptosis, and immune pathway crosstalk for breast cancer stratification and therapy guidance.

**Overview**
This repository contains the full reproducibility package for the manuscript:

"An integrative machine learning approach identifies the centrality of ferroptosis, cuproptosis, and immune pathway crosstalk for breast cancer stratification and therapy guidance".

We present a computational framework integrating:
  Differential expression analysis (DESeq2)
  Pathway-focused Wilcoxon testing (ferroptosis, cuproptosis, immune genes)
  Machine learning classifiers (RF, SVM, XGB, GBM, AdaBoost, KNN, ANN, LASSO, etc.)
  SHAP interpretability and feature importance
  PCA and violin visualization of a robust 4-gene panel (FOXO4, EGFR, FGF2, CDKN2A)
  Immune infiltration analysis using xCell
All analyses are performed on TCGA-BRCA RNA-seq data (public, 1,111 tumors + 113 normals).
No wet-lab experiments were performed.

**Repository Structure**
├── data/                  # TCGA-BRCA raw counts, metadata, VST-normalized matrix
├── scripts/               # R and Python scripts, organized by step
│   ├── R/                 # DESeq2, Wilcoxon, immune infiltration
│   └── Python/            # ML notebooks, SHAP, PCA, plots
├── results/               # Analysis outputs, metrics, tables
├── figures/               # Uncropped original plots (PNG/TIFF)
└── README.md              # This file

**Installation & Requirements**

R (≥ 4.3.0)
install.packages(c("BiocManager","data.table","ggplot2","reshape2"))
BiocManager::install(c("TCGAbiolinks","SummarizedExperiment","DESeq2","apeglm","biomaRt","pheatmap","xCell"))

Python (≥ 3.10)
pip install numpy pandas scikit-learn xgboost imbalanced-learn shap matplotlib seaborn


  **How to Reproduce**
  Download TCGA-BRCA data
  Use 1.BRCA-data_retrieval.R (or provided CSVs in data/).
  
  **Run differential expression (DESeq2)**
  2A_BRCA_deseq2_DEG_analysis.R
  Outputs DEG tables, volcano data, heatmaps.

  **Run Wilcoxon analyses**
  Ferroptosis: 5A_BRCA_Wilcoxon_Ferroptosis_Analysis.R
  Cuproptosis: 5B_BRCA_Wilcoxon_Cuproptosis_Analysis.R
  Immune: 5C_BRCA_wilcoxon_immune_analysis.R

 **Run machine learning panel discovery**
  Open 6.BRCA_4_genes_ML.ipynb
  Check results/ for model metrics and performance plots.
  
 **Run SHAP interpretability & PCA/violin**
  8.SHAP/, 9.Plots.ipynb.

**Run immune infiltration analysis (xCell)**
  10.Immune_Infiltration_Analysis.R


**Key Results**
  A compact 4-gene panel (FOXO4, EGFR, FGF2, CDKN2A) robustly stratifies tumor vs normal.
  The panel integrates ferroptosis, cuproptosis, and immune pathway crosstalk.
  Strong correlations observed with immune cell infiltration (xCell).
  Figures provided in /figures/ correspond to manuscript Figures 2–8.
