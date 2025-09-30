# ===============================
# 📦 Step 0: Install & Load Packages
# ===============================
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("xCell", quietly = TRUE)) remotes::install_github("dviraran/xCell")

library(xCell)
library(readr)
library(dplyr)
library(pheatmap)

# ===============================
# 📁 Step 1: Load VST Expression Matrix
# ===============================
expr_file <- "/Users/apple/Desktop/BRCA/BRCA_VST_Normalized_Matrix.csv"
expr_mat <- read.csv(expr_file, row.names = 1, check.names = FALSE)

# Ensure it's genes x samples
if (ncol(expr_mat) > nrow(expr_mat)) {
  expr_mat <- t(expr_mat)
}

cat("✅ Expression matrix loaded\n")
cat("Num. of genes:", nrow(expr_mat), "\n")

# ===============================
# 🧬 Step 2: Run xCell Analysis
# ===============================
xcell_scores <- xCellAnalysis(as.matrix(expr_mat))  # cell types x samples

# ===============================
# 🔎 Step 3: Extract 4 Minimal Panel Genes
# ===============================
final_genes <- c("FOXO4", "EGFR", "FGF2", "CDKN2A")
gene_expr <- expr_mat[rownames(expr_mat) %in% final_genes, , drop = FALSE]

# ===============================
# 📊 Step 4: Correlation (Spearman)
# ===============================
common_samples <- intersect(colnames(gene_expr), colnames(xcell_scores))
gene_expr <- gene_expr[, common_samples]
xcell_scores <- xcell_scores[, common_samples]

# Transpose xCell scores so: samples x cell types
xcell_scores_t <- t(xcell_scores)

# Correlation: gene expression vs immune cell infiltration
cor_mat <- cor(t(gene_expr), xcell_scores_t, method = "spearman")

# Save correlation matrix
cor_file <- "/Users/apple/Desktop/BRCA/BRCA_4Gene_xCell_Correlation_Matrix.csv"
write.csv(cor_mat, file = cor_file)
cat("✅ Correlation matrix saved to:", cor_file, "\n")

# ===============================
# 🎨 Step 5: Heatmap (600 DPI)
# ===============================
heatmap_file <- "/Users/apple/Desktop/BRCA/BRCA_4Gene_ImmuneCorrelation_Heatmap.png"

png(heatmap_file, width = 12, height = 8, units = "in", res = 600)
pheatmap(cor_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 10,
         main = "Spearman Correlation: 4-Gene Panel vs Immune Cell Infiltration (xCell)",
         border_color = NA)
dev.off()

cat("✅ Heatmap saved to:", heatmap_file, "\n")