# --------------------------------------------
# BRCA_Cuproptosis_Wilcoxon_Analysis.R
# Wilcoxon Test for Cuproptosis Genes + Heatmap + Normalized Matrix
# --------------------------------------------

# Load libraries
library(readr)
library(dplyr)
library(pheatmap)

# --------------------------------------------
# Step 1: Load Normalized VST Data
# --------------------------------------------
vst_mat <- read.csv("/Users/apple/Desktop/BRCA/BRCA_VST_Normalized_Matrix.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("/Users/apple/Desktop/BRCA/BRCA_Metadata_Final.csv", row.names = 1)

# Ensure sample order is correct
vst_mat <- vst_mat[, rownames(meta)]
condition <- factor(ifelse(meta$shortLetterCode == "TP", "Tumor", "Normal"))

# --------------------------------------------
# Step 2: Load Cuproptosis Genes
# --------------------------------------------
cupro_genes <- read_delim("Cuproptosis_GeneList_133.txt", delim = "\t", col_names = FALSE)
gene_list <- cupro_genes$X1

# Match cuproptosis genes with normalized data
matched_genes <- intersect(gene_list, rownames(vst_mat))
vst_cupro <- vst_mat[matched_genes, , drop = FALSE]

# --------------------------------------------
# Step 3: Wilcoxon Test
# --------------------------------------------
wilcox_results <- lapply(matched_genes, function(gene) {
  exp_vals <- as.numeric(vst_cupro[gene, ])
  group1 <- exp_vals[condition == "Tumor"]
  group2 <- exp_vals[condition == "Normal"]
  test <- wilcox.test(group1, group2)
  logfc <- median(group1) - median(group2)
  data.frame(
    Gene = gene,
    Median_Tumor = median(group1),
    Median_Normal = median(group2),
    LogFC = logfc,
    p_value = test$p.value
  )
})

wilcox_df <- bind_rows(wilcox_results)
wilcox_df$adj_p <- p.adjust(wilcox_df$p_value, method = "fdr")
wilcox_df <- wilcox_df[order(wilcox_df$adj_p), ]

# Save full results
write.csv(wilcox_df, "BRCA_Wilcoxon_Cuproptosis.csv", row.names = FALSE)

# --------------------------------------------
# Step 4: Save Significant Genes (Filtered)
# --------------------------------------------
logfc_cutoff <- 1
fdr_cutoff <- 0.05
sig_genes <- wilcox_df %>%
  filter(abs(LogFC) > logfc_cutoff & adj_p < fdr_cutoff)
write.csv(sig_genes, "BRCA_Significant_Cuproptosis_Genes.csv", row.names = FALSE)

# --------------------------------------------
# Step 5: Heatmap of Top 30 Significant Genes
# --------------------------------------------
top_genes <- head(sig_genes$Gene, 30)
write.csv(data.frame(Gene = top_genes), "BRCA_Top30_Cuproptosis_Genes.csv", row.names = FALSE)

annotation_df <- data.frame(condition = condition)
rownames(annotation_df) <- colnames(vst_mat)

jpeg("BRCA_Cuproptosis_Heatmap.jpg", width = 1000, height = 800, quality = 600)
pheatmap(
  vst_cupro[top_genes, , drop = FALSE],
  annotation_col = annotation_df,
  show_colnames = FALSE,
  scale = "row",
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()

# --------------------------------------------
# Final Output Summary
# --------------------------------------------
cat("\n✅ Wilcoxon analysis complete.\n")
cat("• Wilcoxon results: 'BRCA_Wilcoxon_Cuproptosis.csv'\n")
cat("• Significant genes: 'BRCA_Significant_Cuproptosis_Genes.csv'\n")
cat("• Top 30 genes: 'BRCA_Top30_Cuproptosis_Genes.csv'\n")
cat("• Heatmap saved as 'BRCA_Cuproptosis_Heatmap.jpg'\n")