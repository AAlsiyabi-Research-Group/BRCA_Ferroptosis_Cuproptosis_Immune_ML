# --------------------------------------------
# BRCA_Ferroptosis_Wilcoxon_Analysis.R (Updated with output folder)
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

# Ensure sample order matches
vst_mat <- vst_mat[, rownames(meta)]
condition <- factor(ifelse(meta$shortLetterCode == "TP", "Tumor", "Normal"))

# --------------------------------------------
# Step 2: Load Ferroptosis Genes
# --------------------------------------------
ferro_genes <- read_csv("ferroptosis_gene_list.csv")
gene_list <- ferro_genes$GeneSymbol

# Match ferroptosis genes with matrix
matched_genes <- intersect(gene_list, rownames(vst_mat))
vst_ferro <- vst_mat[matched_genes, , drop = FALSE]

# --------------------------------------------
# Step 3: Wilcoxon Test
# --------------------------------------------
wilcox_results <- lapply(matched_genes, function(gene) {
  exp_vals <- as.numeric(vst_ferro[gene, ])
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

# --------------------------------------------
# Step 4: Create Output Folder on Desktop
# --------------------------------------------
out_dir <- file.path("~/Desktop", "BRCA_Ferroptosis_Analysis")
dir.create(out_dir, showWarnings = FALSE)

# --------------------------------------------
# Step 5: Save Results
# --------------------------------------------
write.csv(wilcox_df, file.path(out_dir, "BRCA_Wilcoxon_Ferroptosis.csv"), row.names = FALSE)

logfc_cutoff <- 1
fdr_cutoff <- 0.05

sig_genes <- wilcox_df %>%
  filter(abs(LogFC) > logfc_cutoff & adj_p < fdr_cutoff)

write.csv(sig_genes, file.path(out_dir, "BRCA_Significant_Ferroptosis_Genes.csv"), row.names = FALSE)

# --------------------------------------------
# Step 6: Save Top 30 Genes + Heatmap
# --------------------------------------------
top_genes <- head(sig_genes$Gene, 30)
write.csv(data.frame(Gene = top_genes), file.path(out_dir, "BRCA_Top30_Ferroptosis_Genes.csv"), row.names = FALSE)

annotation_df <- data.frame(condition = condition)
rownames(annotation_df) <- colnames(vst_mat)

tiff(file.path(out_dir, "BRCA_Ferroptosis_Heatmap.tiff"), width = 1000, height = 800, res = 150)
pheatmap(
  vst_ferro[top_genes, , drop = FALSE],
  annotation_col = annotation_df,
  show_colnames = FALSE,
  scale = "row",
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()

# --------------------------------------------
# Final Message
# --------------------------------------------
cat("\n✅ Wilcoxon analysis complete.\n")
cat(paste("• Wilcoxon results:", file.path(out_dir, "BRCA_Wilcoxon_Ferroptosis.csv"), "\n"))
cat(paste("• Significant genes:", file.path(out_dir, "BRCA_Significant_Ferroptosis_Genes.csv"), "\n"))
cat(paste("• Top 30 genes:", file.path(out_dir, "BRCA_Top30_Ferroptosis_Genes.csv"), "\n"))
cat(paste("• Heatmap saved as:", file.path(out_dir, "BRCA_Ferroptosis_Heatmap.tiff"), "\n"))