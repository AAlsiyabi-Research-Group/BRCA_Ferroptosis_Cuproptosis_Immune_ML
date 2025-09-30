# --------------------------------------------
# BRCA_Immune_C7_Wilcoxon_Analysis.R
# Wilcoxon Test for C7 Immune Genes + Heatmap
# --------------------------------------------

# Load libraries
library(GSEABase)
library(readr)
library(dplyr)
library(pheatmap)

# 📂 Set Output Directory
main_dir <- "/Users/apple/Desktop/BRCA"
sub_dir <- "4.Wilcoxon_Immune_C7_BRCA"
output_dir <- file.path(main_dir, sub_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step 1: Load Normalized VST Matrix
vst_mat <- read.csv(file.path(main_dir, "BRCA_VST_Normalized_Matrix.csv"),
                    row.names = 1, check.names = FALSE)
meta <- read.csv("BRCA_Metadata_Final.csv", row.names = 1)
vst_mat <- vst_mat[, rownames(meta)]
condition <- factor(ifelse(meta$shortLetterCode == "TP", "Tumor", "Normal"))

# Step 2: Load Immune Gene Sets (C7)
gmt_file <- file.path(main_dir, "c7.all.v2025.1.Hs.symbols.gmt")
gmt <- getGmt(gmt_file)

# Extract gene sets
gene_sets <- lapply(gmt, geneIds)
names(gene_sets) <- sapply(gmt, function(gs) gs@setName)

# Filter gene sets with matches in vst matrix
matched_sets <- lapply(gene_sets, function(genes) intersect(genes, rownames(vst_mat)))
matched_sets <- Filter(function(x) length(x) > 1, matched_sets)  # keep only sets with >1 gene

# Initialize results
all_results <- list()

# Loop through gene sets
for (set_name in names(matched_sets)) {
  set_genes <- matched_sets[[set_name]]
  vst_subset <- vst_mat[set_genes, , drop = FALSE]
  
  wilcox_results <- lapply(set_genes, function(gene) {
    exp_vals <- as.numeric(vst_subset[gene, ])
    group1 <- exp_vals[condition == "Tumor"]
    group2 <- exp_vals[condition == "Normal"]
    test <- wilcox.test(group1, group2)
    logfc <- median(group1) - median(group2)
    data.frame(
      Gene = gene,
      GeneSet = set_name,
      Median_Tumor = median(group1),
      Median_Normal = median(group2),
      LogFC = logfc,
      p_value = test$p.value
    )
  })
  
  set_df <- bind_rows(wilcox_results)
  set_df$adj_p <- p.adjust(set_df$p_value, method = "fdr")
  all_results[[set_name]] <- set_df
}

# Combine all results
final_df <- bind_rows(all_results)
final_df <- final_df[order(final_df$adj_p), ]
write.csv(final_df, file.path(output_dir, "BRCA_Wilcoxon_ImmuneC7_AllResults.csv"), row.names = FALSE)

# Filter significant genes
logfc_cutoff <- 1
fdr_cutoff <- 0.05
sig_df <- final_df %>%
  filter(abs(LogFC) > logfc_cutoff & adj_p < fdr_cutoff)
write.csv(sig_df, file.path(output_dir, "BRCA_Significant_ImmuneC7_Genes.csv"), row.names = FALSE)

# ✅ Collapse repeated genes to keep only best (lowest adj_p) result per gene
collapsed_sig_df <- sig_df %>%
  group_by(Gene) %>%
  slice_min(order_by = adj_p, n = 1) %>%
  ungroup()

# Save the collapsed version
write.csv(collapsed_sig_df, file.path(output_dir, "BRCA_Significant_ImmuneC7_Genes_Collapsed.csv"), row.names = FALSE)

# ✅ Get top 30 unique genes by significance
top_genes <- head(collapsed_sig_df$Gene, 30)
write.csv(data.frame(Gene = top_genes), file.path(output_dir, "BRCA_Top30_ImmuneC7_Genes.csv"), row.names = FALSE)

# Heatmap
annotation_df <- data.frame(condition = condition)
rownames(annotation_df) <- colnames(vst_mat)

jpeg(file.path(output_dir, "BRCA_ImmuneC7_Heatmap_Top30.jpg"), width = 1000, height = 800, quality = 100)
pheatmap(
  vst_mat[top_genes, , drop = FALSE],
  annotation_col = annotation_df,
  show_colnames = FALSE,
  scale = "row",
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()

# ✅ Final message
cat("\n✅ C7 Immune Wilcoxon analysis complete.\n📁 Files saved in:", output_dir, "\n")
cat("• Full results: BRCA_Wilcoxon_ImmuneC7_AllResults.csv\n")
cat("• Significant genes (raw): BRCA_Significant_ImmuneC7_Genes.csv\n")
cat("• Significant genes (collapsed): BRCA_Significant_ImmuneC7_Genes_Collapsed.csv\n")
cat("• Top 30 genes (unique): BRCA_Top30_ImmuneC7_Genes.csv\n")
cat("• Heatmap: BRCA_ImmuneC7_Heatmap_Top30.jpg\n")