# -------------------------------------
# BRCA_deseq2_DEG_analysis.R
# Step 2: Differential Expression using DESeq2
# -------------------------------------

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "apeglm", "pheatmap", "biomaRt", "ggplot2", "grid"),
                     ask = FALSE, force = TRUE)

library(DESeq2)
library(apeglm)
library(pheatmap)
library(biomaRt)
library(ggplot2)
library(grid)   # for editing grobs in pheatmap

# -------------------------------------
# Load Data
# -------------------------------------
counts <- read.csv("/Users/apple/Desktop/BRCA/TCGA_BRCA_STAR_Counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("/Users/apple/Desktop/BRCA/BRCA_Metadata_Final.csv", row.names = 1)
counts <- counts[, rownames(meta)]  # match sample order

# -------------------------------------
# Filter Low-Count Genes
# -------------------------------------
keep <- rowSums(counts >= 10) >= 0.1 * ncol(counts)
counts_filtered <- counts[keep, ]

# -------------------------------------
# Metadata and Design
# -------------------------------------
if (!"shortLetterCode" %in% colnames(meta)) {
  stop("Column 'shortLetterCode' not found in metadata.")
}
condition <- factor(ifelse(meta$shortLetterCode == "TP", "Tumor", "Normal"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = data.frame(condition = condition),
  design    = ~ condition
)

dds <- DESeq(dds)
res <- lfcShrink(dds, coef = 2, type = "apeglm")
res_ordered <- res[order(res$padj), ]

# -------------------------------------
# Map Ensembl IDs to Gene Symbols
# -------------------------------------
ensembl_ids <- gsub("\\..*", "", rownames(res_ordered))
res_df <- as.data.frame(res_ordered)
res_df$ensembl_gene_id <- ensembl_ids

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = mart
)

res_merged <- merge(res_df, gene_map, by = "ensembl_gene_id", all.x = TRUE)

res_merged$Regulation <- ifelse(
  res_merged$padj < 0.05 & res_merged$log2FoldChange > 1,  "Up regulated",
  ifelse(res_merged$padj < 0.05 & res_merged$log2FoldChange < -1, "Down regulated", "Not significant")
)

res_final <- res_merged[, c("hgnc_symbol", "ensembl_gene_id",
                            setdiff(colnames(res_merged), c("hgnc_symbol", "ensembl_gene_id")))]
deg <- subset(res_final, Regulation %in% c("Up regulated", "Down regulated"))

# -------------------------------------
# Output Directory Setup
# -------------------------------------
output_dir <- "/Users/apple/Desktop/BRCA/BRCA_deseq2_DEG_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------
# Save Results
# -------------------------------------
write.csv(res_final, file.path(output_dir, "DESeq2_BRCA_all_results.csv"), row.names = FALSE)
write.csv(deg,       file.path(output_dir, "DESeq2_BRCA_DEGs_filtered.csv"), row.names = FALSE)

# -------------------------------------
# Heatmap of Top 50 DEGs (High-res JPG, serif + bold)
# -------------------------------------
vsd <- vst(dds, blind = FALSE)

# First 50 DEGs by adjusted p-value from 'deg' table
top50_ids <- head(deg$ensembl_gene_id, 50)

vst_ensembl_ids <- gsub("\\..*", "", rownames(assay(vsd)))
matching_rows   <- which(vst_ensembl_ids %in% top50_ids)

gene_symbols <- deg$hgnc_symbol[match(top50_ids, deg$ensembl_gene_id)]
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- top50_ids[is.na(gene_symbols) | gene_symbols == ""]

heatmap_data <- assay(vsd)[matching_rows, ]
rownames(heatmap_data) <- gene_symbols

annotation_df <- data.frame(condition = condition)
rownames(annotation_df) <- colnames(assay(vsd))

# Build pheatmap silently so we can edit fonts to bold serif
ph <- pheatmap(
  heatmap_data,
  show_rownames  = TRUE,
  show_colnames  = FALSE,
  annotation_col = annotation_df,
  cluster_cols   = FALSE,
  scale          = "row",
  color          = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize       = 12,
  fontsize_row   = 10,
  fontsize_col   = 10,
  silent         = TRUE
)

# --- Helper: set all text grobs to serif + bold safely (avoid 'font'/'fontface' clash)
make_bold_serif <- function(g) {
  if (!is.null(g$gp)) {
    # If a numeric 'font' is present, drop it to avoid conflict with 'fontface'
    if (!is.null(g$gp$font)) g$gp$font <- NULL
    g$gp$fontfamily <- "serif"
    g$gp$fontface   <- "bold"
  }
  if (!is.null(g$children)) {
    g$children <- lapply(g$children, make_bold_serif)
  }
  g
}

# Apply to pheatmap gtable grobs with a safety net
ph$gtable$grobs <- lapply(ph$gtable$grobs, function(gr) {
  tryCatch(make_bold_serif(gr), error = function(e) gr)
})

# Save high-res JPG (BioRender-friendly)
jpeg(file.path(output_dir, "BRCA_DEG_heatmap.jpg"),
     width = 3000, height = 2400, res = 600, quality = 100)
grid.newpage()
grid.draw(ph$gtable)
dev.off()

# -------------------------------------
# Volcano Plot (High-res JPG, serif + bold)
# -------------------------------------
volcano_dir <- file.path(output_dir, "volcano_plot")
dir.create(volcano_dir, showWarnings = FALSE)

# Clean padj
res_final <- res_final[!is.na(res_final$padj), ]
res_final$padj[res_final$padj == 0] <- .Machine$double.xmin

# Assign colors
res_final$color <- ifelse(
  res_final$padj < 0.05 & res_final$log2FoldChange > 1,  "Up",
  ifelse(res_final$padj < 0.05 & res_final$log2FoldChange < -1, "Down", "NS")
)

# Plot
jpeg(file.path(volcano_dir, "BRCA_volcano_plot.jpg"),
     width = 3000, height = 2500, res = 600, quality = 100)
ggplot(res_final, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Up" = "orange", "Down" = "green", "NS" = "black")) +
  theme_minimal(base_size = 16, base_family = "serif") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano Plot - BRCA DESeq2", x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(face = "bold")
  )
dev.off()

cat("\n✅ DESeq2 BRCA analysis complete.",
    "\n• All results: 'DESeq2_BRCA_all_results.csv'",
    "\n• Filtered DEGs: 'DESeq2_BRCA_DEGs_filtered.csv'",
    "\n• Heatmap (JPG): 'BRCA_DEG_heatmap.jpg'",
    "\n• Volcano (JPG): 'volcano_plot/BRCA_volcano_plot.jpg'",
    "\nSaved to:", output_dir, "\n")