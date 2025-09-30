# -------------------------------
# VST Normalization for LUAD with Gene Symbol Mapping
# -------------------------------

library(DESeq2)
library(readr)
library(biomaRt)
library(dplyr)

# Load count matrix (STAR output) and metadata
counts <- read.csv("/Users/apple/Desktop/BRCA/TCGA_BRCA_STAR_Counts.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("/Users/apple/Desktop/BRCA/BRCA_Metadata_Final.csv", row.names = 1)

# Ensure columns in counts match metadata rows
counts <- counts[, rownames(meta)]
condition <- factor(ifelse(meta$shortLetterCode == "TP", "Tumor", "Normal"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(condition = condition),
  design = ~ condition
)

# Estimate size factors and apply VST
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)

# Extract normalized matrix
vst_mat <- assay(vsd)

# Map Ensembl IDs to HGNC gene symbols
ensembl_ids <- gsub("\\..*", "", rownames(vst_mat))  # Remove Ensembl version suffix

# Use biomaRt to get gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Keep only rows with symbols
annot <- annot %>% filter(hgnc_symbol != "")

# Add mapping to matrix
vst_df <- as.data.frame(vst_mat)
vst_df$ensembl_gene_id <- ensembl_ids
vst_merged <- inner_join(annot, vst_df, by = "ensembl_gene_id")

# Collapse duplicates by averaging (if same symbol appears multiple times)
vst_symbol_avg <- vst_merged %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

# Set rownames and clean
vst_final <- as.data.frame(vst_symbol_avg)
rownames(vst_final) <- vst_final$hgnc_symbol
vst_final$hgnc_symbol <- NULL

# Save normalized matrix with gene symbols
write.csv(vst_final, "BRCA_VST_Normalized_Matrix.csv")

cat("✅ VST normalization complete with gene symbol mapping.\n")
cat("📄 File saved: BRCA_VST_Normalized_Matrix.csv\n")