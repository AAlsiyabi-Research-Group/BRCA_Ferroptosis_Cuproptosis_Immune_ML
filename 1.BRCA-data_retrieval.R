# TCGA_BRCA_chunked_download.R (Robust with retry-safe download)

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"), ask = FALSE, force = TRUE)

library(TCGAbiolinks)
library(SummarizedExperiment)

# Set up working directory
project_dir <- "~/Desktop/BRCA-project"
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)

# Step 1: Query metadata (do not download yet)
query_all <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# Step 2: Extract barcodes and split into ~100 sample chunks
barcodes <- getResults(query_all, cols = "cases")
chunks <- split(barcodes, ceiling(seq_along(barcodes) / 100))

# Step 3: Download and prepare each chunk
output_dir <- file.path(project_dir, "BRCA_chunks")
dir.create(output_dir, showWarnings = FALSE)
all_data <- list()

for (i in seq_along(chunks)) {
  cat("\n📦 Processing chunk", i, "of", length(chunks), "\n")
  
  query_chunk <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    barcode = chunks[[i]]
  )
  
  rds_file <- file.path(output_dir, paste0("chunk_", i, ".rds"))
  
  if (!file.exists(rds_file)) {
    tryCatch({
      GDCdownload(query_chunk)
      data_chunk <- GDCprepare(query_chunk)
      saveRDS(data_chunk, rds_file)
      all_data[[i]] <- data_chunk
    }, error = function(e) {
      cat("❌ Error in chunk", i, ":", conditionMessage(e), "\n")
      all_data[[i]] <- NULL
    })
  } else {
    all_data[[i]] <- readRDS(rds_file)
    cat("✅ Chunk", i, "already exists and loaded.\n")
  }
}

# Step 4: Filter out failed (NULL) chunks
all_data <- Filter(Negate(is.null), all_data)

# Step 5: Align metadata across successful chunks
cat("\n✅ Aligning metadata across chunks...\n")
common_cols <- Reduce(intersect, lapply(all_data, function(se) colnames(colData(se))))
all_data_trimmed <- lapply(all_data, function(se) {
  colData(se) <- colData(se)[, common_cols]
  return(se)
})

# Step 6: Combine all into one SummarizedExperiment
cat("\n✅ Combining aligned chunks...\n")
combined_data <- do.call(cbind, all_data_trimmed)

# Step 7: Save expression matrix and metadata
expr_matrix <- as.data.frame(assay(combined_data))
write.csv(expr_matrix, file.path(project_dir, "TCGA_BRCA_STAR_Counts.csv"))

meta_df <- as.data.frame(colData(combined_data))
meta_df[] <- lapply(meta_df, function(x) if (is.list(x)) sapply(x, toString) else x)
write.csv(meta_df, file.path(project_dir, "TCGA_BRCA_Metadata.csv"))

saveRDS(combined_data, file.path(project_dir, "TCGA_BRCA_SummarizedExperiment.rds"))

cat("\n🎉 Done! All output saved to:", project_dir, "\n")