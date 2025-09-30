# ----------------------------------------
# Final Ferroptosis Gene Extraction Script (FerrDb v2 - Symbol column)
# ----------------------------------------

library(dplyr)
library(readr)

# Set subfolder path (relative to LUAD-project)
data_path <- "ferroptosis_database_files"

# Function to read each FerrDb file and label with role
read_and_label <- function(file, role) {
  df <- read_csv(file.path(data_path, file), show_col_types = FALSE)
  if (!"Symbol" %in% colnames(df)) {
    stop(paste("❌ Column 'Symbol' not found in", file))
  }
  df %>%
    select(GeneSymbol = Symbol) %>%
    filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
    mutate(Role = role)
}

# Load and tag genes from all regulatory categories
driver_df       <- read_and_label("driver.csv", "Driver")
suppressor_df   <- read_and_label("suppressor.csv", "Suppressor")
marker_df       <- read_and_label("marker.csv", "Marker")
unclassified_df <- read_and_label("unclassified.reg.csv", "Unclassified")

# Combine, deduplicate, and collapse multiple roles
ferroptosis_genes <- bind_rows(driver_df, suppressor_df, marker_df, unclassified_df) %>%
  group_by(GeneSymbol) %>%
  summarise(Role = paste(sort(unique(Role)), collapse = ", "), .groups = "drop") %>%
  arrange(GeneSymbol)

# Save to LUAD-project directory
write_csv(ferroptosis_genes, file.path("ferroptosis_gene_list.csv"))

# Confirmation message
cat("✅ Ferroptosis gene list saved to 'ferroptosis_gene_list.csv' in LUAD-project.\n")
cat("🧬 Total unique genes:", nrow(ferroptosis_genes), "\n")