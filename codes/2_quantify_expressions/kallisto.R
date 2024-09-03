# Load necessary library
library(dplyr)

# Define the path to Kallisto output directories
expression_path <- "~/MRes.project.2/docs/0_preprocessing/short_reads/transcription_expression"
files <- list.files(path = expression_path, pattern = "abundance.tsv", full.names = TRUE, recursive = TRUE)

expression_data <- list()

# Read each Kallisto abundance file and store the TPM values
for (file in files) {
  df <- read.table(file, header = TRUE, sep = "\t")
  expression_data[[sample_name]] <- df %>% select(target_id, tpm) %>% rename(!!sample_name := tpm)
}

# Merge all data frames by 'target_id'
expression_df <- Reduce(function(x, y) merge(x, y, by = "target_id", all = TRUE), expression_data)

# Rename 'target_id' to 'ID'
colnames(expression_df)[1] <- "ID"

# Save the DataFrame as a tab-separated file
write.table(expression_df, file = "expression_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
