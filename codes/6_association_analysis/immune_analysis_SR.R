source("code/0_pkgs&func.R")

## Load Files ==================================================================
# [1] salmon 
files <- list.files(path = "data/3_salmon/output", pattern = ".sf", full.names = TRUE, 
                    recursive = TRUE)
names(files) <- stringr::str_split(files, pattern = "/", simplify = TRUE)[,4] %>%
  stringr::str_replace("_quant", "")
txi.salmon <- tximport(files, type = "salmon", txIn=TRUE, txOut = TRUE)

# [2] sqanti 
classification <- read.delim("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt") %>% 
  dplyr::filter(filter_result == "Isoform")
isoform_id <- intersect(rownames(txi.salmon$counts), classification$isoform)
count_matrix <- txi.salmon$abundance[isoform_id, ]

# [3] patient sheet 
survival <- read.delim("data/0_samples/hh_ova_surv.txt")
clinical <- read.delim("data/0_samples/hh_ova.txt")
immune <- read.delim("data/0_samples/hh_ova_immune.txt")
mapping <- read.csv("data/0_samples/mapping.csv") %>% 
  dplyr::rename("Row.names" = "ICBRC", "SR.id" = "X.1") %>% 
  dplyr::select(Row.names, SR.id) %>% 
  dplyr::filter(SR.id != "")
metadata <- merge(immune, clinical, by="row.names") 
metadata <- left_join(mapping, metadata, by="Row.names")
rownames(metadata) <- metadata$SR.id

# [4] alignment 
SR_id <- intersect(metadata$SR.id, colnames(count_matrix))
count_matrix <- count_matrix[, SR_id]
metadata <- metadata[SR_id, ]

## Immune Analysis ===========================================================
final_counts <- count_matrix[]
immune_response <- matrix(NA, nrow = nrow(final_counts), ncol = 11) %>% as.data.frame()
colnames(immune_response) <- colnames(metadata[3:10])
rownames(immune_response) <- rownames(final_counts)

# Function to perform correlation calculation for each isoform expression
lm_by_row <- function(immune_fac, return_item = "p_val") {
  # Create temp dataframe
  temp <- matrix(NA, nrow = nrow(final_counts), ncol = 2) %>% as.data.frame()
  rownames(temp) <- rownames(final_counts)
  colnames(temp) <- c("p_val", "p_adj")
  
  # Get immune column
  immune_col <- metadata[, immune_fac]
  names(immune_col) <- rownames(metadata)
  
  # Loop over each isoform
  for (isoform_id in rownames(final_counts)) {
    tryCatch({
      cor_model <- cor.test(final_counts[isoform_id, ], immune_col, method = "spearman", use = "complete.obs")
      temp[isoform_id, "p_val"] <- cor_model$p.value
    }, error = function(e) {
      cat("Error with isoform", isoform_id, ":", conditionMessage(e), "\n")
    })
  }
  
  # Adjust p-values using FDR
  if (return_item == "p_adj") {
    temp$p_adj <- p.adjust(temp$p_val, method = "fdr")
    return(temp$p_adj)
  } else if (return_item == "p_val") {
    return(temp$p_val)
  } else {
    stop("Invalid return_item specified. Use 'p_val' or 'p_adj'.")
  }
}

# Apply the function to each relevant column in the metadata
for (immune_fac in colnames(metadata[3:10])) {
  immune_response[, immune_fac] <- lm_by_row(immune_fac, return_item = "p_adj")
}















