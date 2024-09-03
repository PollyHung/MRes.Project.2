# source the packages 
source("~/Desktop/tappAS/codes/1_pkgs&func.R")


## experiment design =========================================================== 
# read in design and change names 
experiment_design <- read.delim("docs/sample_id/flair_manifest.txt", header = FALSE) %>% 
  dplyr::select(V1, V2) %>% 
  dplyr::rename(sample = V1, group = V2) 

# rename the sample id to remove - and replace it by _
experiment_design$sample <- gsub("-", "_", experiment_design$sample)

# arrange the design matrix to put all normal samples on top and tumour samples on bottom 
experiment_design <- experiment_design %>% dplyr::arrange(group)

# add another column of sample_id for ease of data calling 
experiment_design$sample_id <- c(paste0("normal_", c(1:6)), paste0("tumour_", c(1:19)))




## expression matrix ===========================================================
# read in data frame and sample lsit 
expression_matrix <- read.table(file = "docs/flair/flair_counts_filtered.tsv", 
                                header = TRUE, check.names = FALSE, row.names = 1)
sample_list <- read.table("docs/sample_id/sample_list.txt") %>% unlist %>% unname

# like experiment design matrix, rename the columns by replacing - with _ 
colnames(expression_matrix) <- gsub("-", "_", sample_list)

# turn it into a matrix 
expression_matrix <- as.matrix(expression_matrix)

# reorder the column of the expression matrix in the sequence of experiment matrix 
expression_matrix <- expression_matrix[, experiment_design$sample]
colnames(expression_matrix) == experiment_design$sample

# turn it into a TPM normalized matrix 
library(edgeR)
calculateTPM <- function(counts, lengths) {
  rpk <- counts / lengths
  tpm <- rpk / rowSums(rpk) * 1e6
  return(tpm)
}
mylength <- classifications$length
names(mylength) <- classifications$isoform
mylength <- mylength[rownames(expression_matrix)]
tpm <- calculateTPM(expression_matrix, mylength)
write.table(tpm, file = "docs/processed/tpm_normalized.txt", sep = "\t", quote = FALSE)

tpm <- read.table("~/Desktop/tappAS/docs/processed/tpm_normalized.txt")


## annotation_features =========================================================
annotation_features <- readGFF("docs/sqanti/qc.gff3")
gtf <- readGFF("docs/sqanti/filter_rule.filtered.gtf")
classifications <- read.delim("docs/sqanti/filter_rule_RulesFilter_result_classification.txt")

PB_ids <- classifications %>% dplyr::filter(filter_result == "Isoform") %>% 
  pull(isoform) %>% unique()
write.table(PB_ids, "docs/sample_id/PB.ID.txt", quote = F, sep = "\t", row.names = F, col.names = F)



## list all files and save 
normalize_files <- ls()
save.image("docs/processed/normalized_data.RData")













