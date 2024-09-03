source("code/1_pkgs&func.R")

## Load in data ================================================================
# [1] load in tables 
counts <- read.table("data/3_flair/pacbio/flair.counts.tsv", row.names = 1, 
                     header = TRUE, check.names = FALSE)
copy_number <- read.table("/Volumes/Wild Flower/1_CopyNumberOV/docs/HH_ova/gistic/facets_seg/cval_50/hisens_adj/broad_values_by_arm.txt", 
                          sep = "\t", header = TRUE, row.names = 1)
mapping <- read.csv("data/0_samples/mapping_table.csv", row.names = 1) 

# [2] rename 
colnames(counts) <- gsub("-", "_", colnames(counts))
mapping$long_read <- gsub("-", "_", mapping$long_read)

# [3] mapping 
mapping <- mapping %>% dplyr::filter(long_read %in% colnames(counts)) 
mapping <- mapping[which(mapping$copy_number %in% colnames(copy_number)), ]

# [4] select 
counts <- counts[, mapping$long_read]
copy_number <- copy_number[, mapping$copy_number]

# [5] tpm normalization 
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = mapping, 
                              design = ~1)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = FALSE)
transformed_counts <- assay(vsd)



## Association ================================================================= 
association_matrix <- matrix(data = NA, 
                             nrow = nrow(transformed_counts), 
                             ncol = nrow(copy_number)) %>% as.data.frame()
rownames(association_matrix) <- rownames(transformed_counts)
colnames(association_matrix) <- rownames(copy_number)

for(arm in rownames(copy_number)){
  arm_level_cn <- unlist(copy_number[arm, ])
  for(isoform in rownames(transformed_counts)){
    isoform_expression <- unlist(transformed_counts[isoform, ])
    test <- cor.test(isoform_expression, arm_level_cn)
    association_matrix[isoform, arm] <- test$p.value
  }
}

association_adjusted <- association_matrix %>%
  mutate(across(everything(), ~ p.adjust(., method = "fdr")))
write.csv(association_adjusted, "data/6_association/copy_number_associated_padj.csv")
write.csv(association_matrix, "data/6_association/copy_number_associated_pval.csv")

association_filtered <- association_adjusted %>% filter(rowSums(. < 0.10) > 0)

check_above <- function(column) { all(column > 0.20) }
sapply(association_adjusted, check_above)

save.image("data/6_association/copy_number.RData")




