library(VariantAnnotation)
library(dplyr)
library(magrittr)
library(tidyverse)
library(vcfR)

setwd("/Volumes/Wild_Flower/OV_SpliceVariants/")
sample_list <- read.csv("data/0_samples/mapping_table.csv", row.names = 1)
wes_list <- list.files("data/4_mutect2/output/")
id <- intersect(sample_list$copy_number, wes_list)

mutect2_list <- list()
for(sample_id in id){
  mutect2_path = file.path("data/4_mutect2/output", sample_id, "annovar_annot.hg38_multianno.txt")
  mutect2 <- read.table(mutect2_path, sep = "\t", header = T)
  mutect2$sample <- sample_id
  mutect2_list[[sample_id]] <- mutect2
}

mutations <- do.call(rbind, mutect2_list)
rownames(mutations) <- NULL
save(list = "mutations", file = "data/4_mutect2/mutation.Rd")


map <- read.csv("/Volumes/Wild_Flower/OV_SpliceVariants/data/0_samples/mapping_table.csv", row.names = 1) %>% 
  dplyr::mutate(long_read = gsub("_", "-", long_read)) %>% 
  dplyr::filter(long_read %in% colnames(psi_lr.imp))
load("/Volumes/Wild_Flower/OV_SpliceVariants/data/4_mutect2/mutation.Rd")

gene_list <- c("TP53", "BRCA1", "BRCA2", "CCNE1", "CCND1", "PIK3CA", "RB1", "MYC", 
               "NF1", "CDKN2A", "NOTCH3", "BRD4", "ARID1A", "CDK12", "MUC16", "RAD51", 
               "KRAS", "PTEN", "CSMD3")  

mutation_list <- list()     

for(gene in gene_list) {
  cat("analysing:", gene)
  wes_id <- mutations$sample %>% unique()
  sample_mut <- map$long_read[map$copy_number %in% unique(mutations$sample[mutations$Gene.ensGene == gene])] 
  sample_wt <- map$long_read[map$copy_number %in% setdiff(wes_id, unique(mutations$sample[mutations$Gene.ensGene == gene]))] 
  
  psi_mut <- psi_lr.imp[, sample_mut, drop = FALSE]
  psi_wt <- psi_lr.imp[, sample_wt, drop = FALSE]
  
  tryCatch({
    if (ncol(psi_mut) == 0 || ncol(psi_wt) == 0) {
      stop("No data in psi_mut or psi_wt, skipping this gene")
    }
    
    mutation_dpsi <- data.frame(
      deltaPSI = rowMeans(psi_mut) - rowMeans(psi_wt), 
      absDeltaPSI = abs(rowMeans(psi_mut) - rowMeans(psi_wt)), 
      p_value = rep(NA, nrow(psi_mut))
    )
    
    for(i in rownames(mutation_dpsi)) {
      test <- wilcox.test(psi_mut[i, ], psi_wt[i, ], exact = FALSE)
      mutation_dpsi[i, "p_value"] <-  test$p.value
    }
    
    mutation_dpsi <- mutation_dpsi %>% dplyr::filter(absDeltaPSI > 0.20)
    mutation_dpsi$p_adj <- p.adjust(mutation_dpsi$p_value, method = "fdr")
    mutation_dpsi <- mutation_dpsi %>% dplyr::filter(p_adj < 0.25)
    mutation_dpsi$AS_event <- rownames(mutation_dpsi)
    mutation_dpsi$gene <- gene
    
    mutation_list[[gene]] <- mutation_dpsi
    
  }, error = function(e) {
    message(paste("Skipping gene", gene, "due to error:", e$message))
  })
}

## MUC16 skipped because all samples have MUC16 mutation 
mutation_psi <- do.call(rbind, mutation_list)

save(list = c("mutation_psi", "mutations"), file = "data/4_mutect2/mutation.Rd")














# setwd("/Volumes/Wild_Flower/OV_SpliceVariants/data/")
# 
# load("6_association/sample_data/psi_knn_imp.Rd")
# genes <- c("TP53", "PIK3CA", "RB1", "KRAS", "PIK3R1", # key 
#            "RNF43", "FGFR2", "ERBB2", "MET", "AKT1", "PTEN", "GNAS", "FGFR1", "FGFR3", "SMARCA4", "VHL", "SPOP") # other 
# mapping <- read.table("1_sampleSheet/mapping_complete.txt", header = T, sep = "\t")
# sample_list <- intersect(list.files("4_mutect2/output/"), mapping$copy_number)
# mapping <- mapping %>% dplyr::filter(copy_number %in% sample_list)
# lr_list <- intersect(mapping$long_read, colnames(psi_lr.imp))
# mapping <- mapping %>% dplyr::filter(long_read %in% lr_list)
# 
# mutect2 <- list()
# 
# for(sample_id in sample_list){
#   vcf_file <- file.path("4_mutect2/output", sample_id, paste0(sample_id, ".filtered.somatic.vcf.gz"))
#   vcf <- readVcf(vcf_file, "hg38")
#   vcf_info <- as.data.frame(info(vcf)$FUNCOTATION)
#   vcf_info$value <- iconv(vcf_info$value, from = "UTF-8", to = "UTF-8", sub = "")
#   temp2 <- vcf_info %>%
#     separate(value, into = c("gene_name", "build", "chr", "start", "end", "mutation_type", "mutation_detail", "specific_mutation", "rest"), sep = "\\|", fill = "right", extra = "merge", remove = FALSE) %>%
#     mutate(gene_name = gsub("\\[", "", gene_name))
#   temp2$rest <- sub(".*(g\\.chr)", "\\1", temp2$rest)
#   gene_info <- temp2 %>% 
#     separate(rest, into = c("mut_event", "ensemble_transcript_id"), sep = "\\|", extra = "drop")
#   AS_info <- data.frame(AS_FilterStatus = unlist(info(vcf)$AS_FilterStatus), 
#                         FUNCOTATION = unlist(info(vcf)$FUNCOTATION), 
#                         NLOD = unlist(info(vcf)$NLOD), 
#                         TLOD = unlist(info(vcf)$TLOD))
#   mutation <- merge(AS_info, gene_info, by.x = "FUNCOTATION", by.y = "value")
#   mutation <- mutation %>% dplyr::select(gene_name, ensemble_transcript_id, 
#                                          chr, start, end, mutation_type, 
#                                          mutation_detail, specific_mutation, mut_event, 
#                                          AS_FilterStatus, NLOD, TLOD)
#   mutation_filter <- mutation %>% dplyr::filter(AS_FilterStatus == ".")
#   mutation_filter$sample_id <- sample_id
#   mutect2[[sample_id]] <- mutation_filter
# }
# 
# mutect2_combined <- do.call(rbind, mutect2)
# rownames(mutect2_combined) <- NULL
# mutect2_combined <- merge(mutect2_combined, mapping, by.x = "sample_id", by.y = "copy_number")
# 
# save(list = "mutect2_combined", file = "4_mutect2/mutect2_combined.Rd")
# 
# 
# vcf_file <- file.path("4_mutect2/output", sample_id, paste0(sample_id, ".filtered.somatic.vcf.gz"))




