source("code/0_pkgs&func.R")

## Load Data ===================================================================
# [1] data read in 
gencode <- read.table("data/3_flair/gencode.v39/flair_gene.counts.tsv", 
                      row.names = 1, header = TRUE, check.names = FALSE)
pacbio <- read.table("data/3_flair/pacbio/flair_counts_filtered.tsv", 
                     row.names = 1, header = TRUE, check.names = FALSE)
metadata <- read.table("data/0_samples/flair_manifest.txt") %>% 
  dplyr::select(V1, V2) %>% 
  dplyr::rename(sample_id=V1, condition=V2) 
mapping <- read.csv("data/0_samples/mapping_table.csv", row.names = 1)
immune <- read.table("data/0_samples/hh_ova_immune.txt")
tmb <- read.csv("/Volumes/Wild Flower/1_CopyNumberOV/docs/00_mitéra/clinical/raw/ova/HH_TMB.csv",
                row.names = 1) %>% dplyr::rename(tmb=V1)

# [2] change long read id "-" to "_" 
colnames(gencode) <- gsub("-", "_", colnames(gencode))
colnames(pacbio) <- gsub("-", "_", colnames(pacbio))
metadata$sample_id <- gsub("-", "_", metadata$sample_id)
rownames(metadata) <- metadata$sample_id
mapping$long_read <- gsub("-", "_", mapping$long_read)

# [3] mapping to copy number sample ID 
mapping <- mapping %>% 
  dplyr::filter(long_read %in% colnames(gencode)) %>% # matched to expression matrix
  dplyr::filter(copy_number %in% rownames(immune)) # matched to immune 
  #dplyr::filter(copy_number %in% rownames(tmb))
gencode <- gencode[, mapping$long_read]
pacbio <- pacbio[, mapping$long_read]
metadata <- metadata[mapping$long_read, ]
immune <- immune[mapping$copy_number, ]
tmb <- tmb[mapping$copy_number, ] %>% as.data.frame()
rownames(tmb) <- mapping$copy_number
colnames(tmb) <- "tmb"

# [4] rename long reads to copy number ID 
colnames(gencode) <- mapping$copy_number
colnames(pacbio) <- mapping$copy_number
metadata$sample_id <- mapping$copy_number
rownames(metadata) <- metadata$sample_id

# [5] check: 
colnames(gencode) %>% intersect(colnames(pacbio)) %>% intersect(metadata$sample_id) %>% 
  intersect(rownames(immune)) %>% intersect(rownames(tmb))
# [1] "X1652" "X603"  "X1040" "X1922" "X120"  "X1188" "X1718" "X2212" "X137"  "X2429" "X2438" "X2556"
# [13] "X2554" "X2587" "X2846" "X3025" "X3250" "X3297"


## Data Preprocessing ==========================================================
# [1] tpm calculation 
if(!file.exists("data/6_association/matched_raw_df/gencode_normalized_variance_stabled_tpm.csv")){
  gencode_dds <- DESeqDataSetFromMatrix(countData = gencode, colData = metadata, design = ~1)
  gencode_dds <- gencode_dds[rowSums(counts(gencode_dds)) >= 10, ]
  gencode_dds <- estimateSizeFactors(gencode_dds)
  normalized_counts <- counts(gencode_dds, normalized = TRUE)
  gencode_vsd <- vst(gencode_dds, blind = FALSE)
  gencode_counts <- assay(gencode_vsd)
  gencode_counts <- as.data.frame(gencode_counts)
  write.csv(gencode_counts, "data/6_association/matched_raw_df/gencode_normalized_variance_stabled_tpm.csv")
}
if(!file.exists("data/6_association/matched_raw_df/pacbio_normalized_variance_stabled_tpm.csv")){
  pacbio_dds <- DESeqDataSetFromMatrix(countData = pacbio, colData = metadata, design = ~1)
  pacbio_dds <- pacbio_dds[rowSums(counts(pacbio_dds)) >= 10, ]
  pacbio_dds <- estimateSizeFactors(pacbio_dds)
  normalized_counts <- counts(pacbio_dds, normalized = TRUE)
  pacbio_vsd <- vst(pacbio_dds, blind = FALSE)
  pacbio_counts <- assay(pacbio_vsd)
  pacbio_counts <- as.data.frame(pacbio_counts)
  write.csv(pacbio_counts, "data/6_association/matched_raw_df/pacbio_normalized_variance_stabled_tpm.csv")
}
pacbio_counts <- read.csv("data/6_association/matched_raw_df/pacbio_normalized_variance_stabled_tpm.csv", row.names = 1)
gencode_counts <- read.csv("data/6_association/matched_raw_df/gencode_normalized_variance_stabled_tpm.csv")

# [2] convert ensembl id to hugo symbols
if(!file.exists("data/6_association/matched_raw_df/gencode_normalized_variance_stabled_tpm_hgnc.csv")){
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  enst_ids <- rownames(gencode_counts)
  enst_2_hugo <- getBM(filters = "ensembl_transcript_id_version",
                       attributes = c("ensembl_transcript_id_version", "hgnc_symbol", "external_gene_name"),
                       values = enst_ids,
                       mart = ensembl)
  enst_2_hugo <- enst_2_hugo[!duplicated(enst_2_hugo$ensembl_transcript_id_version), ]
  enst_2_hugo <- enst_2_hugo[!(enst_2_hugo$external_gene_name == ""), ]
  gencode_counts$ENST_ID <- rownames(gencode_counts)
  gencode_counts <- merge(gencode_counts, enst_2_hugo, 
                          by.x = "ENST_ID", by.y = "ensembl_transcript_id_version")
  gencode_counts$ENST_ID <- NULL
  gencode_counts$hgnc_symbol <- NULL
  gencode_counts <- gencode_counts %>% group_by(external_gene_name) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>% as.data.frame()
  rownames(gencode_counts) <- gencode_counts$external_gene_name
  gencode_counts$external_gene_name <- NULL
  write.csv(gencode_counts, "data/6_association/matched_raw_df/gencode_normalized_variance_stabled_tpm_hgnc.csv")
}
if(!file.exists("")){
  files <- list.files(path = "data/3_salmon/gencode.v39/output", pattern = ".sf", full.names = TRUE, 
                      recursive = TRUE)
  names(files) <- stringr::str_split(files, pattern = "/", simplify = TRUE)[, 5] %>%
    stringr::str_replace("_quant", "")
  txi.salmon <- tximport(files, type = "salmon", txIn=TRUE, txOut = TRUE)
  gencode_SR <- txi.salmon$abundance %>% as.data.frame()
  gencode_SR$combined <- rownames(gencode_SR)
  rownames(gencode_SR) <- NULL
  gencode_SR <- gencode_SR %>% mutate(gene_symbol = sapply(strsplit(as.character(combined), "\\|"), "[", 2))
  gencode_SR$combined <- NULL
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  enst_ids <- gencode_SR$gene_symbol
  enst_2_hugo <- getBM(filters = "ensembl_transcript_id_version",
                       attributes = c("ensembl_gene_id_version", "external_gene_name"),
                       values = enst_ids,
                       mart = ensembl)
  enst_2_hugo <- enst_2_hugo[!duplicated(enst_2_hugo$ensembl_transcript_id_version), ]
  enst_2_hugo <- enst_2_hugo[!(enst_2_hugo$external_gene_name == ""), ]
}

rna_seq <- read.csv("data/6_association/matched_raw_df/gencode_normalized_variance_stabled_tpm_hgnc.csv", 
                    row.names = 1)

## Estimate Immune =============================================================
# [1] MCP counter 
mcp_results <- MCPcounter.estimate(expression = rna_seq, featuresType = "HUGO_symbols") 
mcp_results <- as.data.frame(t(mcp_results))

# [2] TLS
tls_signature = c("SKAP1", "PTGDS", "LAMP3")
#tls_signature = c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", 
#                  "PTGDS", "CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3")
tls_rna_seq <- rna_seq[tls_signature, ] %>% na.omit()
tls_rna_seq <- as.matrix(t(tls_rna_seq)) 
tls_rna_seq <- scale(tls_rna_seq, center = F) 
tls_rna_seq <- as.data.frame(tls_rna_seq) 
tls_rna_seq$TLS_mean <- apply(tls_rna_seq, 1, mean) 
tls_rna_seq$Row.names <- rownames(tls_rna_seq)

# [3] HLA gene expressions 
immune_genes <- data.frame(MHC1=colSums(rna_seq[c("HLA-A", "HLA-B", "HLA-C"), ]), 
                           MHC2=colSums(rna_seq[c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DRA"), ]))

## Association Study ===========================================================
# [1] merge all the things together 
combined_immune <- merge(immune, tmb, by="row.names")
combined_immune <- merge(combined_immune, mcp_results, by.x="Row.names", by.y="row.names")
combined_immune <- merge(combined_immune, tls_rna_seq, by="Row.names")
combined_immune <- merge(combined_immune, immune_genes, by.x="Row.names", by.y="row.names")
rownames(combined_immune) <- combined_immune$Row.names
combined_immune$Row.names <- NULL
combined_immune <- combined_immune[colnames(pacbio_counts), ] ## align 
combined_immune$SKAP1 <- NULL
combined_immune$PTGDS <- NULL
combined_immune$LAMP3 <- NULL

# [2] create empty dataframes 
immune_association <- matrix(data = NA, nrow = nrow(pacbio_counts), 
                             ncol = ncol(combined_immune)) %>% as.data.frame()
rownames(immune_association) <- rownames(pacbio_counts)
colnames(immune_association) <- colnames(combined_immune)

# [3] association analysis 
for(i in colnames(combined_immune)){
  immune_exp <- combined_immune[, i]
  for(j in rownames(pacbio_counts)){
    isoform_exp <- pacbio_counts[j, ] %>% unlist()
    test <- cor.test(immune_exp, isoform_exp, method = "spearman")
    immune_association[j, i] <- test$p.value
  }
}
immune_association_adjusted <- immune_association %>% 
  mutate(across(everything(), ~ p.adjust(., method = "fdr")))

# [4] view the results
check_above <- function(column) {all(column < 0.10)}
sapply(immune_association_adjusted, check_above)
view <- immune_association_adjusted %>% filter(rowSums(. < 0.10) > 0)



save.image("data/6_association/immune.RData")





