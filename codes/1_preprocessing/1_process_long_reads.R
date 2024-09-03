## This codes processes long reads transcriptome matrix 
## 1. Creates TPM matrix from count matrix 
## 2. Filter for low expressing isoforms 
## 3. Collapse isoforms into gene expression 


source("code/0_pkgs&func.R")


# read in raw count matrix 
counts <- read.delim("data/3_flair/flair_counts_filtered.tsv", row.names = 1, check.names = FALSE) 
counts$pb_id <- str_replace(rownames(counts), "\\.[^.]*$", "")

# read in classification file 
gff3 <- read.gff("data/2_sqanti/qc.gff3")
gff3 <- gff3 %>% dplyr::filter(type == "gene")
gff3$gene <- str_extract(gff3$attributes, "(?<=ID=)[^;]+")
gff3$pb_id <- str_replace(gff3$seqid, "\\.[^.]*$", "")

# mapping table? 
mapping <- gff3 %>% dplyr::select(pb_id, gene)
collapsed_mapping <- mapping %>%
  group_by(pb_id) %>%
  summarize(most_common_gene = names(sort(table(gene), decreasing = TRUE))[1])

# merge 
counts <- merge(counts, collapsed_mapping, by="pb_id")
counts$pb_id <- NULL
collapsed_counts <- counts %>%
  group_by(most_common_gene) %>%
  summarise(across(everything(), sum))
collapsed_counts <- as.data.frame(collapsed_counts)
rownames(collapsed_counts) <- collapsed_counts$most_common_gene
collapsed_counts$most_common_gene <- NULL

# MCP counter genes, are they in the collapsed dataframe? 
MCP_genes <- read.delim("references/MCPcounter_genes.txt")
which(!MCP_genes$HUGO.symbols %in% rownames(collapsed_counts)) ## 38 in but 73 not in 
 
# use MCP counter regardless
mcp_results <- MCPcounter.estimate(
  expression = collapsed_counts, ## dataframe rna_seq 
  featuresType = "HUGO_symbols")
mcp_results <- as.data.frame(t(mcp_results))
mcp_results <- select(mcp_results, all_of(mcp_pat))



#4. Specific genes: c("CD274", "GZMB", "GNLY", "PRF1", "MHCI"). [MHCI in human](https://www.ncbi.nlm.nih.gov/books/NBK27156/#_A579_) is HLA-A, HLA-B, HLA-C.   

#**Functions:**     
#```{r immune_pheno}
## Separate out the RNA_seq from TCGA PANCAN rna-seq dataset 
immune_pheno <- function(cancer, rna_seq, mcp_pat, tls_sig, tmb_col, gene_pat, write_each = FALSE, write_uni = TRUE){
  ## cancer = cancer type, usually given as "i" in a for (i in TCGA){} loop 
  ## rna_seq = the rna_seq dataframe 
  ## write_result = TRUE or FALSE (default to FALSE), determines if the result is write to environment. 
  ## tls_signature = c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS", "CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3")
  ## tmb_col = c("Diagnosis Age", "Sex", "Fraction Genome Altered", "Mutation Count", "TMB (nonsynonymous)") 
  ## gene_pat = "^HLA.*[ABC]$|^CD274$|^GZMB$|^GNLY$|^PRF1$"
  ## mcp_pat = c("t.cell", "cd8_t.cell", "cytotoxic_lymphocytes", "b_lineage", "nk_cells", "monocytic_lineage", "myeloid_dendritic_cells", "neutrophils", "endothelial_cells", "fibroblasts")
  
  ## Read in samples and data cleaning -----------------------
  sample_names <- get(cancer, envir = .GlobalEnv)
  rna_seq <- get(rna_seq, envir = .GlobalEnv)
  sample_names <- intersect(sample_names, colnames(rna_seq)) ## this will be our new sample_names 
  rna_seq <- rna_seq[, sample_names] ## we will only preserve cancer patients we care about, this is a DATAFRAME!!
  
  ## MCP counter ---------------------------- 
  mcp_results <- MCPcounter.estimate(
    expression = rna_seq, ## dataframe rna_seq 
    featuresType = "HUGO_symbols")
  rownames(mcp_results) <- c("t.cell", "cd8_t.cell", "cytotoxic_lymphocytes", "b_lineage", "nk_cells", 
                             "monocytic_lineage", "myeloid_dendritic_cells", "neutrophils", "endothelial_cells", 
                             "fibroblasts")
  mcp_results <- as.data.frame(t(mcp_results))
  mcp_results <- select(mcp_results, all_of(mcp_pat))
  
  ## TLS status -----------------------
  tls_signature <- tls_sig  ## from literature 
  tls_rna_seq <- rna_seq[tls_signature, ] ## isolate their rna_expression
  tls_rna_seq <- as.matrix(t(tls_rna_seq)) ## transform it into a matrix 
  tls_rna_seq <- scale(tls_rna_seq, center = F) ## scale it 
  tls_rna_seq <- as.data.frame(tls_rna_seq) ## transform it back into a dataframe 
  tls_rna_seq$TLS_mean <- apply(tls_rna_seq, 1, mean) ## calculate the TLS by mean 
  
  ## TMB calculation ----------------------
  tmb <- read_tsv(paste0(MUTATION, "/", cancer,"_tmb.tsv"), show_col_types = FALSE) ## the read in is a tibble 
  tmb <- as.data.frame(tmb)
  rownames(tmb) <- tmb$`Sample ID`
  sample_names <- intersect(sample_names, rownames(tmb)) ## there are fewer patients than the original patient list. 
  tmb <- tmb[sample_names, tmb_col]
  
  ## PDL1 exp, GZMB, GNLY, PRF1, MHC1 expression ----------------------
  genes <- which(grepl(gene_pat, rownames(RNA_seq)))
  spec_gene <- rna_seq[genes, ] 
  MHCI <- colMeans(spec_gene[rownames(spec_gene) %in% c("HLA-A", "HLA-B", "HLA-C"), ], na.rm = TRUE)
  MHCII <- colMeans(spec_gene[which(grepl("^HLA-D", rownames(spec_gene))), ], na.rm = TRUE)
  spec_gene <- rbind(spec_gene, MHCI = MHCI, MHCII = MHCII)
  spec_gene <- as.data.frame(t(spec_gene), )
  spec_gene <- spec_gene[, which(!grepl("^HLA", colnames(spec_gene)))]
  





