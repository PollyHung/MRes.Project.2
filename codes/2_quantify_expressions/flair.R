## library packages 
library(seqinr)
library(dplyr)
library(magrittr)
library(rtracklayer)
library(Biostrings)
library(doParallel)
library(parallel)
library(rhmmer)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(scales)
library(GenomicRanges)
library(tidyr)
library(purrr)
library(data.table)
library(NOISeq)

# prepare manifest file for flair 
# --reads_manifest    Tab delimited file containing sample id, condition, batch,
#                     reads.fq, where reads.fq is the path to the sample fastq file.
# sample1      condition1      batch1  mydata/sample1.fq
# sample2      condition1      batch1  mydata/sample2.fq
# sample3      condition1      batch1  mydata/sample3.fq
# sample4      condition2      batch1  mydata/sample4.fq
# sample5      condition2      batch1  mydata/sample5.fq
# sample6      condition2      batch1  mydata/sample6.fq

## create the manifest table 
if(!file.exists("~/MRes.project.2/codes/0_sample_id/complete_short_reads.txt")){
  mapping <- read.table("~/MRes.project.2/codes/0_sample_id/sample_list.txt") %>% 
    dplyr::mutate(condition = case_when(grepl("FN", V1) ~ "normal",
                                        TRUE ~ "tumour")) %>% 
    dplyr::mutate(batch = "batch1") %>% 
    dplyr::rename(`sample id`=V1)
  mapping$reads.fq <- paste0("/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/", 
                             mapping$`sample id`, 
                             "/isoseq3/ccs.fq")
  write.table(x = mapping, 
              file = "~/MRes.project.2/codes/0_sample_id/flair_manifest.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE)
  
  short_reads <- read.table("/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/short_reads/RNA-seq/G72_run990/G72_run990_TPU_md5sum.txt")
  short_reads$comp_dir <- paste0("/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/short_reads/RNA-seq/G72_run990/", 
                                 short_reads$V2)
  paired_paths <- data.frame(R1 = character(), R2 = character(), stringsAsFactors = FALSE)
  file_paths <- short_reads$comp_dir
  for (i in seq(1, length(file_paths), by = 2)) {
    r1_file <- file_paths[i]
    r2_file <- file_paths[i + 1]
    paired_paths <- rbind(paired_paths, data.frame(R1 = r1_file, R2 = r2_file, stringsAsFactors = FALSE))
  }
  
  write.table(paired_paths, "~/MRes.project.2/codes/0_sample_id/complete_short_reads.txt", 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
}

## perform flair quantify in bash script 
system("./MRes.project.2/codes/3_differential_expression/flair.sh")

## tappAS ======================================================================
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"

# expression matrix
set.seed(50)
flair_count <- read.table(file.path(tappAS, "flair/flair_counts_filtered.tsv"), 
                          header = TRUE, check.names = FALSE, row.names = 1)
sample_list <- read.table("~/MRes.project.2/codes/0_sample_id/sample_list.txt") %>% unlist %>% unname
colnames(flair_count) <- sample_list
flair_count <- as.matrix(flair_count)

# annotation gff3
annotation_gff <- readGFF(file.path(tappAS, "sqanti3/qc.gff3"))

# experiment design 
experiment_resign <- read.delim("~/MRes.project.2/codes/0_sample_id/flair_manifest.txt", header = FALSE) %>% 
  dplyr::select(V1, V2) %>% 
  dplyr::rename(sample = V1, group = V2)

# filtering was done by the SQANTI3 step, now we will normalize by TMM normalization from NOISeq package 
normalized_exp <- NOISeq::tmm(datos = flair_count, 
                              long = nrow(flair_count), 
                              logratioTrim = 0.3, 
                              sumTrim = 0.05, 
                              k = 0, 
                              lc = 0, 
                              refColumn = 1)














