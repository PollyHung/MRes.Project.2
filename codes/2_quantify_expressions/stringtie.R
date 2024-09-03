library(dplyr)
library(tidyverse)

gtf_dir <- "/Volumes/Wild_Flower/OV_SpliceVariants/data/3_stringtie/stringtie_result"
gtf_files <- list.files(gtf_dir, pattern = "*.gtf", full.names = TRUE)

extract_tpm <- function(gtf_file, filter = "PacBio") { # StringTie
  gtf_data <- read.table(gtf_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  sample_id = gsub("_stringtie", "", tools::file_path_sans_ext(basename(gtf_file)))
  
  gtf_data <- gtf_data %>% 
    dplyr::filter(V2 == filter, V3 == "transcript") %>% 
    separate(V9, sep = "; ", into = c("gene_id", "transcript_id", "coverage", "FPKM", "TPM")) %>% 
    dplyr::mutate(TPM = gsub("TPM |;", "", TPM)) %>% 
    dplyr::mutate(transcript_id = gsub("transcript_id ", "", transcript_id)) %>% 
    dplyr::select(transcript_id, TPM) %>% 
    dplyr::rename(!!sample_id := "TPM")
  
  return(gtf_data)
}

pacbio <- lapply(gtf_files, extract_tpm)
pacbio_mtx <- Reduce(function(x, y) full_join(x, y, by = "transcript_id"), pacbio)
pacbio_mtx[2:ncol(pacbio_mtx)] <- apply(pacbio_mtx[2:ncol(pacbio_mtx)], 2, function(x){as.numeric(x)})
rownames(pacbio_mtx) <- pacbio_mtx$transcript_id
pacbio_mtx$transcript_id <- NULL

