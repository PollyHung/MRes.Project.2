## This code is to summarise the findings together to one dataframe 
## Focusing on functional annotations 
## https://services.healthtech.dtu.dk/
source("code/0_pkgs&func.R")


## SQANTI3 
classification <- read.delim("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt") %>% 
  dplyr::filter(filter_result == "Isoform")

## TransDecoder: what is the corresponding protein sequence? 
if(!file.exists("data/4_transcoder/transdecoder_map.csv")){
  # Load data
  pep <- read.fasta(file = file.path("data/4_transcoder/qc_corrected.fasta.transdecoder.pep"), seqtype = "AA", forceDNAtolower = FALSE)
  blastp_out <- read.table(file = file.path("data/4_transcoder/blastp.out"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  gencode <- import.gff("data/0_references/gencode.v39.chr_patch_hapl_scaff.annotation.sorted.gtf")
  uniprot2ENSG <- read.csv("data/0_references/mapping_uniprot_to_ensembl.csv", header = TRUE)
  
  # Extract peptide names and sequences
  pep_names <- getName(pep)
  pep_sequence <- seqinr::getSequence(pep)
  pb_ids <- gsub(".p[0-9]", "", pep_names, perl = TRUE)
  
  # Filter sequences passing QC
  pass_QC <- pb_ids %in% classification$isoform
  pep_QCpass_transcripts <- pep[pass_QC]
  annotation <- getAnnot(pep_QCpass_transcripts)
  pb_ids_pep_QCpass <- gsub(".p[0-9]", "", getName(pep_QCpass_transcripts), perl = TRUE)
  
  # Modify peptide sequences
  pep_modified <- lapply(pep_sequence[pass_QC], function(seq) {
    seq_a <- paste(toupper(as.character(unlist(seq))), collapse = "")
    gsub("\\*", "", seq_a)
  })
  
  # Write modified sequences to file
  write.fasta(sequences = pep_modified, names = getName(pep_QCpass_transcripts), file.out = file.path("data/4_transcoder/QC_pass_transcripts.transdecoder.fasta.pep"), as.string = TRUE)
  
  # Determine ORF types
  has_blastp <- getName(pep_QCpass_transcripts) %in% blastp_out$V1
  ORF_type <- rep("NA", length(pep_QCpass_transcripts))
  annotation_types <- c("type:complete", "type:5prime_partial", "type:3prime_partial", "type:internal")
  type_labels <- c("complete", "5prime_partial", "3prime_partial", "internal")
  for (i in seq_along(annotation_types)) {
    ORF_type[grep(annotation_types[i], annotation, ignore.case = TRUE)] <- type_labels[i]
  }
  
  # Extract annotation details
  uniprot_id <- uniprot_id <- lapply(annotation, function(annot){ # e.g. Q9H3E8
    res <- strsplit(annot[[1]], ",", perl = TRUE) 
    res2 <- strsplit(res[[1]][3], "\\|")
    uid <- res2[[1]][1]
  })
  
  pident <- lapply(annotation, function(annot){
    res <- strsplit(annot[[1]], ",", perl = TRUE)
    res2 <- strsplit(res[[1]][3], "\\|")
    res2[[1]][2]
  })
  
  evalue <- as.numeric(sapply(annotation, function(annot) {
    res <- strsplit(annot[[1]], ",", perl = TRUE)
    res2 <- strsplit(res[[1]][3], "\\|")
    res3 <- strsplit(res2[[1]][3], " ")
    res3[[1]][1]
  }))
  
  # Extract other details
  tscore <- lapply(annotation, function(annot){
    res <- strsplit(annot[[1]], "score=")
    res2 <- strsplit(res[[1]][2], ",")
    as.numeric(res2[[1]][1])
  })
  
  tlength <- lapply(annotation, function(annot){
    res <- strsplit(annot[[1]], "len:")
    res2 <- strsplit(res[[1]][2], " ")
    as.numeric(res2[[1]][1])
  })
  
  # Determine strand orientation
  tstrand <- sapply(annotation, function(annot) {
    res <- strsplit(annot[[1]], ",")
    if (grepl("\\(\\+\\)", res[[1]][1], perl = TRUE)) "+" else "-"
  })
  
  # Map SQANTI associated genes
  idx_sq <- match(pb_ids_pep_QCpass, classification$isoform)
  Sqanti_associated_gene <- classification$associated_gene[idx_sq]
  Sqanti_associated_gene_simp <- sapply(Sqanti_associated_gene, function(id) strsplit(as.character(id), "\\.")[[1]][1])
  
  # Annotate with Gencode
  gencode_genes <- gencode[gencode$type == "gene", ]
  idx_gencode <- match(Sqanti_associated_gene, gencode_genes$gene_name)
  Gencode_gene_type <- gencode_genes$gene_type[idx_gencode]
  Gencode_gene_name <- gencode_genes$gene_name[idx_gencode]
  Gencode_gene_id <- gencode_genes$gene_id[idx_gencode]
  novel_genes_count <- sum(is.na(Gencode_gene_id))
  
  # Annotate with Uniprot to Ensembl Gene ID
  uniprot_id_modified <- sapply(uniprot_id, function(id) strsplit(as.character(id), "-")[[1]][1])
  id_mapping <- match(uniprot_id_modified, uniprot2ENSG$uniprot_gn_id)
  Ensembl_geneID <- uniprot2ENSG$ensembl_gene_id[id_mapping]
  
  # Merging results into one table
  transdecoder_map <- data.frame(
    PBid = pb_ids_pep_QCpass,
    PeptideId = getName(pep_QCpass_transcripts),
    Blastp_match = unlist(uniprot_id),
    Blastp_evalue = evalue,
    Blastp_Uniprot2Ensembl = Ensembl_geneID,
    Transdecoder_ORF_type = ORF_type,
    Transdecoder_score = unlist(tscore),
    Transdecoder_ORF_length = unlist(tlength),
    Transdecoder_ORF_strand = tstrand,
    Gencode_gene_name = Gencode_gene_name,
    Gencode_gene_type = Gencode_gene_type,
    stringsAsFactors = FALSE
  )
  
  write.csv(transdecoder_map, "data/4_transcoder/transdecoder_map.csv")
  rm(list = ls())
} else {
  transdecoder <- read.csv("data/4_transcoder/transdecoder_map.csv", row.names = 1)
}

## TMHMM: is this protein surface protein? 
tmhmm <- read.table("data/4_transcoder/out_tmhmm.txt", header = F) %>% 
  dplyr::mutate(isoform = gsub("\\.[^.]*$", "", V1)) %>%  ## removes the .p1 peptide number 
  dplyr::select(isoform, V1, V3) %>% 
  dplyr::rename(peptide_id = V3, tmhmm_region = V3) 

## PFAM domains: what protein domain is this open reading frame likely be? 
pfam <- read_domtblout("data/4_transcoder/pfam.domtblout") %>% 
  as.data.frame() %>% 
  dplyr::filter(domain_ievalue < 10E-10) %>% 
  dplyr::mutate(isoform = gsub("\\.[^.]*$", "", query_name)) %>% 
  dplyr::select(isoform, query_name, domain_name, domain_score, description) %>% 
  dplyr::rename(peptide_id = query_name, pfam_name = domain_name, 
                pfam_score = domain_score, pfam_descrption = description)

## Deeploc2: where is this protein located subcellular? 



## Combine them together 

































