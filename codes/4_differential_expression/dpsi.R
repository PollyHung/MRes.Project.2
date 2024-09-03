library(dplyr)
library(tidyr)
library(tidyverse)

setwd("/Volumes/Wild_Flower/OV_SpliceVariants/")

## read in result 
dpsi_gc <- read.table("data/5_differential_expression/DAS/long_reads/diffSplice_gc.dpsi", sep="\t") 
colnames(dpsi_gc) <- c("dPSI","p_value")
dpsi_gc_sig <- dpsi_gc %>% dplyr::filter(p_value < 0.05)
dpsi_gc_sig$rownames <- rownames(dpsi_gc_sig)

## clean up the dpsi event table 
dpsi_gc_sig <- dpsi_gc_sig %>% separate(col = rownames, sep = ";", 
                                        into = c("pacbio_id", "event"), 
                                        remove = FALSE)
dpsi_gc_sig <- dpsi_gc_sig %>% separate(col = event, 
                                        sep = ":", 
                                        into = c("AS_type", "chromosome", "coordinate"), 
                                        extra = "merge")
dpsi_gc_sig$strand <- str_extract(dpsi_gc_sig$coordinate, "[+-]$")
dpsi_gc_sig$coordinate <- gsub(":[+-]$", "", dpsi_gc_sig$coordinate)

## map the pacbio_id to gene names 
classifications <- read.table("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt", header = T) %>% 
  dplyr::filter(filter_result == "Isoform")
pacbio_to_gene <- data.frame(isoform_id = sub("\\.[^.]*$", "", classifications$isoform), 
                             gene_id = classifications$associated_gene) %>% distinct()
pacbio_to_gene <- pacbio_to_gene %>% dplyr::filter(isoform_id %in% dpsi_gc_sig$pacbio_id)
dpsi_gc_sig <- merge(dpsi_gc_sig, pacbio_to_gene, by.x = "pacbio_id", by.y = "isoform_id")
dpsi_gc_sig$event_id <- paste0(dpsi_gc_sig$gene_id, ":",dpsi_gc_sig$AS_type)

## add average psi values 
PSI <- read.table("data/5_differential_expression/DAS/long_reads/diffSplice.psivec")
PSI <- PSI %>% 
  dplyr::mutate(PSI_tumour_avg = rowMeans(dplyr::select(., contains("tumour")), na.rm = TRUE)) %>% 
  dplyr::mutate(PSI_normal_avg = rowMeans(dplyr::select(., contains("normal")), na.rm = TRUE)) 
colnames(PSI)[1:25] <- c(colnames(read.table("data/5_differential_expression/DAS/long_reads/psi_tumour.psi", check.names = F)), 
                         colnames(read.table("data/5_differential_expression/DAS/long_reads/psi_normal.psi", check.names = F)))
dpsi_gc_sig <- merge(dpsi_gc_sig, PSI[, c("PSI_tumour_avg", "PSI_normal_avg")], by.x = "rownames", by.y = "row.names")

## plot
dpsi_gc_sig$`-log10.pval` <- -log10(dpsi_gc_sig$p_value)
#write.table(dpsi_gc_sig, "data/5_differential_expression/DAS/long_reads/result_dpsi_gc_sig.txt", 
#            sep = "\t", quote = F, row.names = F, col.names = T)
dpsi_gc_sig$color <- ifelse(dpsi_gc_sig$dPSI > 0, "+", "-")
dpsi_gc_sig$event_id <- paste(dpsi_gc_sig$event_id, dpsi_gc_sig$chromosome, 
                              dpsi_gc_sig$coordinate, dpsi_gc_sig$strand, sep = ":")
dpsi_gc_sig$event_id <- factor(dpsi_gc_sig$event_id, 
                               levels = dpsi_gc_sig$event_id[order(abs(dpsi_gc_sig$dPSI), decreasing = FALSE)])

# Create the plot
p <- ggplot(dpsi_gc_sig, aes(x = event_id, y = abs(dPSI), fill = color)) +
  geom_bar(stat = "identity", width = 0.7) +   
  coord_flip() +                               
  scale_fill_manual(values = c("-" = "skyblue", "+" = "salmon"), name = "Sign") +  
  labs(x = "AS Event", y = "Absolute dPSI", title = NULL) +  
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12))
ggsave("plot/DAS/01.png", p, width = 9, height = 3.5, units = "in", dpi = 600)









