source("code/0_pkgs&func.R")

## Analyse SUPPA2 ==============================================================
# [1] load files ------------
# differential expression p-values
dpsi <- read.table("data/5_differential_expression/SUPPA/diffSplice.dpsi", sep = "\t") %>% 
  dplyr::rename(dPSI = cancer.normal_dPSI, p_value = cancer.normal_p.val)

# psi of the events 
psi_events <- read.table("data/5_differential_expression/SUPPA/psi_local.psi", 
                         sep = "\t", check.names = FALSE) 

# tpm of the events
event_TPM <- read.table("data/5_differential_expression/SUPPA/diffSplice_avglogtpm.tab", 
                        sep = "\t", header = FALSE) %>% 
  dplyr::rename(event = V1, mean_TPM = V2)

# [2] merging dataframes ------------
merge1 <- merge(dpsi,psi_events,by="row.names")
final_table <- merge(merge1,event_TPM,by.x="Row.names",by.y="event")
rownames(final_table) <- final_table$Row.names
final_table <- final_table[,-1]
final_table <- final_table[!is.nan(final_table$dPSI), ]
final_table$cpval <- p.adjust(final_table$p_value, method = "bonferroni")
final_table$log10pval <- -log10(final_table$p_value)
final_table$logRNAc <-final_table$mean_TPM
final_table$rownames <- rownames(final_table)
final_table <- final_table %>%
  separate(rownames, into = c("PB_ID", "rest"), sep = ";") %>%
  separate(rest, into = c("splicing_event", "chromosome_location"), sep = ":", extra = "merge")

pie(table(final_table$splicing_event))

## use rMATs ===================================================================



