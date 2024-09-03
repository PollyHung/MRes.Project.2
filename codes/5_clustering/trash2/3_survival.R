library(tidyverse)
library(dplyr)
library(survival)

setwd("/Volumes/Wild_Flower/OV_SpliceVariants/data/")
mapping <- read.table("1_sampleSheet/mapping_complete.txt", header = T, sep = "\t")
survival <- read.csv("1_sampleSheet/hh_ova_surv.csv", row.names = 1) 
clinical <- read.table("1_sampleSheet/hh_ova.txt", row.names = NULL) %>% 
  dplyr::select(row.names, age, stage) %>% dplyr::rename(copy_number = row.names)
load("6_clustering/output/clust_lr.RData")
load("6_association/sample_data/psi_knn_imp.Rd")

AS_events <- clust$ts_clusters$Junction#[clust$ts_clusters$Counts > 10]


survival <- merge(mapping, survival, by.x = "copy_number", by.y = "sample_id")
rownames(survival) <- survival$long_read
sample_id <- intersect(survival$long_read, colnames(psi_lr.imp))
psi_lr.imp <- psi_lr.imp[AS_events, sample_id] %>% as.data.frame()
survival <- survival[sample_id, ]
clinical <- merge(mapping, clinical, by = "copy_number")
rownames(clinical) <- clinical$long_read
clinical <- clinical[sample_id, ]
 
pfs_obj <- Surv(time = survival$pfs_time, event = survival$pfs_event)
pfs_results <- matrix(NA, nrow = nrow(psi_lr.imp), ncol = 4) %>% as.data.frame()
colnames(pfs_results) <- c("HR", "upper.95", "lower.95", "p_val")
rownames(pfs_results) <- rownames(psi_lr.imp)

for(i in rownames(psi_lr.imp)){
  AS_exp <- unlist(psi_lr.imp[i, ])
  AS_exp_di <- ifelse(AS_exp > mean(AS_exp), "high", "low")
  test <- coxph(pfs_obj~AS_exp)
  test <- summary(test)
  test2 <- survdiff(pfs_obj~AS_exp_di)
  pfs_results[i, "HR"] <- test$coefficients[1, 1]
  pfs_results[i, "p_val"] <- test2$pvalue
  pfs_results[i, "lower.95"] <- test$conf.int[1, 3]
  pfs_results[i, "upper.95"] <- test$conf.int[1, 4]
}

pfs_results$p_adj <- p.adjust(pfs_results$p_val, method = "fdr")
pfs_results$pfs_prognosis <- ifelse(pfs_results$p_adj > 0.25, "ns", 
                                    ifelse(pfs_results$HR > 1, "unfavourable", "favourable"))


os_obj <- Surv(time = survival$os_time, event = survival$os_event)
os_results <- matrix(NA, nrow = nrow(psi_lr.imp), ncol = 4) %>% as.data.frame()
colnames(os_results) <- c("HR", "upper.95", "lower.95", "p_val")
rownames(os_results) <- rownames(psi_lr.imp)

for(i in rownames(psi_lr.imp)){
  AS_exp <- unlist(psi_lr.imp[i, ])
  AS_exp_di <- ifelse(AS_exp > mean(AS_exp), "high", "low")
  test <- coxph(os_obj~AS_exp)
  test <- summary(test)
  test2 <- survdiff(os_obj~AS_exp_di)
  os_results[i, "HR"] <- test$coefficients[1, 1]
  os_results[i, "p_val"] <- test2$pvalue
  os_results[i, "lower.95"] <- test$conf.int[1, 3]
  os_results[i, "upper.95"] <- test$conf.int[1, 4]
}

os_results$p_adj <- p.adjust(os_results$p_val, method = "fdr")
os_results$os_prognosis <- ifelse(os_results$p_adj > 0.25, "ns", 
                                  ifelse(os_results$HR > 1, "unfavourable", "favourable"))

survival <- data.frame(AS_event = rownames(os_results), 
                       os_prognosis = os_results$os_prognosis, 
                       pfs_prognosis = pfs_results$pfs_prognosis)
save(list = c("survival", "os_results", "pfs_results"), 
     file = "6_association/survival_association.Rd")
