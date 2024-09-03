source("code/0_pkgs&func.R")

## Load Files ==================================================================
# [1] salmon 
files <- list.files(path = "data/3_salmon/output", pattern = ".sf", full.names = TRUE, 
                    recursive = TRUE)
names(files) <- stringr::str_split(files, pattern = "/", simplify = TRUE)[,4] %>%
  stringr::str_replace("_quant", "")
txi.salmon <- tximport(files, type = "salmon", txIn=TRUE, txOut = TRUE)

# [2] sqanti 
classification <- read.delim("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt") %>% 
  dplyr::filter(filter_result == "Isoform")
isoform_id <- intersect(rownames(txi.salmon$counts), classification$isoform)
count_matrix <- txi.salmon$abundance[isoform_id, ]

# [3] patient sheet 
survival <- read.delim("data/0_samples/hh_ova_surv.txt")
clinical <- read.delim("data/0_samples/hh_ova.txt")
immune <- read.delim("data/0_samples/hh_ova_immune.txt")
mapping <- read.csv("data/0_samples/mapping.csv") %>% 
  dplyr::rename("Row.names" = "ICBRC", "SR.id" = "X.1") %>% 
  dplyr::select(Row.names, SR.id) %>% 
  dplyr::filter(SR.id != "")
metadata <- merge(immune, clinical, by="row.names") 
metadata <- left_join(mapping, metadata, by="Row.names")
rownames(metadata) <- metadata$SR.id

# [4] alignment 
SR_id <- intersect(metadata$SR.id, colnames(count_matrix))
count_matrix <- count_matrix[, SR_id]
metadata <- metadata[SR_id, ]

## Survival Analysis ===========================================================
# [1] overall survival
isoform_id <- read.csv("data/6_survival/os.csv") %>% 
  dplyr::filter(p_adj < 0.05) %>% 
  dplyr::select(X) %>% unlist() %>% unname()
final_counts <- count_matrix[isoform_id, ]

os <- matrix(NA, nrow = nrow(final_counts), ncol = 4) %>% as.data.frame()
colnames(os) <- c("HR", "upper_95", "lower_95", "p_val")
rownames(os) <- rownames(final_counts)
surv_obj <- Surv(metadata$os_time, metadata$os_event)

for(isoform_id in rownames(final_counts)){
  tryCatch({
    # survival model 
    coxph_model <- coxph(surv_obj ~ final_counts[isoform_id, ]+age+stage, data = metadata)
    model <- summary(coxph_model) 
    # record result 
    os[isoform_id, "HR"] <- model$coefficients[1, "coef"]
    os[isoform_id, "p_val"] <- model$coefficients[1, "Pr(>|z|)"]
    os[isoform_id, "upper_95"] <- model$conf.int[1, "upper .95"]
    os[isoform_id, "lower_95"] <- model$conf.int[1, "lower .95"]
  }, error = function(e) {
    # Handle the error (e.g., print a message)
    cat("Error with isoform", isoform_id, ":", conditionMessage(e), "\n")
  })
}
os$p_adj <- p.adjust(os$p_val, method = "fdr")
os_sig <- os %>% dplyr::filter(p_adj < 0.25)


# [2] progression free survival 
isoform_id <- read.csv("data/6_survival/pfs_sig.csv") %>% dplyr::select(X) %>% unlist() %>% unname()
final_counts <- count_matrix[isoform_id, ]

pfs <- matrix(NA, nrow = nrow(final_counts), ncol = 4) %>% as.data.frame()
colnames(pfs) <- c("HR", "upper_95", "lower_95", "p_val")
rownames(pfs) <- rownames(final_counts)
surv_obj <- Surv(metadata$pfs_time, metadata$pfs_event)

# loop over 
for(isoform_id in rownames(final_counts)){
  tryCatch({
    # survival model 
    coxph_model <- coxph(surv_obj ~ final_counts[isoform_id, ]+age+stage, data = metadata)
    model <- summary(coxph_model) 
    # record result 
    pfs[isoform_id, "HR"] <- model$coefficients[1, "coef"]
    pfs[isoform_id, "p_val"] <- model$coefficients[1, "Pr(>|z|)"]
    pfs[isoform_id, "upper_95"] <- model$conf.int[1, "upper .95"]
    pfs[isoform_id, "lower_95"] <- model$conf.int[1, "lower .95"]
  }, error = function(e) {
    # Handle the error (e.g., print a message)
    cat("Error with isoform", isoform_id, ":", conditionMessage(e), "\n")
  })
}
pfs$p_adj <- p.adjust(pfs$p_val, method = "fdr")
pfs_sig <- pfs %>% dplyr::filter(p_adj < 0.25)


