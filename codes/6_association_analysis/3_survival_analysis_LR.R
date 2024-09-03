source("code/1_pkgs&func.R")

## DEXSeq2 =====================================================================
# read in -----------
counts <- read.delim("data/3_flair/flair_counts_filtered.tsv", row.names = 1, check.names = FALSE) 
mapping <- read.csv("data/0_samples/mapping_table.csv", row.names = 1) %>% 
  dplyr::mutate(long_read = gsub("_", "-", long_read))
classification <- read.delim("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt") %>% 
  dplyr::filter(filter_result == "Isoform")
survival <- read.delim("data/0_samples/hh_ova_surv.txt")
immune <- read.delim("data/0_samples/hh_ova_immune.txt")
clinical <- read.delim("data/0_samples/hh_ova.txt")

# mutate -----------
indices <- match(mapping$long_read, colnames(counts))
colnames(counts)[indices] <- mapping$copy_number
counts <- counts[mapping$copy_number]
metadata <- merge(survival[mapping$copy_number, ], immune[mapping$copy_number, ], by="row.names")
metadata[1, "Row.names"] <- setdiff(colnames(counts), metadata$Row.names)
rownames(metadata) <- metadata$Row.names
metadata$Row.names <- NULL
clinical <- clinical[mapping$copy_number, c("age", "stage")]
rownames(clinical)[19] <- "X2316"
metadata <- merge(clinical, metadata, by="row.names")
rownames(metadata) <- metadata$Row.names
metadata <- metadata[, 2:ncol(metadata)]
metadata <- metadata[colnames(counts), ]

# construct -------------
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~1)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# normalize -------------
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
      # dimension reduced from 26625 to 11892

# variance stabilizing transformation ------------
vsd <- vst(dds, blind = FALSE)
transformed_counts <- assay(vsd)

# matrix of transformed counts --------------
final_counts <- transformed_counts[, rownames(metadata)]

## survival analysis ===========================================================
# overall survival ------------
# empty table 
os <- matrix(NA, nrow = nrow(final_counts), ncol = 4) %>% as.data.frame()
colnames(os) <- c("HR", "upper_95", "lower_95", "p_val")
rownames(os) <- rownames(final_counts)

# survival object 
surv_obj <- Surv(metadata$os_time, metadata$os_event)

# loop over 
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

# multiple testing adjustments 
os$p_adj <- p.adjust(os$p_val, method = "fdr")
os_sig <- os %>% dplyr::filter(p_adj < 0.05)
write.csv(os, "data/6_survival/os.csv")
write.csv(os_sig, "data/6_survival/os_sig.csv")

# progression free survival ---------------
pfs <- matrix(NA, nrow = nrow(final_counts), ncol = 4) %>% as.data.frame()
colnames(pfs) <- c("HR", "upper_95", "lower_95", "p_val")
rownames(pfs) <- rownames(final_counts)

# survival object 
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

# multiple testing adjustments 
pfs$p_adj <- p.adjust(pfs$p_val, method = "fdr")
pfs_sig <- pfs %>% dplyr::filter(p_adj < 0.05)
write.csv(pfs, "data/6_survival/pfs.csv")
write.csv(pfs_sig, "data/6_survival/pfs_sig.csv")



## Immune response correlation =================================================
# without selection: ----------------
survival_select <- final_counts
immune_response <- matrix(NA, nrow = nrow(survival_select), ncol = 11) %>% as.data.frame()
colnames(immune_response) <- colnames(metadata[7:ncol(metadata)])
rownames(immune_response) <- rownames(survival_select)

# Function to perform correlation calculation for each isoform expression
lm_by_row <- function(immune_fac, return_item = "p_val") {
  # Create temp dataframe
  temp <- matrix(NA, nrow = nrow(survival_select), ncol = 2) %>% as.data.frame()
  rownames(temp) <- rownames(survival_select)
  colnames(temp) <- c("p_val", "p_adj")
  
  # Get immune column
  immune_col <- metadata[, immune_fac]
  names(immune_col) <- rownames(metadata)
  
  # Loop over each isoform
  for (isoform_id in rownames(survival_select)) {
    tryCatch({
      cor_model <- cor.test(survival_select[isoform_id, ], immune_col, method = "spearman", use = "complete.obs")
      temp[isoform_id, "p_val"] <- cor_model$p.value
    }, error = function(e) {
      cat("Error with isoform", isoform_id, ":", conditionMessage(e), "\n")
    })
  }
  
  # Adjust p-values using FDR
  if (return_item == "p_adj") {
    temp$p_adj <- p.adjust(temp$p_val, method = "fdr")
    return(temp$p_adj)
  } else if (return_item == "p_val") {
    return(temp$p_val)
  } else {
    stop("Invalid return_item specified. Use 'p_val' or 'p_adj'.")
  }
}

# Apply the function to each relevant column in the metadata
for (immune_fac in colnames(metadata)[7:ncol(metadata)]) {
  immune_response[, immune_fac] <- lm_by_row(immune_fac, return_item = "p_adj")
}

# Print the results
write.csv(immune_response, file = "data/6_survival/immune_feature_p_adj.csv")

# with selection: ----------------
survival_select <-  final_counts[unique(c(rownames(os_sig), rownames(pfs_sig))), ] 
immune_response <- matrix(NA, nrow = nrow(survival_select), ncol = 8) %>% as.data.frame()
colnames(immune_response) <- colnames(metadata[7:ncol(metadata)])
rownames(immune_response) <- rownames(survival_select)

# Apply the function to each relevant column in the metadata
for (immune_fac in colnames(metadata)[7:ncol(metadata)]) {
  immune_response[, immune_fac] <- lm_by_row(immune_fac, return_item = "p_adj")
}

temp <- immune_response %>% rowwise() %>% dplyr::filter(any(c_across(where(is.numeric)) < 0.10))

# Print the results
write.csv(immune_response, file = "data/6_survival/selected_isoform_immune_feature_p_adj.csv")

## Plots =======================================================================
load("data/6_survival/survival.RData")

# survival
metadata$sample <- rownames(metadata)
survival_data <- metadata %>%
  select(sample, os_event, os_time, pfs_event, pfs_time) %>%
  mutate(os_time = as.numeric(os_time))
survival_data <- survival_data[match(colnames(expression_data), survival_data$sample), ]

# function
plot_kaplan_meier <- function(isoform) {
  median_expr <- mean(expression_data[isoform, ])
  group <- ifelse(expression_data[isoform, ] >= median_expr, "High", "Low")
  fit <- survfit(surv_object ~ group)
  
  p <- ggsurvplot(fit, data = survival_data, 
                  pval = TRUE, 
                  risk.table = TRUE, 
                  title = paste("Kaplan-Meier Plot for", isoform),
                  xlab = "Time (days)", 
                  ylab = "Survival Probability",
                  legend.title = "Expression",
                  legend.labs = c("Low", "High"))
  return(p)
}

# plot 
surv_sig <- rownames(pfs_sig)
expression_data <- final_counts[surv_sig, ]
surv_object <- Surv(time = survival_data$os_time, event = survival_data$os_event)
km_plots_os <- list()
for (isoform in surv_sig) {
  km_plots_os[[isoform]] <- plot_kaplan_meier(isoform)
}

surv_sig <- rownames(pfs_sig)
expression_data <- final_counts[surv_sig, ]
surv_object <- Surv(time = survival_data$pfs_time, event = survival_data$pfs_event)
km_plots_pfs <- list()
for (isoform in surv_sig) {
  km_plots_pfs[[isoform]] <- plot_kaplan_meier(isoform)
}



for (isoform in names(km_plots_os)) {
  print(km_plots_os[[isoform]])
}





save.image("data/6_survival/survival.RData")

