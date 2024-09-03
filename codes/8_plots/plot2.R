library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(magrittr)
library(survminer)
library(survival)
library(ggpattern)

colours = c("#FF8C9E", "#EF9C66", "#FFDA76", "#9CA777", "#41B3A2", "#378CE7", "#40679E", "#667BC6", "#D7C3F1", "#944E63")
color.heavy = c("#E07F80", "#95AAD3", "#F0CF7F", "#C195C4", "#94C47D", "#8785BA", "#E69965", "#72BEB7")
color.light = c("#EDAEAE", "#B9DBF4", "#F4E6AA", "#F0D8E9", "#C8DBC8", "#BAB4D8", "#F6BF93", "#B8E2DE")

source("~/Desktop/project/codes/code/8_plots/function.R")
setwd("~/Desktop/project/plots/8_plots/")
load("mapping.Rd")
load("psi.Rd")
load("metatable.Rd")
load("cor_mtx.Rd")
load("mutation.Rd")
load("copy_number.Rd")
load("immune_mtx.Rd")
load("gsea.Rd")

## Plot 1: Patients can be clustered by PSI pattern ----------------------------
# png("plot/plot.2/p1.png", width = 6.8, height = 6.8, units = "in", res = 600)
# heatmap(cor_mtx, margins = c(10, 10))
# dev.off()

heatmap_result <- heatmap(cor_mtx, margins = c(10, 10))
samples <- colnames(cor_mtx)[heatmap_result$colInd]

## Plot 2: Cluster 2 is associated with worse overall survival -----------------
## Overall Survival 
fit <- survfit(Surv(os_time, os_event) ~ as.factor(group), data = metatable)
p2 <- ggsurvplot_customised(dataframe = metatable, 
                            pval_coord = c(200, 0.25), 
                            break_time_by = 500, 
                            conf_int = TRUE, 
                            title = "Overall Survival", 
                            legend_title = "", 
                            subtitle = "Sample Size = 18", 
                            legend_labels = c("Group 1", "Group 2")) 
# png("plot/plot.2/p2.png", width = 2.3, height = 5, units = "in", res = 600)
# print(p2)
# dev.off()

## Progression Free Survival
fit <- survfit(Surv(pfs_time, pfs_event) ~ as.factor(group), data = metatable)
p3 <- ggsurvplot_customised(dataframe = metatable, 
                            pval_coord = c(200, 0.25), 
                            break_time_by = 500, 
                            conf_int = TRUE, 
                            title = "Progression Free Survival", 
                            legend_title = "", 
                            subtitle = "Sample Size = 17", 
                            legend_labels = c("Group 1", "Group 2")) 
# png("plot/plot.2/p3.png", width = 2.3, height = 5, units = "in", res = 600)
# print(p3)
# dev.off()


## Plot 3: Cluster 2 is associated with higher patient stage ------------
p4 <- ggplot(metatable, aes(x = group, y = stage)) +
  geom_boxplot(fill = "#378CE7", color = "black", linewidth = 0.25, outlier.size = 0.25) +
  labs(x = "Group", y = "Patient Stage") +
  theme_bw(base_size = 7) + 
  theme(axis.title = element_text(size = 6), 
        legend.key.size = unit(0.07, "in"), 
        legend.title = element_text(size = 6))
# ggsave("plot/plot.2/p7.png", p7, units = 'in',height = 1, width = 1, dpi = 600)

## Plot 4: Cluster 2 is not associated with any pre-identified common mutations in HGSOC -------------------------
p5.df <- mutations %>% 
  dplyr::filter(gene %in% c("BRCA1", "BRCA2", "CCNE1", "TP53", "CSMD3", 
                            "NF1", "CDK12", "FAT3", "GABRA6", "RB1", "PTEN")) %>% 
  dplyr::select(long_read, gene) %>% 
  distinct() %>% dplyr::rename(sample = long_read) %>% 
  dplyr::mutate(sample = gsub("_", "-", sample))
p5.df <- left_join(p5.df, condition_mtx[, c("sample", "group")], by = "sample")
p5.df <- p5.df %>% mutate(mutation = 1) %>% 
  pivot_wider(names_from = gene, values_from = mutation, values_fill = list(mutation = 0)) %>% 
  pivot_longer(cols = BRCA1:TP53, names_to = "gene", values_to = "mutation") 
p5.df$group <- ifelse(p5.df$group == "group_1", "G1", "G2")
p5.df$mutation <- ifelse(p5.df$mutation == "0", "wt", "mut")
p5 <- ggplot(p5.df, aes(x = group, fill = factor(mutation))) + 
  geom_bar(position = "fill", stat = "count") + 
  facet_wrap(~gene, nrow = 1) + 
  scale_fill_manual(values = c("wt" = "#40679E", "mut" = "#FF8C9E"), name = "Mutation") +
  labs(x = "Group", y = "Proportion") + 
  theme_bw(base_size = 7) + 
  theme(axis.title = element_text(size = 6), 
        legend.key.size = unit(0.07, "in"), 
        legend.title = element_text(size = 6))
# ggsave("plot/plot.2/p4.png", p4, units = "in", width = 3, height = 1.2, dpi = 600)

## Tile Plot 
p5.df$sample <- factor(p5.df$sample, levels = unique(p5.df$sample))
missing_samples <- setdiff(samples, p5.df$sample)
missing_df <- data.frame(gene = rep(unique(p5.df$gene), each = length(missing_samples)),
                         sample = rep(missing_samples, times = length(unique(p5.df$gene))),
                         mutation = "no data")
p6.df <- rbind(p5.df[, c("sample", "gene", "mutation")], missing_df)
p6.df$sample <- factor(p6.df$sample, levels = samples)
p6.df <- p6.df %>% dplyr::arrange(sample)
p6 <- ggplot(p6.df, aes(x = gene, y = sample, fill = factor(mutation))) + 
  geom_tile(color = "white") + 
  coord_flip() +
  scale_fill_manual(values = c("wt" = "#40679E", "mut" = "#FF8C9E", "no data" = "white"), name = "Mutation") +
  labs(y = "Sample") +
  theme_bw(base_size = 6) + 
  theme(axis.title = element_text(size = 6), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.key.size = unit(0.07, "in"), 
        legend.title = element_text(size = 6))
# ggsave("plot/plot.2/p5.png", p5, units = "in", width = 3, height = 1.2, dpi = 600)

# temp <- filter(p4.df, gene == "CCNE1")
# fisher.test(x = as.factor(temp$group), y = as.factor(temp$mutation))

## 


## Plot 

## Plot 5: Cluster 2 is associated with increase AS of metabolic related pathways ---------------------
p7.df <- gsea %>% dplyr::filter(`FDR q-val` <= 0.25)
p7.df$direction <- ifelse(p7.df$NES > 0, "Upregulated", "Downregulated")
p7.df$database <- unlist(lapply(p7.df$NAME, function(x){ str_split(x, pattern = "_")[[1]][1] }))
p7.df$alpha <- ifelse(p7.df$`FDR q-val` < 0.05, 1, 0.5)
p7.df$label <- c("ECM Receptor Interaction",
                 "ERK1 and ERK2 Cascade",
                 "Positive Regulation of MAPK Cascade",
                 "Regulation of Immune Effector Process",
                 "Positive Regulation of Immune Effector Process",
                 "Response to Molecule of Bacterial Origin",
                 "Regulation of Body Fluid Levels",
                 "MAPK Cascade",
                 "Immune Effector Process",
                 "Negative Regulation of Transport",
                 "Cellular Response to Organic Cyclic Compound",
                 "Cellular Response to Steroid Hormone Stimulus",
                 "Response to Steroid Hormone",
                 "Response to Ketone",
                 "Regulation of Microtubule-Based Process",
                 "Regulation of Microtubule Cytoskeleton Organization",
                 "Nucleoside Triphosphate Metabolic Process",
                 "Regulation of Canonical WNT Signaling Pathway",
                 "ATP Metabolic Process",
                 "Canonical WNT Signaling Pathway",
                 "Cytokinesis",
                 "Cellular Response to Hormone Stimulus",
                 "Response to Elevated Platelet Cytosolic CA2",
                 "Toll-Like Receptor Cascades",
                 "Platelet Activation, Signaling, and Aggregation",
                 "Nervous System Development",
                 "Signaling by Interleukins",
                 "Regulation of Expression of SLITs and ROBOs",
                 "Hemostasis",
                 "Nonsense-Mediated Decay (NMD)",
                 "Extracellular Matrix Organization",
                 "Vesicle-Mediated Transport",
                 "Selenoamino Acid Metabolism",
                 "Metabolism of Lipids")
p7 <- ggplot(p7.df, aes(x = reorder(label, -NES), y = NES, fill = direction, pattern = database, alpha = alpha)) +
  geom_bar_pattern(stat = "identity", 
                   position = "dodge", 
                   pattern_angle = 45,
                   pattern_color = "black", 
                   pattern_fill = "black", 
                   color = "black", 
                   pattern_spacing = 0.05) +
  scale_fill_manual(values = c("Upregulated" = "#FF8C9E", "Downregulated" = "#40679E")) +
  scale_pattern_manual(values = c("REACTOME" = "stripe", "KEGG" = "none", "GOBP" = "circle")) +
  scale_alpha_identity() +
  labs(x = "Pathway", y = "NES", fill = "Direction", pattern = "Database") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(angle = 90, vjust = 0.5))
ggsave("~/Desktop/project/figures/plot2/gsea.png", p7, width = 6.5, height = 4.5, units = "in", dpi = 600)


## Plot 6: Cluster 2 is associated with 






## Plot 4: Cluster 2 is associated with increase copy number in chromosome 3q -------------------------
## Copy Number 
# copy_no.arm$group <- ifelse(copy_no.arm$group == "group_1", "G1", "G2")
p6 <- ggplot(copy_no.arm, aes(x = group, y = `3q`)) +
  geom_boxplot(fill = "#667BC6", color = "black", linewidth = 0.25, outlier.size = 0.25) +
  labs(x = "Group", y = "3q copy no.") +
  theme_bw(base_size = 7) + 
  theme(axis.title = element_text(size = 6), 
        legend.key.size = unit(0.07, "in"), 
        legend.title = element_text(size = 6))
# ggsave("plot/plot.2/p6.png", p6, units = 'in',height = 1, width = 1, dpi = 600)


## Plot 6 ----------------------------------------------------------------------
## Association with Stage 
# metatable$group <- ifelse(metatable$group == "group_1", "G1", "G2")



## Plot 8 ----------------------------------------------------------------------
## Copy Number of these genes 
p8.df <- copy_no.gene[c("BRCA1", "BRCA2", "CCNE1", "TP53", "CSMD3", 
                        "NF1", "CDK12", "FAT3", "GABRA6", "RB1", "PTEN"), 2:14]
p8.df <- as.data.frame(t(p8.df)) %>% rownames_to_column(var = "sample")
p8.df <- left_join(p8.df, condition_mtx, by = "sample")
# p8.df <- as.data.frame(p8.df) %>% rownames_to_column(var = "gene")
# p8.df <- pivot_longer(p8.df, cols = -gene, names_to = "sample", values_to = "copy_number")
# p8.df$group <- ifelse(p8.df$group == "group_1", "G1", "G2")

# BRCA1  BRCA2  CCNE1   TP53  CSMD3    NF1  CDK12   FAT3 GABRA6    RB1   PTEN 
p8 <- ggplot(p8.df, aes(x = group, y = PTEN)) +
  geom_boxplot(fill = "#40679E", color = "black", linewidth = 0.25, outlier.size = 0.25) +
  labs(x = "Group", y = "PTEN copy no.") +
  theme_bw(base_size = 7) + 
  theme(axis.title = element_text(size = 6), 
        legend.key.size = unit(0.07, "in"), 
        legend.title = element_text(size = 6))
# ggsave("plot/plot.2/p5.png", p5, units = 'in',height = 1.2, width = 1.2, dpi = 600)

wilcox.test(p8.df$PTEN[p8.df$group == "group_1"], p8.df$PTEN[p8.df$group == "group_2"])
ggsave("plot/plot.2/p8.png", p8, units = 'in',height = 1, width = 1, dpi = 600)

## Plot 7 ----------------------------------------------------------------------
## Immune 
immune_features <- setdiff(colnames(immune_mtx), c("long_read", "Row.names", "group"))
group_1 <- filter(immune_mtx, group == "group_1") %>% pull(Row.names)
group_2 <- filter(immune_mtx, group == "group_2") %>% pull(Row.names)
for(i in immune_features){
  test <- wilcox.test(immune_mtx[group_1, i], immune_mtx[group_2, i])
  print(test$p.value)
}








dge <- data.frame(genes = rownames(log_tpm_mtx), 
                  log2FC = rep(1, nrow(log_tpm_mtx)), 
                  p_value = rep(1, nrow(log_tpm_mtx)))
rownames(dge) <- dge$genes
for(i in rownames(psi.data)){
  group_1 = psi.data[i, condition_mtx$sample[condition_mtx$group == "group_1"]] %>% as.numeric()
  group_2 = psi.data[i, condition_mtx$sample[condition_mtx$group == "group_2"]] %>% as.numeric()
  
  dge[i, "log2FC"] <- log2((mean(group_2, na.rm = TRUE) + 0.0001)/(mean(group_1, na.rm = TRUE)+0.0001))
  dge[i, "p_value"] <- wilcox.test(group_1, group_2)$p.value
}


dge$p.adj = p.adjust(dge$p_value, method = "fdr")
dge <- dplyr::arrange(dge, -log2FC)
write.table(dge[, c("genes", "log2FC")], "~/Desktop/dge_true.rnk", sep = "\t", quote = F, col.names = F, row.names = F)


setwd("~/Desktop/project/docs/data/5_differential_expression/GSEA_dge/")
gsea_dge <- rbind(read.delim("gsea_report_for_na_neg_1725134561410.tsv", check.names = F), 
                  read.delim("gsea_report_for_na_pos_1725134561410.tsv", check.names = F))
gsea_dge <- readxl::read_xlsx("~/Desktop/gsea.xlsx")

gsea_dge <- gsea_dge[, 1:11]
p7.df <- gsea_dge %>% dplyr::filter(`FDR q-val` <= 0.05)
p7.df <- gsea_dge %>%
  dplyr::filter(grepl("METABO", NAME)) %>%
  dplyr::filter(`FDR q-val` <= 0.25)
p7.df$direction <- ifelse(p7.df$NES > 0, "Upregulated", "Downregulated")
p7.df$database <- unlist(lapply(p7.df$NAME, function(x){ str_split(x, pattern = "_")[[1]][1] }))
p7.df$alpha <- ifelse(p7.df$`FDR q-val` < 0.05, 1, 0.5)
p7 <- ggplot(p7.df, aes(x = reorder(NAME, -NES), y = NES, fill = direction, pattern = database, alpha = alpha)) +
  geom_bar_pattern(stat = "identity", 
                   position = "dodge", 
                   pattern_angle = 45, 
                   pattern_color = "black", 
                   pattern_fill = "black", 
                   color = "black", 
                   pattern_spacing = 0.05) +
  scale_fill_manual(values = c("Upregulated" = "#FF8C9E", "Downregulated" = "#40679E")) +
  scale_pattern_manual(values = c("REACTOME" = "stripe", "KEGG" = "none", "GOBP" = "circle")) +
  #scale_alpha_identity() +
  labs(x = "Pathway", y = "NES", fill = "Direction", pattern = "Database") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(angle = 90, vjust = 0.5))
ggsave("~/Desktop/project/figures/plot2/gsea_dge.png", p7, width = 15, height = 15, units = "in", dpi = 600)

setwd("~/Desktop/project/docs/data/5_differential_expression/GSEA_SUPPA2/")
gsea_suppa <- rbind(read.delim("gsea_report_for_na_neg_1725135606691.tsv", check.names = F), 
                    read.delim("gsea_report_for_na_pos_1725135606691.tsv", check.names = F))
gsea_suppa  <- gsea_suppa[, 1:11]
p8.df <- gsea_suppa %>% dplyr::filter(`FDR q-val` <= 0.25)
p8.df$label <- c("Regulation of Immune Effector Process",              
                 "ERK1 and ERK2 Cascade",                              
                 "ECM Receptor Interaction",                           
                 "Positive Regulation of MAPK Cascade",                
                 "Response to Elevated Platelet Cytosolic Ca2+",    
                 "Response to Molecule of Bacterial Origin",           
                 "Negative Regulation of Transport",                   
                 "Positive Regulation of Immune Effector Process",     
                 "Toll-Like Receptor Cascades",                    
                 "Immune Effector Process",                            
                 "Regulation of Body Fluid Levels",                    
                 "MAPK Cascade",                                       
                 "Cellular Response to Organic Cyclic Compound",       
                 "Response to Steroid Hormone",                        
                 "Cellular Response to Steroid Hormone Stimulus",      
                 "Response to Ketone",                                 
                 "Regulation of Microtubule-Based Process",            
                 "Regulation of Canonical Wnt Signaling Pathway",      
                 "Nucleoside Triphosphate Metabolic Process",          
                 "TCF-Dependent Signaling in Response to Wnt",     
                 "ATP Metabolic Process",                              
                 "Regulation of Microtubule Cytoskeleton Organization",
                 "Canonical Wnt Signaling Pathway")
p8.df$direction <- ifelse(p8.df$NES > 0, "Upregulated", "Downregulated")
p8.df$database <- unlist(lapply(p8.df$NAME, function(x){ str_split(x, pattern = "_")[[1]][1] }))
p8.df$alpha <- ifelse(p8.df$`FDR q-val` < 0.05, 1, 0.5)
p8 <- ggplot(p8.df, aes(x = reorder(label, -NES), y = NES, fill = direction, pattern = database, alpha = alpha)) +
  geom_bar_pattern(stat = "identity", 
                   position = "dodge", 
                   pattern_angle = 45,
                   pattern_color = "black", 
                   pattern_fill = "black", 
                   color = "black", 
                   pattern_spacing = 0.05) +
  scale_fill_manual(values = c("Upregulated" = "#FF8C9E", "Downregulated" = "#40679E")) +
  scale_pattern_manual(values = c("REACTOME" = "stripe", "KEGG" = "none", "GOBP" = "circle")) +
  scale_alpha_identity() +
  labs(x = "Pathway", y = "NES", fill = "Direction", pattern = "Database") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(angle = 90, vjust = 0.5))
ggsave("~/Desktop/project/figures/plot2/gsea_suppa.png", p8, width = 5, height = 4.5, units = "in", dpi = 600)


gencode39 <- readGFF("~/Desktop/project/docs/data/0_references/gencode_v39.gtf")










