# ---
# title: "AS proteins"
# author: "You"
# date: "2024-08-18"
# output: html_document
# ---

library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(DescTools)

setwd("/Volumes/Wild_Flower/OV_SpliceVariants")

# load files 
load("data/6_association/sample_data/psi_knn_imp.Rd")
load("data/1_sampleSheet/hh_ova_surv.Rd")
load("data/6_association/AS_regulators.Rd")
load("data/6_association/metatable.Rd")
load("data/4_mutect2/mutation.Rd")
load("data/6_association/lr_sr_wes_map.Rd")
load("data/6_association/tpm_matrix.Rd")
load("data/6_association/condition_mtx.Rd")
load("data/6_association/copy_number.Rd")


# Which splicing factor's expression differentiates between group 1 and group 2? 
# differential expression estimated by log2(mean(group.2)/mean(group.1))
SF_exp <- log_tpm_mtx[intersect(AS_genes, rownames(log_tpm_mtx)), ]
result.1 <- matrix(nrow = nrow(SF_exp), ncol = 3) %>% as.data.frame()
rownames(result.1) <- rownames(SF_exp)
colnames(result.1) <- c("log2FC", "p_value", "p_adj")

for(i in rownames(SF_exp)){
  group.1 <- SF_exp[i, condition_mtx$sample[condition_mtx$group == "group_1"]] %>% unlist()
  group.2 <- SF_exp[i, condition_mtx$sample[condition_mtx$group == "group_2"]] %>% unlist()
  
  test <- wilcox.test(group.1, group.2)
  e = 0.0001
  result.1[i, "log2FC"] <- log2((mean(group.2)+e)/(mean(group.1)+e))
  result.1[i, "p_value"] <- test$p.value
}

result.1 <- result.1 %>% dplyr::filter(abs(log2FC) > 0.5)
result.1$p_adj <- p.adjust(result.1$p_value, method = "fdr")
result.1 %>% dplyr::filter(p_adj < 0.25) %>% dplyr::arrange(p_adj)


# Which splicing factor's copy number differentiates between group 1 and group 2? 
# difference in copy number evaluated by mean(group.2) - mean(group.1)
# ∆CN > 0.5, threshold, so that group 2's copy number is 0.5 unit more 
copy_no.gene <- copy_no.gene[AS_genes, ] %>% na.omit()
result.2 <- matrix(nrow = nrow(copy_no.gene), ncol = 4) %>% as.data.frame()
rownames(result.2) <- rownames(copy_no.gene)
colnames(result.2) <- c("group.1.avg", "group.2.avg", "p_value", "p_adj")

for(i in rownames(copy_no.gene)){
  g1 <- intersect(colnames(copy_no.gene), condition_mtx$sample[condition_mtx$group == "group_1"])
  g2 <- intersect(colnames(copy_no.gene), condition_mtx$sample[condition_mtx$group == "group_2"])
  group.1 <- copy_no.gene[i, g1] %>% unlist()
  group.2 <- copy_no.gene[i, g2] %>% unlist()
  
  test <- wilcox.test(group.1, group.2)
  result.2[i, "group.1.avg"] <- mean(group.1)
  result.2[i, "group.2.avg"] <- mean(group.2)
  result.2[i, "p_value"] <- test$p.value
}

result.2$diff_cn <- result.2$group.2.avg - result.2$group.1.avg
result.2$abs.diff_cn <- abs(result.2$diff_cn)
result.2 <- result.2 %>% dplyr::filter(abs.diff_cn > 0.5) 
result.2$p_adj <- p.adjust(result.2$p_value, method = "fdr")
result.2 %>% dplyr::filter(p_adj < 0.25) %>% dplyr::arrange(p_adj)


# Which arm level copy number differentiates between group 1 and group 2?
result.3 <- matrix(nrow = nrow(copy_no.arm), ncol = 4) %>% as.data.frame()
rownames(result.3) <- rownames(copy_no.arm)
colnames(result.3) <- c("group.1.avg", "group.2.avg", "p_value", "p_adj")

for(i in rownames(copy_no.arm)){
  g1 <- intersect(colnames(copy_no.arm), condition_mtx$sample[condition_mtx$group == "group_1"])
  g2 <- intersect(colnames(copy_no.arm), condition_mtx$sample[condition_mtx$group == "group_2"])
  group.1 <- copy_no.arm[i, g1] %>% unlist()
  group.2 <- copy_no.arm[i, g2] %>% unlist()
  
  test <- wilcox.test(group.1, group.2)
  result.3[i, "group.1.avg"] <- mean(group.1)
  result.3[i, "group.2.avg"] <- mean(group.2)
  result.3[i, "p_value"] <- test$p.value
}

result.3$diff_cn <- result.3$group.2.avg - result.3$group.1.avg
result.3$abs.diff_cn <- abs(result.3$diff_cn)
result.3 <- result.3 %>% dplyr::filter(abs.diff_cn > 0.5)
result.3$p_adj <- p.adjust(result.3$p_value, method = "fdr")
result.3 %>% dplyr::arrange(p_adj)



# Which AS events does each SF genes regulate? 
result.4 <- matrix(nrow = nrow(psi_lr.imp), ncol = nrow(SF_exp)) %>% as.data.frame()
rownames(result.4) <- rownames(psi_lr.imp)
colnames(result.4) <- rownames(SF_exp)

for(SF in rownames(SF_exp)){
  gene_exp <- SF_exp[SF, condition_mtx$sample] %>% unlist()
  for(AS in rownames(psi_lr.imp)) {
    AS_psi <- psi_lr.imp[AS, condition_mtx$sample] %>% unlist()
    test <- cor.test(AS_psi, gene_exp)
    result.4[AS, SF] <- test$p.value
  }
}
result.4 <- apply(result.4, 2, function(x){ p.adjust(x, method = "fdr")} )
result.4 <- as.data.frame(result.4) %>% rownames_to_column(var = "AS_events")
result.4.long <- result.4 %>% 
  pivot_longer(cols = -AS_events, names_to = "SF_genes", values_to = "p_adjusted") %>% 
  dplyr::filter(p_adjusted < 0.25) %>% 
  dplyr::arrange(p_adjusted) #%>% 
  #dplyr::filter(SF_genes %in% rownames(result.2))

# other genes of interest: CCNE1, BRCA1, BRCA2
gene_exp <- log_tpm_mtx[c("CCNE1", "BRCA1", "BRCA2"), ]
result.5 <- matrix(nrow = nrow(gene_exp), ncol = 3) %>% as.data.frame()
rownames(result.5) <- rownames(gene_exp)
colnames(result.5) <- c("log2FC", "p_value", "p_adj")

for(i in rownames(gene_exp)){
  group.1 <- gene_exp[i, condition_mtx$sample[condition_mtx$group == "group_1"]] %>% unlist()
  group.2 <- gene_exp[i, condition_mtx$sample[condition_mtx$group == "group_2"]] %>% unlist()
  
  test <- wilcox.test(group.1, group.2)
  e = 0.0001
  result.5[i, "log2FC"] <- log2((mean(group.2)+e)/(mean(group.1)+e))
  result.5[i, "p_value"] <- test$p.value
}
result.5$p_adj <- p.adjust(result.5$p_value, method = "fdr")
# gene_exp$gene <- rownames(gene_exp)
# plot <- gene_exp %>% pivot_longer(cols = -gene, names_to = "sample", values_to = "gene_exp") 
# plot <- left_join(plot, condition_mtx)
# ggplot(plot, aes(x = group, y = gene_exp)) +
#   geom_boxplot() +
#   facet_wrap(~ gene, scales = "free_y", ncol = 1) + 
#   theme_bw()

# other genes: CCNE1, BRCA1, BRCA2
load("data/6_association/copy_number.Rd")
result.6 <- matrix(nrow = 3, ncol = 4) %>% as.data.frame()
rownames(result.6) <- c("CCNE1", "BRCA1", "BRCA2")
colnames(result.6) <- c("group.1.avg", "group.2.avg", "p_value", "p_adj")

for(i in c("CCNE1", "BRCA1", "BRCA2")){
  g1 <- intersect(colnames(copy_no.gene), condition_mtx$sample[condition_mtx$group == "group_1"])
  g2 <- intersect(colnames(copy_no.gene), condition_mtx$sample[condition_mtx$group == "group_2"])
  group.1 <- copy_no.gene[i, g1] %>% unlist()
  group.2 <- copy_no.gene[i, g2] %>% unlist()
  
  test <- wilcox.test(group.1, group.2)
  result.6[i, "group.1.avg"] <- mean(group.1)
  result.6[i, "group.2.avg"] <- mean(group.2)
  result.6[i, "p_value"] <- test$p.value
}
result.6$p_adj <- p.adjust(result.6$p_value, method = "fdr")
# copy_no <- copy_no.gene[c("CCNE1", "BRCA1", "BRCA2"), 2:ncol(copy_no.gene)]
# copy_no$gene <- rownames(copy_no)
# copy_no <- copy_no %>% pivot_longer(cols = -gene, names_to = "sample", values_to = "copy_number")
# copy_no <- left_join(copy_no, condition_mtx)
# ggplot(copy_no, aes(x = group, y = copy_number)) +
#   geom_boxplot() +
#   facet_wrap(~ gene, scales = "free_y", ncol = 1) + 
#   theme_bw()


# Compare proportion of mutated samples in group 1 vs group 2 ===================
# dplyr::filter(grepl("exonic|splicing", Func.ensGene)) ## filter by exonic mutations 
# three mutation groups: "exonic|splicing", "UTR|upstream|downstream", and "intronic|intergenic"
SF_mut <- mutations %>% 
  dplyr::filter(Gene.ensGene %in% AS_genes) %>% 
  dplyr::select(Chr, Ref, Alt, Func.ensGene, Gene.ensGene, SIFT_score, sample) 
SF_mut <- merge(SF_mut, metatable, by.x = "sample", by.y = "copy_number")

temp <- metatable %>% dplyr::filter(copy_number %in% unique(SF_mut$sample))
filter <- c("exonic|splicing", "UTR|upstream|downstream", "intronic|intergenic")

# fisher exact test to test proportion 
result.7 <- matrix(nrow = length(unique(SF_mut$Gene.ensGene)), ncol = 4) %>% as.data.frame()
colnames(result.7) <- c("exonic|splicing", "UTR|upstream|downstream", "intronic|intergenic", "all_mutation")
rownames(result.7) <- unique(SF_mut$Gene.ensGene)

for(SF in unique(SF_mut$Gene.ensGene)){
  for(filt in filter){
    mut_gene <- SF_mut %>% 
      dplyr::filter(Gene.ensGene == SF) %>% 
      dplyr::filter(grepl(filt, Func.ensGene))
    temp[temp$copy_number %in% mut_gene$sample, "mut"] <- 1
    temp[!temp$copy_number %in% mut_gene$sample, "mut"] <- 0
    mut_contingency <- table(temp$group, temp$mut)
    if(length(unique(temp[["mut"]])) == 1){
      print("no difference between group 1 and 2")
    } else {
      test <- fisher.test(mut_contingency)
      result.7[SF, filt] <- test$p.value
    }
  }
  
  mut_gene <- SF_mut %>% dplyr::filter(Gene.ensGene == SF)
  temp[temp$copy_number %in% mut_gene$sample, "mut"] <- 1
  temp[!temp$copy_number %in% mut_gene$sample, "mut"] <- 0
  mut_contingency <- table(temp$group, temp$mut)
  if(length(unique(temp[["mut"]])) == 1){
    print("no difference between group 1 and 2")
  } else {
    test <- fisher.test(mut_contingency)
    result.7[SF, "all_mutation"] <- test$p.value
  }
}

result.7.adj <- apply(result.7, 2, function(x){ p.adjust(x, method = "fdr") })


# Differentially expressed AS events between group 1 and group 2? 
result.8 <- matrix(nrow = nrow(psi_lr.imp), ncol = 3) %>% as.data.frame()
rownames(result.8) <- rownames(psi_lr.imp) 
colnames(result.8) <- c("dPSI", "log2FC", "p_value")

for(AS in rownames(psi_lr.imp)){
  group.1 <- psi_lr.imp[AS, condition_mtx$sample[condition_mtx$group == "group_1"]]
  group.2 <- psi_lr.imp[AS, condition_mtx$sample[condition_mtx$group == "group_2"]]
  
  test <- wilcox.test(group.1, group.2)
  result.8[AS, "p_value"] <- test$p.value
  result.8[AS, "dPSI"] <- mean(group.2) - mean(group.1)
  result.8[AS, "log2FC"] <- log2((mean(group.2)+e)/(mean(group.1)+e)) #mean(group.2) - mean(group.1)
}
result.8.adj <- result.8 %>% dplyr::filter(abs(dPSI) > 0.20)
result.8.adj$p_adj <- p.adjust(result.8.adj$p_value, method = "fdr")

result.8.adj %>% dplyr::filter(p_adj < 0.25) 
result.8 <- result.8 %>% dplyr::arrange(desc(dPSI))
gene <- data.frame(gene <- sub(";.*", "", rownames(result.8)), 
                   dPSI <- result.8$dPSI) 
write.table(gene, "data/6_association/gsea.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# MSLN alone 
MSLN <- psi_lr.imp[grepl("MSLN", rownames(psi_lr.imp)), ]
result.9 <- matrix(nrow = nrow(MSLN), ncol = 2) %>% as.data.frame()
rownames(result.9) <- rownames(rownames(MSLN))
colnames(result.9) <- c("dPSI", "p_value")

for(AS in rownames(MSLN)){
  group.1 <- MSLN[AS, condition_mtx$sample[condition_mtx$group == "group_1"]]
  group.2 <- MSLN[AS, condition_mtx$sample[condition_mtx$group == "group_2"]]
  
  test <- wilcox.test(group.1, group.2)
  result.9[AS, "p_value"] <- test$p.value
  result.9[AS, "dPSI"] <- mean(group.2)-mean(group.1)
}


# GSEA 
source("/Volumes/Wild_Flower/Psoriasis_ImmuneFeatures/skin/codes/scripts/scRNA_functions.R")
reactome <- rbind(read.delim("data/7_gsea/reactome/gsea_report_for_na_neg_1724075321737.tsv", check.names = F), 
                  read.delim("data/7_gsea/reactome/gsea_report_for_na_pos_1724075321737.tsv", check.names = F))
reactome <- reactome[, 1:11]
p1 <- gsea_barplot(gsea_df = reactome, FDR = 0.25, top.n = 20, title = "log2FC PSI Reactome")
ggsave(filename = "plot/plot_1/log2FC PSI reactome.png", plot = p1, width = 3.25, height = 2.21, units = "in", dpi = 600)

kegg <- rbind(read.delim("data/7_gsea/kegg/gsea_report_for_na_neg_1724075489398.tsv", check.names = F), 
              read.delim("data/7_gsea/kegg/gsea_report_for_na_pos_1724075489398.tsv", check.names = F))
kegg <- kegg[, 1:11]
p2 <- gsea_barplot(gsea_df = kegg, FDR = 0.25, top.n = 20, title = "log2FC PSI Kegg")
ggsave(filename = "plot/plot_1/log2FC PSI kegg.png", plot = p2, width = 2.2, height = 1, units = "in", dpi = 600)

gobp <- rbind(read.delim("data/7_gsea/gobp/gsea_report_for_na_neg_1724075503250.tsv", check.names = F), 
              read.delim("data/7_gsea/gobp/gsea_report_for_na_pos_1724075503250.tsv", check.names = F))
gobp <- gobp[, 1:11]
p3 <- gsea_barplot(gsea_df = gobp, FDR = 0.25, top.n = 20, title = "log2FC PSI GoBP")
ggsave(filename = "plot/plot_1/log2FC PSI gobp.png", plot = p3, width = 3.55, height = 3.5, units = "in", dpi = 600)


# dichotimise chr arm 
result.10 <- matrix(nrow = nrow(copy_no.arm), ncol = 2) %>% as.data.frame()
rownames(result.10) <- rownames(copy_no.arm)
colnames(result.10) <- c("p_value", "p_adj")

copy_no.arm.di <- ifelse(copy_no.arm > 0.5, "amp", "wt")
copy_no.arm.di <- as.data.frame(t(copy_no.arm.di))
copy_no.arm.di$sample <- rownames(copy_no.arm.di)
copy_no.arm.di <- left_join(copy_no.arm.di, condition_mtx, by="sample")

for(i in rownames(copy_no.arm)){
  if(length(unique(copy_no.arm.di[, i])) == 1){
    print("no difference between group 1 and 2")
  } else {
    contigency_mtx <- table(copy_no.arm.di["group"], copy_no.arm.di[i])
    test <- fisher.test(contigency_mtx)
    result.10[i, "p_value"] <- test$p.value
  }
}
copy_no.arm <- as.data.frame(t(copy_no.arm))
copy_no.arm$sample <- rownames(copy_no.arm)
copy_no.arm <- left_join(copy_no.arm, condition_mtx, by="sample")

