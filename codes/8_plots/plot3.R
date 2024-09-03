library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(magrittr)
library(survminer)
library(survival)
library(karyoploteR)
library(rtracklayer)

colours = c("#FF8C9E", "#EF9C66", "#FFDA76", "#9CA777", "#41B3A2", "#378CE7", "#40679E", "#667BC6", "#D7C3F1", "#944E63")
color.heavy = c("#E07F80", "#95AAD3", "#F0CF7F", "#C195C4", "#94C47D", "#8785BA", "#E69965", "#72BEB7")
color.light = c("#EDAEAE", "#B9DBF4", "#F4E6AA", "#F0D8E9", "#C8DBC8", "#BAB4D8", "#F6BF93", "#B8E2DE")

source("~/Desktop/project/codes/code/8_plots/function.R")
setwd("~/Desktop/project/plots/8_plots/")
load("mapping.Rd")
load("psi.Rd")
load("metatable.Rd")
load("cor_mtx.Rd")
load("copy_number.Rd")
load("immune_mtx.Rd")
load("gsea.Rd")
load("tpm_lr.Rd")
load("diff_SF.Rd")
load("SF_genes.Rd")
load("SF_families.Rd")
load("mutation.Rd")

condition_mtx$cluster <- ifelse(condition_mtx$group == "group_1", "C1", "C2")
metatable$cluster <- ifelse(metatable$group == "group_1", "C1", "C2")

cluster <- condition_mtx$sample

# ## Plot 1: overall splicing factor's expression --------------------------------
p1.df <- log_tpm_mtx[intersect(SF, rownames(log_tpm_mtx)), ]
p1.df$gene <- rownames(p1.df)
p1.df <- p1.df %>% pivot_longer(cols = `0700055A`:`T18-158-FT1`, names_to = "sample", values_to = "expression")
p1.df <- left_join(p1.df, condition_mtx, by = "sample")

wilcox.test(p1.df$expression[p1.df$group == "group_1"], p1.df$expression[p1.df$group == "group_2"])

p1 <- ggplot(p1.df, aes(x = expression, y = cluster))+
  geom_boxplot(fill = "#40679E", color = "black", linewidth = 0.25, outlier.size = 0.25) +
  coord_flip() +
  theme_bw(base_size = 7) +
  labs(x = "SF expression", y = "") +
  theme(axis.title = element_text(size = 6))
# ggsave("plot/plot.3/p1.png", p1, width = 1, height = 1, units = "in", dpi = 600)

p2 <- ggplot(p1.df, aes(x = expression, fill = cluster)) +
  geom_density(alpha = 0.5, color = "black", linewidth = 0.25) +
  theme_bw(base_size = 7) +
  labs(x = "SF expression", y = "Density") +
  theme(axis.title = element_text(size = 6), legend.title = element_blank(),
        legend.key.size = unit(0.07, "in"))
# ggsave("plot/plot.3/p2.png", p2, width = 1.6, height = 1, units = "in", dpi = 600)

## Plot 2: group by splicing factor family and check again --------------------------------
diff.exp.SF.family$SF_fam = rownames(diff.exp.SF.family)
SF_fam <- rownames(filter(diff.exp.SF.family, p_adj < 0.20))
SF_fam <- SF_families[SF_fam]
SF_fam <- do.call(rbind, lapply(names(SF_fam), function(fam) {
  data.frame(SF_fam = fam, gene = SF_fam[[fam]], stringsAsFactors = FALSE)
}))
p2.df <- p1.df %>% dplyr::filter(gene %in% SF_fam$gene)
p2.df <- left_join(p2.df, SF_fam, by = "gene")
p2.df$label <- gsub("_Family", "", p2.df$SF_fam)
p2.df <- left_join(p2.df, diff.exp.SF.family, by = "SF_fam")
p2.df <- p2.df %>% dplyr::mutate(label = forcats::fct_reorder(label, p_adj, .fun = min))


p3 <- ggplot(p2.df, aes(x = expression, y = cluster)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.2, fill = "#667BC6", color = "black") + 
  coord_flip() + 
  labs(y = "", x = "Copy Number") + 
  theme_minimal(base_size = 7) + 
  theme(axis.title.y = element_text(size = 6)) + 
  facet_wrap(~ label, scales = "free_x", nrow = 1)
ggsave("plot/plot.3/p3.png", p3, width = 2.8, height = 1.2, units = "in", dpi = 600)


## Plot 3: Differential Expressed Splicing Factors between Cluster 1 and Cluster 2 ---------------------------
p3.df <- diff.exp.SF %>% 
  mutate(significance = case_when(p_adj < 0.25 & abs(log2FC) >= 0.5 ~ "Significant", TRUE ~ "Not Significant"))

p4 <- ggplot(p3.df, aes(x = log2FC, y = -log10(p_adj), color = significance)) +
  geom_point(alpha = 0.8, size = 1) +
  scale_color_manual(values = c("Significant" = "#E07F80", "Not Significant" = "lightgrey")) +
  geom_vline(xintercept = c(-0.5, 0.5), 
             linetype = "dashed", 
             color = "grey") +
  geom_hline(yintercept = -log10(0.25), 
             linetype = "dashed", 
             color = "grey") +
  labs(x = "log2(Fold Change)", 
       y = "-log10(p-value)", 
       color = "Significance") +
  theme_bw(base_size = 9) + 
  theme(legend.position = "none", 
        axis.title = element_text(size = 8)) + 
  geom_text_repel(data = subset(p3.df, significance == "Significant"), 
                  aes(label = rownames(subset(p3.df, significance == "Significant"))), 
                  max.overlaps = 10, 
                  size = 2.5, 
                  color = "black") 
ggsave("~/Desktop/project/figures/plot3/volcano.png", p4, units = "in", width = 2.4, height = 2.1, dpi = 600)



## Plot 4: Are any of these gene mutated between Cluster 1 and Cluster 2? ---------------------------
## We selected for genes with log2 fold change above 0.5 to continue 
gene <- c("SRSF2", "SRPK3", "PRPF40B", "PPIL2", "KHSRP", "HNRNPA1", "FRG1", "EIF3A", "CELF2", "BCAS1")
temp <- filter(diff.mut.SF, abs(pct_diff) > 40)
temp$p_adj <- p.adjust(temp$p_value, method = "fdr")
SF_mut <- rownames(temp)
p5.df <- mutations %>% dplyr::filter(gene %in% SF_mut) %>% 
  dplyr::select(gene, long_read, group) %>% 
  dplyr::rename(sample = long_read) %>% 
  dplyr::mutate(group = ifelse(group == "group_1", "C1", "C2")) %>% 
  dplyr::mutate(mutation = 1) %>% 
  distinct() %>% 
  pivot_wider(names_from = gene, values_from = mutation, values_fill = 0) %>% 
  pivot_longer(cols = -c(sample, group), names_to = "gene", values_to = "mutation") %>% 
  dplyr::mutate(label = ifelse(mutation == 1, "mut", "wt"))
p5.df$sample <- factor(p5.df$sample, levels = intersect(cluster, unique(p5.df$sample)))
p5.df <- p5.df %>% dplyr::arrange(sample)

p5 <- ggplot(p5.df, aes(x = gene, y = sample, fill = factor(label))) + 
  geom_tile(color = "white") + 
  coord_flip() +
  scale_fill_manual(values = c("wt" = "#40679E", "mut" = "#FF8C9E", "no data" = "white"), name = "Mutation") +
  labs(y = "Sample") +
  theme_bw(base_size = 6) + 
  theme(axis.title = element_text(size = 6), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.key.size = unit(0.07, "in"), 
        legend.title = element_text(size = 6))
ggsave("~/Desktop/project/figures/plot3/mutation_tile.png", p5, units = "in", width = 2.4, height = 1.6, dpi = 600)


## Plot 5: Are any of these gene's copy number different between Cluster 1 and Cluster 2? ---------------------------
temp <- filter(diff.cn.SF, abs.diff_cn > 0.5)
temp$p_adj <- p.adjust(temp$p_value, method = "fdr")
temp <- merge(temp, copy_no.gene[rownames(temp), ], by = "row.names")

temp <- temp %>% pivot_longer(cols = `0700055A`:`T18-158-FT1`, names_to = "sample", values_to = "copy_number")
temp <- left_join(temp, condition_mtx, by = "sample")

p6.df <- temp %>% dplyr::filter(grepl("3q", Cytoband))
p6 <- ggplot(p6.df, aes(x = copy_number, y = cluster)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.2, fill = "#41B3A2", color = "black") + 
  coord_flip() + 
  labs(y = "", x = "Copy Number") + 
  theme_bw(base_size = 7) + 
  theme(axis.title.y = element_text(size = 6)) + 
  facet_wrap(~ Row.names, scales = "free_x", nrow = 1)
# ggsave("plot/plot.3/p6.png", p6, width = 2.5, height = 1, units = "in", dpi = 600)



temp <- diff.cn.SF %>% dplyr::filter(abs.diff_cn > 0.5)
temp$p_adj <- p.adjust(temp$p_value)


