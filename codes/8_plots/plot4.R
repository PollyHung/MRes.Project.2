library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(magrittr)
library(ggtranscript)
library(rtracklayer)

colours = c("#FF8C9E", "#EF9C66", "#FFDA76", "#9CA777", "#41B3A2", "#378CE7", "#40679E", "#667BC6", "#D7C3F1", "#944E63")
color.heavy = c("#E07F80", "#95AAD3", "#F0CF7F", "#C195C4", "#94C47D", "#8785BA", "#E69965", "#72BEB7")
color.light = c("#EDAEAE", "#B9DBF4", "#F4E6AA", "#F0D8E9", "#C8DBC8", "#BAB4D8", "#F6BF93", "#B8E2DE")

setwd("~/Desktop/project/plots/8_plots/")
source("~/Desktop/project/codes/code/8_plots/function.R")
load("mapping.Rd")
load("psi.Rd")
load("metatable.Rd")
load("cor_mtx.Rd")
load("copy_number.Rd")
load("immune_mtx.Rd")
load("gsea.Rd")
load("tpm_lr.Rd")
load("SF_genes.Rd")
load("mutation.Rd")
load("surface_proteins.Rd")


## Plot 1 Differentially Expressed Alternative Splicing Events between Clusters 
# dPSI.filter %>% dplyr::filter(p.adj < 0.25)
psi_long <- as.data.frame(psi.data) %>% 
  dplyr::mutate(AS_event = rownames(psi.data)) %>% 
  pivot_longer(cols = `0700055A`:`T18-158-FT1`, names_to = "sample", values_to = "PSI")
psi_long$gene_id <- unlist(lapply(psi_long$AS_event, function(x){str_split(x, ";")[[1]][1]}))
psi_long <- left_join(psi_long, condition_mtx, by = "sample")
psi_long$cluster <- ifelse(psi_long$group == "group_1", "C1", "C2")

p1.df <- psi_long %>% dplyr::filter(AS_event == "MSLN;PB.8038;A3:chr16:768379:768648:+")
p1 <- ggplot(p1.df, aes(y = PSI, x = cluster)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.5, fill = "#667BC6") + 
  theme_bw(base_size = 7) + 
  theme(axis.title = element_text(size = 6))

p2.df <- psi_long %>% dplyr::filter(AS_event == "MSLN;PB.8038;RI:chr16:768565-768648:+")
p2 <- ggplot(p2.df, aes(y = PSI, x = cluster)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.5, fill = "#667BC6") + 
  theme_bw(base_size = 7) + 
  theme(axis.title = element_text(size = 6))


## Plot 2: Isoforms 
setwd("~/Desktop/project/docs/data/")
classification <- read.table("2_sqanti/filter_rule_RulesFilter_result_classification.txt", header = T)
gtf <- readGFF("2_sqanti/filter_rule.filtered.gtf")
gtf <- merge(gtf, classification[, c("isoform", "structural_category", "associated_gene", "subcategory")], 
             by.x = "transcript_id", by.y = "isoform")

MSLN_exons <- gtf %>% 
  dplyr::filter(type == "exon", seqid == "chr16", gene_id == "PB.8038") %>% 
  dplyr::mutate(label_name = gsub("PB.8038", "MSLN", transcript_id)) %>% 
  dplyr::arrange(structural_category) %>% 
  dplyr::mutate(label_ordered = factor(paste(structural_category, label_name),
                                levels = unique(paste(structural_category, label_name)))) %>% 
  dplyr::mutate(label_ordered = fct_rev(label_ordered)) %>% 
  dplyr::mutate(structural_category = gsub("_", " ", structural_category)) 

p2 <- ggplot(MSLN_exons, aes(xstart = start,
                             xend = end,
                             y = label_ordered)) +
  geom_range(aes(fill = structural_category)) +
  geom_intron(data = to_intron(MSLN_exons, "transcript_id"),
              aes(strand = strand), 
              arrow.min.intron.length = 300) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 13), 
        legend.title = element_text(size = 13), 
        legend.key.size = unit(0.2, "in"))
ggsave("~/Desktop/project/figures/plot4/p2.png", p2, width = 7, height = 5, units = "in", dpi = 600) 



# ## MSLN =========================================================================================================
setwd("~/Desktop/project/docs/data/")
sample_list <- list.files("1_preprocess/")
MSLN_list <- list()
for(sample_id in sample_list){
  gff <- readGFF(file.path("1_preprocess", sample_id,"sqanti_output/filter_rule.filtered.gtf"))
  MSLN <- dplyr::filter(gff, seqid == "chr16", start > 759178, end < 770797)
  if(nrow(MSLN) != 0){
    MSLN$sample <- sample_id
    MSLN_list[[sample_id]] <- MSLN
  }
}

MSLN <- do.call(rbind, MSLN_list)
MSLN$sample <- gsub("_", "-", MSLN$sample)
MSLN <- left_join(MSLN, condition_mtx, by = "sample") %>%
  dplyr::select(seqid, source, type, start, end, strand, gene_id, transcript_id, sample, group) %>%
  dplyr::filter(type == "exon")

MSLN_C1 <- MSLN %>% dplyr::filter(type == "exon", group == "group_1")
MSLN_C2 <- MSLN %>% dplyr::filter(type == "exon", group == "group_2")

p3 <- ggplot(MSLN_C1, aes(xstart = start, xend = end, y = sample)) +
  geom_range(fill = "#FF8C9E") +
  geom_intron(data = to_intron(MSLN_C1, "sample"), aes(strand = strand),
              arrow.min.intron.length = 300) +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 13), legend.title = element_text(size = 13),
        legend.key.size = unit(0.2, "in"))

p4 <- ggplot(MSLN_C2, aes(xstart = start, xend = end, y = sample)) +
  geom_range(fill = "#378CE7") +
  geom_intron(data = to_intron(MSLN_C2, "sample"), aes(strand = strand),
              arrow.min.intron.length = 300) +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 13), legend.title = element_text(size = 13),
        legend.key.size = unit(0.2, "in"))

reference <- readGFF("~/Desktop/project/docs/data/0_references/gencode_v39.gtf")
reference.MSLN <- reference %>%
  dplyr::filter(hgnc_id == "HGNC:7371") %>%
  dplyr::filter(type == "exon") %>%
  dplyr::select(seqid, source, type, start, end, score, strand, phase, gene_id, transcript_id) %>%
  # dplyr::rename(gene_id = hgnc_id) %>%
  # dplyr::mutate(sample = "gencode v39") %>%
  dplyr::mutate(group = "reference") # %>% 
  # dplyr::filter(start < 768379)

p5 <- ggplot(reference.MSLN, aes(xstart = start,
                                 xend = end,
                                 y = gene_id)) +
  geom_range(fill = "lightgrey") +
  geom_intron(data = to_intron(reference.MSLN, "gene_id"),
              aes(strand = strand),
              arrow.min.intron.length = 300) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.2, "in"))

MSLN_A3 <- readGFF("~/Desktop/project/docs/data/5_differential_expression/DAS/long_reads/MSLN_A3.gtf")
MSLN_RI <- readGFF("~/Desktop/project/docs/data/5_differential_expression/DAS/long_reads/MSLN_RI.gtf")
# MSLN_A3$name <- "MSLN.A3"
# MSLN_RI$name <- "MSLN.RI"
# MSLN_A3 <- MSLN_A3 %>% 
#   dplyr::select(seqid, source, type, start, end, strand, gene_id, transcript_id) %>%
#   dplyr::mutate(sample = gene_id) %>%
#   dplyr::mutate(group = "AS_event") 
# MSLN_RI <- MSLN_RI %>%
#   dplyr::select(seqid, source, type, start, end, strand, gene_id, transcript_id) %>%
#   dplyr::mutate(sample = gene_id) %>%
#   dplyr::mutate(group = "AS_event") 

MSLN.reference <- reference %>%
  dplyr::filter(hgnc_id == "HGNC:7371") %>%
  dplyr::filter(type == "exon") %>%
  dplyr::select(seqid, source, type, start, end, score, strand, phase, gene_id, transcript_id) %>%
  dplyr::filter(start < 768379)

MSLN_A3 <- rbind(MSLN.reference, MSLN_A3)
MSLN_RI <- rbind(MSLN.reference, MSLN_RI)
MSLN_RI$group <- "intron retention"
MSLN_RI$sample <- "RI:chr16:768565-768648:+"
MSLN_A3$group <- "alternative 3'"
MSLN_A3$sample <- "A3:chr16:768379:768648:+"
  
MSLN_mut <- rbind(MSLN_A3, MSLN_RI)

# p5 <- ggplot(MSLN_mut, aes(xstart = start,
#                            xend = end,
#                            y = sample)) +
#   geom_range(aes(fill = sample)) +
#   geom_intron(data = to_intron(MSLN_mut, "sample"),
#               aes(strand = strand),
#               arrow.min.intron.length = 300) +
#   theme_minimal(base_size = 14) +
#   theme(axis.title = element_text(size = 13),
#         legend.title = element_text(size = 13),
#         legend.key.size = unit(0.2, "in"))
# MSLN_plot <- rbind(reference.MSLN, MSLN_mut, MSLN)
MSLN_mut$group <- factor(MSLN_mut$group, levels = c("intron retention", "alternative 3'"))
MSLN_mut
p6 <- ggplot(MSLN_mut, aes(xstart = start, xend = end, y = sample)) +
  geom_range(aes(fill = group)) +
  geom_intron(data = to_intron(MSLN_mut, "sample"),
              aes(strand = strand),
              arrow.min.intron.length = 300) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.2, "in")) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +  # Adjust facet heights based on the number of rows
  coord_cartesian(clip = "off") +
  theme(panel.spacing.y = unit(0.5, "lines"),  # Adjust spacing between facets if needed
        strip.text.y = element_text(angle = 0))
# ggsave("~/Desktop/project/figures/plot4/p3.png", p6, width = 12, height = 30, units = "in", dpi = 600)

# MSLN_A3.simp <- MSLN_A3 %>% dplyr::select(seqid, start, end, strand, gene_id, transcript_id)
# 
# write.table(MSLN_A3.simp, file = "~/Desktop/MSLN_A3.simp.gff", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(MSLN_RI, file = "~/Desktop/MSLN_RI.gff", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 








MSLN_A3 <- readGFF("~/Desktop/project/docs/data/5_differential_expression/DAS/long_reads/MSLN_A3.gtf")
MSLN_RI <- readGFF("~/Desktop/project/docs/data/5_differential_expression/DAS/long_reads/MSLN_RI.gtf")
MSLN_A3.A1 <- MSLN_A3[1, ]
MSLN_A3.A2 <- MSLN_A3[2, ]
MSLN_RI.A1 <- MSLN_RI[1, ]
MSLN_RI.A2 <- MSLN_RI[2, ]
MSLN_RI.A3 <- MSLN_RI[3, ]

MSLN.reference <- reference %>%
  dplyr::filter(hgnc_id == "HGNC:7371") %>%
  dplyr::filter(type == "exon") %>%
  dplyr::select(seqid, source, type, start, end, score, strand, phase, gene_id, transcript_id) %>%
  dplyr::filter(start < 768379)

MSLN_A3.A1 <- rbind(MSLN.reference, MSLN_A3.A1)
MSLN_A3.A2 <- rbind(MSLN.reference, MSLN_A3.A2)
MSLN_RI.A1 <- rbind(MSLN.reference, MSLN_RI.A1)
MSLN_RI.A2 <- rbind(MSLN.reference, MSLN_RI.A2)
MSLN_RI.A3 <- rbind(MSLN.reference, MSLN_RI.A3)

MSLN_A3.A1$group <- "A3:chr16:768379:768648:+:alternative1"
MSLN_A3.A2$group <- "A3:chr16:768379:768648:+:alternative2"
MSLN_RI.A1$group <- "RI:chr16:768565-768648:+:alternative2"
MSLN_RI.A2$group <- "RI:chr16:768565-768648:+:alternative2"
MSLN_RI.A3$group <- "RI:chr16:768565-768648:+:alternative1"

MSLN_mut <- rbind(MSLN_A3.A1, MSLN_A3.A2, MSLN_RI.A1, MSLN_RI.A2, MSLN_RI.A3)

# p5 <- ggplot(MSLN_mut, aes(xstart = start,
#                            xend = end,
#                            y = sample)) +
#   geom_range(aes(fill = sample)) +
#   geom_intron(data = to_intron(MSLN_mut, "sample"),
#               aes(strand = strand),
#               arrow.min.intron.length = 300) +
#   theme_minimal(base_size = 14) +
#   theme(axis.title = element_text(size = 13),
#         legend.title = element_text(size = 13),
#         legend.key.size = unit(0.2, "in"))
# MSLN_plot <- rbind(reference.MSLN, MSLN_mut, MSLN)
MSLN_mut$group <- factor(MSLN_mut$group, levels = c("intron retention", "alternative 3'"))
MSLN_mut <- rbind(reference.MSLN, MSLN_mut)
p6 <- ggplot(MSLN_mut, aes(xstart = start, xend = end, y = group)) +
  geom_range(aes(fill = group)) +
  geom_intron(data = to_intron(MSLN_mut, "group"),
              aes(strand = strand),
              arrow.min.intron.length = 300) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.2, "in")) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +  # Adjust facet heights based on the number of rows
  coord_cartesian(clip = "off") +
  theme(panel.spacing.y = unit(0.5, "lines"),  # Adjust spacing between facets if needed
        strip.text.y = element_text(angle = 0))
# ggsave("~/Desktop/project/figures/plot4/p3.png", p6, width = 12, height = 30, units = "in", dpi = 600)

# MSLN_A3.simp <- MSLN_A3 %>% dplyr::select(seqid, start, end, strand, gene_id, transcript_id)
MSLN_mut.simp <- MSLN_mut %>% dplyr::select(seqid, source, type, start, end, score, strand, phase, group)
colnames(MSLN_mut.simp) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "group")
MSLN_mut.simp$score <- "."
MSLN_mut.simp$frame <- "."

write.table(MSLN_mut.simp, file = "~/Desktop/MSLN_mut.gff", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(MSLN_RI, file = "~/Desktop/MSLN_RI.gff", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 



















