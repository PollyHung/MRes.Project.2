library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(magrittr)

colours = c("#FF8C9E", "#EF9C66", "#FFDA76", "#9CA777", "#41B3A2", "#378CE7", "#40679E", "#667BC6", "#D7C3F1", "#944E63")
color.heavy = c("#E07F80", "#95AAD3", "#F0CF7F", "#C195C4", "#94C47D", "#8785BA", "#E69965", "#72BEB7")
color.light = c("#EDAEAE", "#B9DBF4", "#F4E6AA", "#F0D8E9", "#C8DBC8", "#BAB4D8", "#F6BF93", "#B8E2DE")


load("data/8_plots/mapping.Rd")
load("data/8_plots/colours.Rd")
load("data/8_plots/plot.1/IsoformFiltered.Rd")
load("data/8_plots/plot.1/ShortRead_TPM.Rd")
load("data/8_plots/plot.1/LongRead_TPM.Rd")
load("data/8_plots/plot.1/kallisto_tpm.Rd")
load("data/8_plots/plot.1/psi.Rd")

## Plot 1 ----------------------------------------------------------------------
p1.df <- sqanti_filtered %>% 
  dplyr::count(label) %>% 
  dplyr::mutate(percentage = n / sum(n) * 100) %>% 
  dplyr::mutate(label = paste0(label, " (", round(percentage, digits = 2), "%)"))
p1 <- ggplot(p1.df, aes(x = "", y = percentage, fill = label)) +
  geom_bar(width = 0.25, stat = "identity", color = "white", linewidth = 0.1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colours) +
  theme_void(base_size = 7) + 
  theme(legend.title = element_blank(), 
        legend.key.size = unit(0.1, units = "in"))

## Plot 2 ----------------------------------------------------------------------
p2.df <- sqanti_filtered %>% 
  dplyr::group_by(label) %>% 
  dplyr::count(subcategory) %>% 
  dplyr::mutate(percentage = n / sum(n) * 100) %>% 
  dplyr::mutate(subcategory = gsub("_", " ", subcategory))
p2.df <- p2.df %>% 
  dplyr::mutate(label_sub = case_when(subcategory == "3prime fragment" ~ "3' fragment", 
                          subcategory == "5prime fragment" ~ "5' fragment", 
                          subcategory == "alternative 3end5end" ~ "alt 3'5'", 
                          subcategory == "alternative 5end" ~ "alt 5'", 
                          subcategory == "at least one novel splicesite" ~ "≥1 novel splicesite", 
                          subcategory == "combination of known junctions" ~ "merged junctions", 
                          subcategory == "combination of known splicesites" ~ "merged splicesites", 
                          TRUE ~ subcategory))
p2 <- ggplot(p2.df, aes(x = label, y = percentage, fill = label_sub)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.05) +
  scale_fill_manual(values = colours) +
  labs(x = NULL, y = "Percentage") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "right", 
        legend.direction = "horizontal", 
        legend.key.size = unit(0.1, units = "in")) + 
  guides(fill = guide_legend(ncol = 2))

## Plot 3 ----------------------------------------------------------------------
sr_df <- short_reads %>% dplyr::select(any_of(map$short_read))
sr_cov <- ifelse(sr_df == 0, 0, 1) %>% as.data.frame()
sr_cov$isoform <- rownames(sr_cov)
sr_cov <- merge(sr_cov, sqanti_filtered[, c("isoform", "label")], by = "isoform")
p3.df <- sr_cov %>%
  pivot_longer(cols = starts_with("G72"), names_to = "sample", values_to = "detected") %>%
  group_by(label, sample) %>%
  summarise(percentage_detected = mean(detected) * 100)
p3 <- ggplot(p3.df, aes(x = label, y = percentage_detected, fill = label)) +
  geom_violin(trim = FALSE, linewidth = 0.2) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.2) +
  scale_fill_manual(values = colours) + # Adjust colors as needed
  labs(x = "Structural Category", 
       y = "Isoforms (%)") +
  theme_bw(base_size = 7) + ylim(0, 100) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 6)) # Remove legend if not needed

## Plot 4 ----------------------------------------------------------------------
p4.df <- ifelse(long_reads == 0, 0, 1) %>% as.data.frame()
p4.df$isoform <- rownames(p4.df)
p4.df <- merge(p4.df, sqanti_filtered[, c("isoform", "label")], by = "isoform")
p4.df <- p4.df %>% pivot_longer(cols = colnames(p4.df)[2:26], names_to = "sample", values_to = "detected") %>%
  group_by(label, sample) %>%
  summarise(percentage_detected = mean(detected) * 100)
p4 <- ggplot(p4.df, aes(x = label, y = percentage_detected, fill = label)) +
  geom_violin(trim = FALSE, linewidth = 0.2) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.2) +
  scale_fill_manual(values = colours) + # Adjust colors as needed
  labs(x = "Structural Category", 
       y = "Isoforms (%)") +
  theme_bw(base_size = 7) + ylim(0, 100) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 6)) # Remove legend if not needed

## Plot 5 ----------------------------------------------------------------------
p5.df <- psi.metadata %>% dplyr::group_by(AS_type) %>%
  summarise(gene_count = n_distinct(gene), pacbio_count = n_distinct(rownames)) %>%
  pivot_longer(cols = c(gene_count, pacbio_count), names_to = "Type", values_to = "Count")
p5 <- ggplot(p5.df, aes(x = AS_type, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("pacbio_count" = "#FF8C66", "gene_count" = "#7F8C66"), 
                    labels = c("pacbio_count" = "AS event", "gene_count" = "Gene")) +
  labs(x = "Types of AS", y = "Counts of cases") +
  theme_bw(base_size = 7) + 
  ylim(0, 550) + 
  theme(axis.title = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.1, units = "in")) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 1.8)

## Plot 6 ----------------------------------------------------------------------
psi_sr <- read.table("data/5_differential_expression/DAS/short_reads/psi_event.psi")
psi_sr[is.na(psi_sr)] <- 0
p6.df <- data.frame(AS_events = rownames(psi_sr), 
                    mean_PSI = rowMeans(psi_sr))
p6 <- ggplot(p5.df, aes(x = AS_type, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("pacbio_count" = "#FF8C66", "gene_count" = "#7F8C66"), 
                    labels = c("pacbio_count" = "AS event", "gene_count" = "Gene")) +
  labs(x = "Types of AS", y = "Counts of cases") +
  theme_bw(base_size = 7) + 
  ylim(0, 550) + 
  theme(axis.title = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.1, units = "in")) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 1.8)




























## Save Plots 
# ggsave("plot/plot.1/p1.png", p1, width = 1.5, height = 1, units = "in", dpi = 600)
# ggsave("plot/plot.1/p2.png", p2, width = 2.8, height = 1.1, units = "in", dpi = 600)
# ggsave("plot/plot.1/p3.png", p3, width = 1, height = 1.1, units = "in", dpi = 600)
# ggsave("plot/plot.1/p4.png", p4, width = 1, height = 1.1, units = "in", dpi = 600)
# ggsave("plot/plot.1/p5.png", p5, width = 2.5, height = 1.1, units = "in", dpi = 600)