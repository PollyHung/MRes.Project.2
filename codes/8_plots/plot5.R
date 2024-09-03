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


setwd("~/Desktop/project/plots/8_plots/")
load("diff_surfProt.Rd")
load("mapping.Rd")
load("psi.Rd")
load("metatable.Rd")
load("cor_mtx.Rd")
load("mutation.Rd")
load("surface_proteins.Rd")

diff_surfProt <- separate(diff_surfProt, col = "AS_event", sep = ";", 
                          into = c("gene", "pb_id", "coordinates"), remove = FALSE)
diff_surfProt <- separate(diff_surfProt, col = "coordinates", sep = ":", 
                          into = c("AS_type", "coordinate"), extra = "merge")


p1.df <- data.frame(AS_event = diff_surfProt$AS_event,
                    C1_avg = rowMeans(psi.data[diff_surfProt$AS_event, condition_mtx$sample[condition_mtx$group == "group_1"]]),
                    C2_avg = rowMeans(psi.data[diff_surfProt$AS_event, condition_mtx$sample[condition_mtx$group == "group_2"]]),
                    dPSI = diff_surfProt$dpsi,
                    abs.dPSI = diff_surfProt$abs.dpsi,
                    FDR = diff_surfProt$p.adj)
p1.df <- p1.df %>% pivot_longer(cols = c(C1_avg, C2_avg), names_to = "Condition", values_to = "Avg_Value")

p1 <- ggplot(p1.df, aes(x = reorder(AS_event, Avg_Value), y = Avg_Value, color = Condition)) +
  geom_point(size = 2) +
  geom_line(aes(group = AS_event), color = "black", size = 0.5) +
  coord_flip() +
  labs(x = "Gene ID", y = "Average PSI") +
  theme_minimal(base_size = 7) + 
  theme(axis.title = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "in"))
ggsave("~/Desktop/project/figures/plot4/p1.png", p1, width = 5.5, height = 3.5, units = "in", dpi = 600)

MLSN <- data.frame(sample = names(psi.data["MSLN;PB.8038;A3:chr16:768379:768648:+", ]), 
                   psi = unname(psi.data["MSLN;PB.8038;A3:chr16:768379:768648:+", ]))
MLSN <- left_join(MLSN, condition_mtx)
ggplot(MLSN, aes(x = psi, fill = group)) +
  geom_density(alpha = 0.7) + 
  xlim(0, 1) +
  labs(x = "PSI", y = "Density", title = "Density Plot by Group") +
  theme_minimal() +
  scale_fill_manual(values = c("group_1" = "coral", "group_2" = "steelblue")) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_blank()
  )




temp <- readxl::read_xlsx("~/Desktop/temp.xlsx")
p1 <- ggplot(temp, aes(x = ))


copy_no.gene[grepl("SNRPA1|JUP|SNRPC|ZMAT5|PQBP1|WDR83|ISY1|PRPF40B|PRPF4B", rownames(copy_no.gene)), ]
copy_no.gene[grepl("SRSF2|SRPK3|PRPF40B|PPIL2|KHSRP|HNRNPA1|FRG1|EIF3A|CELF2|BCAS1", rownames(copy_no.gene)), ]












