folder = "figure2"


grid <- clust$ts_clusters %>% 
  dplyr::filter(os_prognosis != "ns" | pfs_prognosis != "ns") 
grid <- separate(grid, col = "Junction", sep = ";", into = c("pacbio_id", "others"), remove = F, extra = "merge")
grid <- separate(grid, col = "others", sep = ":", into = c("AS_event", "chr_location"), extra = "merge")
grid$strand <- sub(".*:", "", grid$chr_location)
grid$label <- paste(paste(grid$Gene, grid$AS_event, sep= ":"), paste0("n=",grid$Counts), sep = " ")


grid$log_pval <- -log10(grid$Wilcox_pval)
p1 <- ggplot(grid, aes(x = reorder(label, -log_pval), y = log_pval, fill = deltaPSI_TvsN)) +
  geom_bar(stat = "identity", linewidth = 0.05) + # Use geom_bar to create bars with height based on log_pval
  scale_fill_gradient2(low = "#604CC3", mid = "white", high = "#FFB534", name = "∆PSI") + # Adjust colors as needed
  labs(x = "AS event", y = expression(-log[10](p-value)), title = "Survival association and prevalence") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                          axis.text.y = element_text(angle = 90, hjust = 1))
ggsave(file.path("/Volumes/Wild_Flower/OV_SpliceVariants/plot", folder, "p1.png"), p1, 
       width = 2.3, height = 6, units = "in")

grid <- grid %>% dplyr::arrange(-log_pval)

plot_df <- grid %>% pivot_longer(cols = c(os_prognosis, pfs_prognosis), 
                                 names_to = "Prognosis_Type", 
                                 values_to = "Value")
plot_df$Prognosis_Type <- gsub("_prognosis", "", plot_df$Prognosis_Type)
p2 <- ggplot(plot_df, aes(x = reorder(label, -log_pval), y = Prognosis_Type, fill = factor(Value))) +
  geom_tile(color = "white", linewidth = 0.4) + theme_minimal(base_size = 12) + coord_fixed() + #coord_flip() + 
  scale_fill_manual(values = c("favourable" = "#FF4E88", "unfavourable" = "#536493", "ns" = "grey")) +
  labs(x = "alternative splicing event", y = "", fill = "prognosis") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path("/Volumes/Wild_Flower/OV_SpliceVariants/plot", folder, "p2.png"), p2, 
       width = 3, height = 7, units = "in")


grid <- grid[, colSums(grid != "ns") > 0]
columns <- intersect(colnames(immune_mtx), colnames(grid))
plot_df <- grid %>% pivot_longer(cols = columns, 
                                 names_to = "immune_type", 
                                 values_to = "Value")
plot_df$Value <- ifelse(plot_df$Value == "ns", NA, as.numeric(plot_df$Value))
p3 <- ggplot(plot_df, aes(x = reorder(label, -log_pval), y = immune_type, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.4) + 
  theme_minimal(base_size = 12) + 
  coord_fixed() + 
  scale_fill_gradient2(low = "#004225", mid = "white", high = "#E3651D", midpoint = 0, na.value = "grey") + # Gradient for numeric values
  labs(x = "Alternative Splicing Event", y = "", fill = "immune") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path("/Volumes/Wild_Flower/OV_SpliceVariants/plot", folder, "p3.png"), p3, 
       width = 3.6, height = 7, units = "in") 


columns <- intersect(rownames(copy_number), colnames(grid))
plot_df <- grid %>% pivot_longer(cols = columns, 
                                 names_to = "chromosome_arms", 
                                 values_to = "Value")
plot_df$Value <- ifelse(plot_df$Value == "ns", NA, as.numeric(plot_df$Value))
p4 <- ggplot(plot_df, aes(x = reorder(label, -log_pval), y = chromosome_arms, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.4) + 
  theme_minimal(base_size = 12) + 
  coord_fixed() + 
  scale_fill_gradient(low = "#378CE7", mid = "white", high = "#FF4E88", midpoint = 0, na.value = "grey") + # Gradient for numeric values
  labs(x = "Alternative Splicing Event", y = "", fill = "Copy Number") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path("/Volumes/Wild_Flower/OV_SpliceVariants/plot", folder, "p4.png"), p4, 
       width = 8, height = 8.5, units = "in") 









