# Assuming you have a variable `total_samples_lr` and `total_samples_sr` representing the total number of samples
total_samples_lr <- ncol(psi_lr)  # Replace with the actual total number of long-read samples
total_samples_sr <- ncol(psi_sr)  # Replace with the actual total number of short-read samples

# Create bins
psi_na <- psi_na %>%
  mutate(
    lr_bin = case_when(
      na_in_lr / total_samples_lr > 0.25 & na_in_lr / total_samples_lr <= 0.50 ~ "25~50% missing",
      na_in_lr / total_samples_lr > 0.50 & na_in_lr / total_samples_lr <= 0.75 ~ "50~75% missing",
      na_in_lr / total_samples_lr > 0.75 ~ "> 75% missing", 
      TRUE ~ as.character(na_in_lr)
    ),
    sr_bin = case_when(
      na_in_sr / total_samples_sr > 0.25 & na_in_sr / total_samples_sr <= 0.50 ~ "25~50% missing",
      na_in_sr / total_samples_sr > 0.50 & na_in_sr / total_samples_sr <= 0.75 ~ "50~75% missing",
      na_in_sr / total_samples_sr > 0.75 ~ "> 75% missing", 
      TRUE ~ as.character(na_in_lr)
    )
  )
level_lr <- c(0:7, "25~50% missing", "50~75% missing", "> 75% missing") %>% as.character()
level_sr <- c(0:25, "25~50% missing", "50~75% missing", "> 75% missing") %>% as.character()

psi_na <- psi_na %>%
  mutate(
    lr_bin = factor(lr_bin, levels = level_lr),
    sr_bin = factor(sr_bin, levels = level_sr)
  )

p1 <- ggplot(psi_na, aes(x = lr_bin)) +
  geom_bar() +
  labs(title = "NA in Long Reads",
       x = "% of missing PSI",
       y = "number of isoforms") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        title = element_text(size = 10))
ggsave(plot = p1, filename = "/Volumes/Wild_Flower/OV_SpliceVariants/plot/other_stuff/NA_in_lr.png", 
       width = 3, height = 3, dpi = 600)

p2 <- ggplot(psi_na, aes(x = sr_bin)) +
  geom_bar() +
  labs(title = "NA in Short Reads",
       x = "% of missing PSI",
       y = "number of isoforms") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        title = element_text(size = 10))
ggsave(plot = p2, filename = "/Volumes/Wild_Flower/OV_SpliceVariants/plot/other_stuff/NA_in_sr.png", 
       width = 4.5, height = 3, dpi = 600)

cor_value <- data.frame(AS_event = names(diag(cor_result)), 
                        cor_value = diag(cor_result))
cor_value2 <- data.frame(sample_id = rownames(cor_result2), 
                         cor_value = diag(cor_result2))
p3 <- ggplot(cor_value, aes(x = cor_value)) + 
  geom_density() + geom_vline(xintercept = 0.6) + 
  labs(x = "correlation value", 
       title = "Density Plot of Correlation Between 
  Long-Read and Short-Read PSI Values 
  for the Same Alternative Splicing Events")
ggsave(plot = p3, filename = "/Volumes/Wild_Flower/OV_SpliceVariants/plot/other_stuff/cor_ASevent.png", 
       width = 4, height = 3, dpi = 600)

p4 <- ggplot(cor_value2, aes(x = cor_value)) + 
  geom_density() + 
  labs(x = "correlation value", 
       title = "Density Plot of Correlation Between 
  Long-Read and Short-Read PSI Values 
  for the Same Sample") 
ggsave(plot = p4, filename = "/Volumes/Wild_Flower/OV_SpliceVariants/plot/other_stuff/cor_sample.png", 
       width = 4, height = 3, dpi = 600)

temp <- clust$ts_clusters
p5 <- ggplot(temp, aes(x = Counts)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +  # Adjust binwidth as needed
  stat_bin(geom = "text", aes(label = ..count..), vjust = -0.5, binwidth = 1, color = "black") +
  scale_x_continuous(breaks = unique(temp$Counts)) +
  ylim(c(0, 30)) + 
  labs(title = "Distribution of Counts",
       x = "Counts",
       y = "Frequency") +
  theme_bw()
ggsave(plot = p5, filename = "/Volumes/Wild_Flower/OV_SpliceVariants/plot/other_stuff/dist_of_counts.png", 
       width = 4, height = 3, dpi = 600)


