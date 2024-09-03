# ---
# title: "Clustering Patients"
# author: "You"
# date: "2024-08-10"
# output: html_document
# ---

set.seed(20020208)

## methods: 
# NMF clustering 
# hierarchical 
# k means 

library(pheatmap)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stats)
library(cluster)
library(caret)
library(reshape2)
library(ConsensusClusterPlus)

## hierarchical ======================
setwd("/Volumes/Wild_Flower/OV_SpliceVariants/data/")
load("6_association/sample_data/psi_knn_imp.Rd")

cor_mtx <- cor(psi_lr.imp, method = "spearman")
heatmap(cor_mtx, margins = c(10, 10))

group_1 <- c("HH12000124FN2", "HH12000157FN4", "T16-088-FT2", "T16-002FN1", 
             "T16-106-FT1", "T15-162-FT2", "T15-036-FN2", "T14-042FN3", 
             "T15-022-FT2", "090061A", "T16-035-FN2", "10-149A", "T15-051-FT4")
group_2 <- setdiff(colnames(psi_lr.imp), group_1)
condition_mtx <- data.frame(
  sample = c(group_1, group_2), 
  group = c(rep("group_1", length(group_1)), rep("group_2", length(group_2))), 
  condition = ifelse(grepl("FN", c(group_1, group_2)), "normal", "tumour"))
rownames(condition_mtx) <- condition_mtx$sample
print(condition_mtx)
normal <- condition_mtx$sample[grepl("FN", condition_mtx$sample)]
tumour <- condition_mtx$sample[!grepl("FN", condition_mtx$sample)]


## Jaccard pairwise clustering 
threshold <- 0  # Example threshold, adjust according to your data
binary_data <- psi.data > threshold
binary_data <- as.data.frame(binary_data)  # Convert to a data frame if necessary
jaccard_sim_matrix <- as.matrix(vegdist(t(binary_data), method = "jaccard", binary = TRUE))
hc <- hclust(as.dist(1 - jaccard_sim_matrix), method = "complete")
plot(hc)
clusters <- cutree(hc, k = 2)
heatmap(jaccard_sim_matrix, Rowv = as.dendrogram(hc), symm = TRUE)


## AS-based clustering was associated with clinical characteristics and immune features
# setwd("6_clustering/")
# results = ConsensusClusterPlus(psi_sr.imp, 
#                                maxK=3, 
#                                reps=10, 
#                                pItem=0.8, 
#                                pFeature=1, 
#                                title="consensusClustering", 
#                                clusterAlg="hc", 
#                                distance="spearman",
#                                seed=55,
#                                plot="png", 
#                                writeTable = TRUE, 
#                                verbose = TRUE)

# # here we can see that there is roughly 2 clusters happening by raw clustering
# 
# ## K-means clustering ======================
# # row = samples, columns = AS events 
# melt_psi <- melt(psi_lr.imp) %>% dplyr::rename(Event = Var1, sample = Var2, PSI = value)
# median <- melt_psi %>% group_by(sample) %>% summarize(median_PSI = median(PSI, na.rm = TRUE))
# melt_psi <- melt_psi %>% mutate(sample = factor(sample, levels = median$sample[order(median$median_PSI)]))
# 
# ggplot(melt_psi, aes(x = sample, y = PSI)) +
#   geom_boxplot() + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x = "Sample", y = "PSI Value") + 
#   coord_flip()
# 
# # normalise and scale 
# arcsine_transform <- function(p) {
#   asin(sqrt(p))
# }
# psi <- apply(psi_lr.imp, 2, arcsine_transform)
# 
# # k-means clustering 
# set.seed(20020208)  
# ## optimal number of clusters? -----------------------
# # elbow plot
# wcss <- sapply(2:10, function(k){
#   kmeans_result <- kmeans(as.data.frame(psi_lr.imp), centers = k, nstart = 25)
#   kmeans_result$tot.withinss
# })
# plot(2:10, wcss, type = "b", pch = 19, frame = FALSE,
#      xlab = "Number of Clusters (k)",
#      ylab = "Total Within-Cluster Sum of Squares",
#      main = "Elbow Plot", 
#      ylim = c(2440, 3560))
# 
# # silhouette method 
# avg_silhouette <- sapply(2:10, function(k){
#   km_res <- kmeans(as.data.frame(psi_lr.imp), centers = k, nstart = 25)
#   ss <- silhouette(km_res$cluster, dist(psi_lr.imp))
#   mean(ss[, 3])
# })
# plot(2:10, avg_silhouette, type = "b", pch = 19, frame = FALSE,
#      xlab = "Number of Clusters (k)",
#      ylab = "Average Silhouette Width",
#      main = "Silhouette Plot", 
#      ylim = c(0.1, 0.5))
# 
# # gap statistics 
# gap_stat <- clusGap(as.data.frame(psi_lr.imp), FUN = kmeans, nstart = 25, K.max = 10, B = 50)
# plot(gap_stat, main = "Gap Statistic")
# 
# # between_ss/total_ss score 
# betweenss_ratio <- sapply(2:10, function(k){
#   kmeans_result <- kmeans(as.data.frame(psi_lr.imp), centers = k, nstart = 25)
#   kmeans_result$betweenss / kmeans_result$totss
# })
# plot(2:10, betweenss_ratio, type = "b", pch = 19, frame = FALSE,
#      xlab = "Number of Clusters (k)",
#      ylab = "Between_SS / Total_SS Ratio",
#      main = "Between_SS / Total_SS Ratio", 
#      ylim = c(0.5, 0.7))
# 
# ## 4 cluster seems the best 
# k_clust <- kmeans(as.data.frame(psi_lr.imp), centers = 4)
# table(k_clust$cluster)
# #   1   2   3   4 
# # 535 627 427 473 
# 
# subset_events <- names(k_clust$cluster)[k_clust$cluster == 4]
# subset_mtx <- cor(psi_lr.imp[variant_AS, ], method = "spearman")
# heatmap(subset_mtx, margins = c(10, 10))
# 
# 
# ## remove less variance rows =====================
# row_variances <- apply(psi_lr.imp, 1, var, na.rm = TRUE)
# variant_AS <- names(row_variances)[row_variances > quantile(unname(row_variances))[[2]]]
# subset_mtx <- cor(psi_lr.imp[variant_AS, ], method = "spearman")
# heatmap(subset_mtx, margins = c(10, 10))
# 
# k_clust <- kmeans(as.data.frame(t(psi_lr.imp[variant_AS, ])), centers = 2)
# condition_mtx <- cbind(condition_mtx[names(k_clust$cluster), ], 
#                        k_cluster = unname(k_clust$cluster))
# 
# 
# ## try distance matrix ============================
# hc <- hclust(dist(as.data.frame(t(psi_lr.imp)), method = "euclidean"), method = "complete")
# plot(hc)
# gmm <- Mclust(t(psi_lr.imp))
# plot(gmm$BIC)
# summary(gmm)
# 
# 
# ## rowsum? ===========================
# upper_names <- colnames(psi_lr.imp)[colMeans(psi_lr.imp) > median.default(colMeans(psi_lr.imp))]
# upper <- rep("upper", length(upper_names))
# names(upper) <- upper_names
# lower_names <- colnames(psi_lr.imp)[colMeans(psi_lr.imp) <= median.default(colMeans(psi_lr.imp))]
# lower <- rep("lower", length(lower_names))
# names(lower) <- lower_names
# byMedian <- c(upper, lower)
# byMedian <- byMedian[condition_mtx$sample]
# 
# condition_mtx$byMedian <- unname(byMedian)
# 
# 
