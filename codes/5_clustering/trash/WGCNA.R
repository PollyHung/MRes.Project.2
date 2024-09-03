#---
# name: Cluster Isoforms Based on Their Expression Profile 
# date: 2024-July-31
#---

library(caret)
library(magrittr)
library(DEXSeq)
library(dynamicTreeCut)
library(gplots)
library(WGCNA)
library(survival)
library(MCPcounter)
set.seed(50)

## [1] tpm normalized isoform expression matrix ================================
isoform_expressions <- read.csv("data/3_flair/pacbio/normalized_tpm.csv", 
                                check.names = FALSE, row.names = 1)

## [2] calculate correlation matrix ============================================
if(!file.exists("data/8_clustering/isoform_tpm_cor.Rd")){
  cor_matrix <- cor(isoform_expressions, method = "pearson")
  dist.meta <- stats::as.dist(1 - cor_matrix)
  h.meta <- stats::hclust(dist.meta, method = "complete")
  
  plot_heatmap(matrix = cor_matrix, plot_name = "01_correlation_heatmap.png")
  plot_dendrogram(matrix = h.meta, plot_name = "02_dendrogram_by_distmtx.png")
  
  save(list = c("cor_matrix", "dist.meta", "h.meta"), file = "data/8_clustering/isoform_tpm_cor.Rd")
}

## [3] Find variable features (feature selection prior to WGCNA)? ==============
rf_df <- as.data.frame(t(isoform_expressions))
rf_df$labels <- ifelse(grepl("FN", rownames(rf_df)), "Normal", "Tumour")

if(!file.exists("data/8_clustering/rf_featureSelection.rds")){
  trainControl <- trainControl(method = "cv", 
                               number = 10, 
                               verboseIter = TRUE)
  model <- train(labels ~ ., 
                 data = rf_df, 
                 method = "rf", 
                 trControl = trainControl, 
                 importance = TRUE)
  feature_imp <- varImp(model, scale = FALSE)
  imp_df <- as.data.frame(feature_imp$importance)
  imp_isoforms <- imp_df %>% dplyr::arrange(desc(Normal)) %>% dplyr::filter(Normal >= 1, Tumour >= 1)
  
  save(list = c("model", "feature_imp", "imp_df", "imp_isoforms"), file = "data/8_clustering/rf_featureSelection.rds")
}

if(!file.exists("data/8_clustering/isoform_featSelect_cor.Rd")){
  subset_df <- isoform_expressions[rownames(imp_isoforms), ]
  cor_matrix <- cor(subset_df, method = "pearson")
  dist.meta <- stats::as.dist(1 - cor_matrix)
  h.meta <- stats::hclust(dist.meta, method = "complete")
  
  plot_heatmap(matrix = cor_matrix, plot_name = "03_correlation_heatmap_featSelect.png")
  plot_dendrogram(matrix = h.meta, plot_name = "04_dendrogram_by_distmtx_featSelect.png")
  
  save(list = c("subset_df", "cor_matrix", "dist.meta", "h.meta"), file = "data/8_clustering/isoform_featSelect_cor.Rd")
}

## [4] try clustering isoforms =================================================
if(!file.exists("data/8_clustering/isoform_clustering.Rd")){
  cor_matrix <- cor(t(subset_df), method = "pearson")
  dist.meta <- stats::as.dist(1 - cor_matrix)
  h.meta <- stats::hclust(dist.meta, method = "complete")
  
  plot_heatmap(matrix = cor_matrix, plot_name = "05_correlation_heatmap_featSelect_isoforms.png")
  plot_dendrogram(matrix = h.meta, plot_name = "06_dendrogram_by_distmtx_featSelect_isoforms.png")
  plot_heatmap(matrix = as.matrix(t(subset_df)), plot_name = "07_correlation_heatmap_isoforms_by_sample.png")
  
  save(list = c("cor_matrix", "dist.meta", "h.meta"), file = "data/8_clustering/isoform_clustering.Rd")
}

## [5] use WGCNA to identify clusters ==========================================
# Rationale: https://smorabit.github.io/hdWGCNA/articles/isoform_pbmcs.html#isoform-co-expression-network-analysis
if(!file.exists("data/8_clustering/wgcna.rds")){
  input_mat = t(subset_df)
  normal_samples <- rownames(input_mat)[grepl("FN", rownames(input_mat))]
  reorder_samples <- c(normal_samples, setdiff(rownames(input_mat), normal_samples))
  input_mat <- input_mat[reorder_samples, ]
  #allowWGCNAThreads(nThreads = 60)
  
  # pick soft threshold, signed network 
  powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
  sft <- pickSoftThreshold(input_mat, 
                           powerVector = powers, 
                           corFnc = cor, 
                           networkType = "signed")
  #png("~/MRes.project.2/plots/02_clustering/08_soft_thresholding.png", width = 7, height = 4, units = "in", res = 600)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[ ,1], 
       -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2", 
       type = "n", 
       main = paste("Scale independence"), 
       ylim = c(0, 1))
  text(sft$fitIndices[ ,1], -sign(sft$fitIndices[ ,3])*sft$fitIndices[, 2], labels = powers, cex = 0.9, col = "red")
  abline(h = c(0.6, 0.7, 0.8), col = "red")
  plot(sft$fitIndices[, 1], 
       sft$fitIndices[, 5], 
       xlab = "Soft Threshold (power)", 
       ylab = "Mean connectivity",
       type = "n", 
       main = paste("Mean connectivity"))
  text(sft$fitIndices[ ,1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
  #dev.off()
  
  net = blockwiseModules(input_mat, 
                         power = 14,
                         TOMType = "signed", 
                         minModuleSize = 20,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  table(net$colors)
  mergedColors = labels2colors(net$colors)
  plotDendroAndColors(net$dendrograms[[1]],
                      mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05)
  module_df <- data.frame(gene_id = names(net$colors),
                          colors = labels2colors(net$colors))
  save(list = c("input_mat", "powers", "sft", "net", "module_df", "mergedColors"), 
       file = "data/8_clustering/wgcna.rds")
}

