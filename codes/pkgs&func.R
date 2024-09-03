library(NOISeq)
library(maSigPro)
library(edgeR)
library(DESeq2)
library(DEXSeq)
library(goseq)
library(ggplot2)
library(VennDiagram)
library(MASS)
library(plyr)
library(ggrepel)
library(seqinr)
library(dplyr)
library(magrittr)
library(rtracklayer)
library(Biostrings)
library(doParallel)
library(parallel)
library(rhmmer)
library(RColorBrewer)
library(ggsci)
library(scales)
library(GenomicRanges)
library(tidyr)
library(purrr)
library(data.table)
library(DRIMSeq)
library(IsoformSwitchAnalyzeR)
library(Biostrings)
library(DGEobj.utils)
library(survival)
library(pheatmap)
library(survminer)
library(tximport)
library(limma)
library(ape)
library(biomaRt)
library(MCPcounter)
library(easier)
library(genekitr)
library(mclust)
library(FNN)
library(caret)
library(kernlab)


## Functions: 
plot_heatmap <- function(matrix, 
                         plot_dir = "plot/clustering/", 
                         plot_name){
  custom_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  png(file.path(plot_dir, plot_name), width = 6, height = 6, units = "in", res = 600)
  heatmap.2(cor_matrix, margins = c(9, 9), key = TRUE, key.title = "Color Key", 
            key.xlab = "Correlation", key.ylab = "", trace = "none", 
            density.info = "none", col = custom_palette,
            breaks = seq(0, 1, length.out = 101))
  dev.off()
}

plot_dendrogram <- function(matrix, 
                            plot_dir = "plot/clustering/", 
                            plot_name){
  png(file.path(plot_dir, plot_name), width = 6, height = 6, units = "in", res = 600)
  plot(h.meta)
  dev.off()
}

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}
