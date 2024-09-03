setwd("/Volumes/Wild_Flower/OV_SpliceVariants/data/")
load("6_association/sample_data/immune_mtx.Rd")
load("6_clustering/output/clust_lr.RData")
load("6_association/sample_data/psi_knn_imp.Rd")

## association with immune? 
p_value <- matrix(NA, nrow = nrow(clust$ts_clusters), ncol = ncol(immune_mtx)) %>% as.data.frame()
rownames(p_value) <- clust$ts_clusters$Junction
colnames(p_value) <- colnames(immune_mtx)

cor_val <- matrix(NA, nrow = nrow(clust$ts_clusters), ncol = ncol(immune_mtx)) %>% as.data.frame()
rownames(cor_val) <- clust$ts_clusters$Junction
colnames(cor_val) <- colnames(immune_mtx)

sample_id <- intersect(rownames(immune_mtx), clust$cluster_exp$Sample)

for(i in colnames(immune_mtx)){
  immune <- immune_mtx[sample_id, i] %>% unlist()
  for(j in clust$ts_clusters$Junction){
    PSI <- psi_lr.imp[j, sample_id] %>% unlist()
    test <- cor.test(PSI, immune, method = "pearson")
    p_value[j, i] <- test$p.value
    cor_val[j, i] <- test$estimate
  }
}

p_adjust <- apply(p_value, 2, function(x){p.adjust(x, method = "fdr")})

cor_val <- apply(p_adjust, 2, function(x){ifelse(x < 0.10, "sig", "ns")}) %>% as.data.frame()
cor_pval <- as.data.frame(mapply(function(x, y) ifelse(x == "sig", y, x), cor_pval, cor_val))
cor_pval$`Endothelial cells` <- NULL
cor_pval$Fibroblasts <- NULL

save(list = c("cor_pval", "p_adjust", "p_value", "cor_val"), 
     file = "6_association/immune_association.Rd")
