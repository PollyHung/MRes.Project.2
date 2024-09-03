

setwd("/Volumes/Wild_Flower/OV_SpliceVariants/data/")
load("6_association/sample_data/psi_knn_imp.Rd")
load("6_clustering/output/clust_lr.RData")
mapping <- read.table("1_sampleSheet/mapping_complete.txt", header = T, sep = "\t")
mapping <- mapping %>% dplyr::filter(long_read %in% colnames(psi_lr.imp)) %>% 
  dplyr::select("long_read", "copy_number")
copy_number <- read.table("/Volumes/Wild_Flower/OV_CopyNumber/docs/HH_ova/gistic/input_2/broad_values_by_arm.txt", 
                          sep = "\t", row.names = 1) %>% t() %>% as.data.frame() %>% 
  dplyr::rename(copy_number = `Chromosome Arm`)
copy_number <- merge(copy_number, mapping, by="copy_number")
rownames(copy_number) <- copy_number$long_read
copy_number$long_read <- NULL
copy_number$copy_number <- NULL
psi <- psi_lr.imp[clust$ts_clusters$Junction, rownames(copy_number)]

p_value <- matrix(NA, nrow = nrow(clust$ts_clusters), ncol = ncol(copy_number)) %>% as.data.frame()
rownames(p_value) <- clust$ts_clusters$Junction
colnames(p_value) <- colnames(copy_number)

cor_val <- matrix(NA, nrow = nrow(clust$ts_clusters), ncol = ncol(copy_number)) %>% as.data.frame()
rownames(cor_val) <- clust$ts_clusters$Junction
colnames(cor_val) <- colnames(copy_number)

sample_id <- intersect(rownames(copy_number), clust$cluster_exp$Sample)

for(i in colnames(copy_number)){
  cn <- copy_number[sample_id, i] %>% unlist() %>% as.numeric()
  for(j in rownames(psi)){
    psi_val <- psi[j, sample_id] %>% unlist()
    test <- cor.test(cn, psi_val, method = "pearson")
    p_value[j, i] <- test$p.value
    cor_val[j, i] <- test$estimate
  }
}

p_adjust <- apply(p_value, 2, function(x){p.adjust(x, method = "fdr")})

cor_pval <- apply(p_adjust, 2, function(x){ifelse(x < 0.10, "sig", "ns")}) %>% as.data.frame()
cor_pval <- as.data.frame(mapply(function(x, y) ifelse(x == "sig", y, x), cor_pval, cor_val))
cor_pval <- cor_pval[, colSums(cor_pval != "ns") > 0]

save(list = c("p_value", "cor_val", "p_adjust", "cor_pval"), 
     file = "6_association/copynumber_association.Rd")









