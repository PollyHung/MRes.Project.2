setwd("/Volumes/Wild_Flower/OV_SpliceVariants/data/")

load("6_clustering/output/clust_lr.RData")
load("6_association/survival_association2.Rd")
clust$ts_clusters <- merge(clust$ts_clusters, survival, by.x = "Junction", by.y = "AS_event")

load("6_association/immune_association.Rd")
cor_pval$Junction <- rownames(cor_val)
clust$ts_clusters <- merge(clust$ts_clusters, cor_pval, by = "Junction")

load("6_association/copynumber_association.Rd")
cor_pval$Junction <- rownames(cor_val)
clust$ts_clusters <- merge(clust$ts_clusters, cor_pval, by = "Junction")

load("4_mutect2/mutect2_combined.Rd")
load("6_association/sample_data/immune_mtx.Rd")















