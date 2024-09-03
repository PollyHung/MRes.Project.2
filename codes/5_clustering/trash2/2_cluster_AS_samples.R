# ---
# title: "Clustering Patients Based on AS Events"
# author: "You"
# date: "2024-08-10"
# output: html_document
# ---

## Library packages 
library(rtracklayer)
library(dplyr)
library(tidyverse)
library(mclust)
require(doParallel)
require(parallel)

### Clustering Patients by AS Expression: --------------------------------------
# HGSOC cancer can be classified into different subtypes based on gene expression 
#   and alternative splicing (AS) levels.
# What Was Clustered?: Clustering samples, the clustering was based on the 
#   expression patterns of AS events.       
# Simultaneous Grouping: The Gaussian Mixture Model (GMM) simultaneously grouped 
#   tumor and normal samples based on these AS expression patterns. This means they 
#   used the expression levels of AS events to determine how similar or different 
#   the samples were, leading to the formation of clusters of samples (patients).           
## Identification of Tumor-Specific AS Events:       
# Tumor Subpopulations: Within each cluster (sample subpopulation), they identified 
#   AS events with two key features:        
#     1. High Frequency of Tumor Samples: The cluster was enriched with tumor samples.        
#     2. Significant Differential Splicing: The AS events showed significant differential 
#        splicing (ΔPSI) compared to normal tissues.     
# Tumor-Specific AS Events: By examining the AS events within these tumor-enriched 
#   clusters, they identified 3,059 AS events in breast cancer that had a |ΔPSI| ≥ 20% 
#   in subpopulations of at least 50 patients. This means these AS events were 
#   consistently altered in a significant portion of the patient subpopulation.        


setwd("/Volumes/Wild_Flower/OV_SpliceVariants/data/")

gtf <- import.gff("2_sqanti/filter_rule.filtered.gtf") 
classification <- read.table("2_sqanti/filter_rule_RulesFilter_result_classification.txt", header = T)
classification$gene_id <- sub("\\.[^.]*$", "", classification$isoform)
load("6_association/sample_data/psi_knn_imp.Rd")    
psi_lr.imp <- as.data.frame(psi_lr.imp)
psi_lr.imp$rownames <- rownames(psi_lr.imp)
psi_lr.imp <- separate(psi_lr.imp, col = "rownames", into = c("pacbio_id", "events"), sep = ";")
psi_lr.imp <- separate(psi_lr.imp, col = "events", into = c("AS_event", "coordinate"), sep = ":", extra = "merge")
sample_id <- colnames(psi_lr.imp)[1:25]

## select which event? 
event = c("SE", "A5", "A3", "AF", "AL", "RI", "MX")
event_file = "5_differential_expression/DAS/long_reads/sqanti_filter.events.ioe"
clust_dir <- "6_clustering/clust_output/"

## isolate event files 
psi <- psi_lr.imp %>% filter(AS_event %in% event) %>% dplyr::select(all_of(sample_id))
events <- read.table(file = event_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
events <- events %>% dplyr::filter(event_id %in% rownames(psi))

## add gene names 
idx.gene <- match(events$gene_id, classification$gene_id)
events$gene_symbol <- classification$associated_gene[idx.gene]

## separate tumour and normal 
psi_tumor <- psi[ , !grepl("FN", sample_id)]
psi_control <- psi[ , grepl("FN", sample_id)]

## Sort PSI_tumor and PSI_control using order from events
idx.events <- match(events$event_id, rownames(psi_tumor))
psi_tumor <- psi_tumor[idx.events, ]
psi_control <- psi_control[idx.events, ]

#Parameters for clustering
min_cluster = 10 # clusters below this threshold are merged
min_tumor = 5 # min number of tumor samples to be clustered after removing NA samples
min_controls = 2 # min number of control samples to be clustered after removing NA

## create directory if not exist 
if(dir.exists(clust_dir)){
  unlink(clust_dir, force = T, recursive = T)
}
dir.create(clust_dir)

## remove events based on coverage and number of samples 
flags <- data.frame(Junction.id = events$event_id,
                    Remove = rep(FALSE, nrow(events)))

for(i in 1:nrow(events)){
  se.id <- events$event_id[i]
  gene <- events$gene_symbol[i]
  
  idx.case <- which(!is.na(psi_tumor[se.id, ]))
  idx.control <- which(!is.na(psi_control[se.id, ]))
  
  exp.case <- psi_tumor[se.id, idx.case, drop = FALSE]
  exp.control <- psi_control[se.id, idx.control, drop = FALSE]
  
  ## if number of the event is 0
  if(length(exp.case) < min_tumor || length(exp.control) < min_controls){
    flags$Remove[i] <- TRUE
    next
  }
  ## if standardGeneric of the event is 0
  if(sd(exp.case)==0 || sd(exp.control)==0){
    flags$Remove[i] <- TRUE
    next
  }
}

## events for GMM clustering 
id.keep <- flags$Junction.id[!flags$Remove]
events_filt <- events %>% filter(event_id %in% id.keep)
cat("INFO: Events to be clustered: ", nrow(events_filt), "\n")

cluster_DF <- matrix(NA, nrow = 1, ncol = 25)
colnames(cluster_DF) <- c(colnames(psi_tumor), colnames(psi_control))
rownames(cluster_DF) <- "skip"

## run mixture in parallel 
nmodels <- nrow(events_filt) # number of models to fit 
for(i in 1:nmodels){
  se.id <- events_filt$event_id[i]
  gene <- events_filt$gene_symbol[i]
  
  cat("Processing ", i, " ", se.id, "\n")
  
  idx.case <- which(!is.na(psi_tumor[se.id, ]))
  idx.control <- which(!is.na(psi_control[se.id, ]))
  
  exp.case <- psi_tumor[se.id, idx.case, drop = FALSE]
  exp.control <- psi_control[se.id, idx.control, drop = FALSE]
  
  df.exp <- data.frame(Junction = se.id,
                       Gene = gene,
                       PSI = c(as.numeric(exp.case), as.numeric(exp.control)),
                       Source = c(rep("Tumor", length(exp.case)),
                                  rep("Normal", length(exp.control))),
                       Sample = c(colnames(exp.case), colnames(exp.control)),
                       stringsAsFactors = FALSE)
  
  ## error function, if an error exist, abort 
  an.error.occured <- FALSE
  tryCatch( {mod.mixed <- mclust::densityMclust(df.exp$PSI, G=1:3) },
            error = function(e) {an.error.occured <<- TRUE} )
  
  if(!an.error.occured){
    nClust <- mod.mixed$G
    clust.assigned <- levels(factor(mod.mixed$classification))
    nAssignClust <- length(clust.assigned)
    
    if(length(clust.assigned)<nClust){ #empty clusters
      newClass <- 1:length(clust.assigned)
      newAssigClass <- rep(NA,nrow(df.exp))
      for(c in 1:length(newClass)){
        newAssigClass[ mod.mixed$classification == clust.assigned[c] ] = c
      }
      df.exp$Classification <- factor(newAssigClass)
    }else{
      df.exp$Classification <- factor(mod.mixed$classification)
    }
    save(df.exp, file = paste0(clust_dir, "/", se.id, ".RData"))
  }
  store_df <- as.data.frame(t(data.frame(df.exp[, c("Classification")])))
  rownames(store_df) <- unique(df.exp$Junction)
  colnames(store_df) <- df.exp$Sample
  cluster_DF <- rbind(cluster_DF, store_df)
}

cluster_DF <- cluster_DF[2:nrow(cluster_DF), ]

## compute cluster composition 
computeClusterComposition <- function(clust_dir = "./clust_output", events){
  
  models <- list.files(clust_dir, pattern = ".RData")
  junction_ids <- gsub(".RData", "", models)
  
  main.counts <- mclapply(seq_along(models), function(i){
    
    jannot <- dplyr::filter(events, event_id == junction_ids[i])
    load(file = file.path(clust_dir, models[i]))
    nClust <- max(as.vector(df.exp$Classification))
    
    #Compute clusters composition in each model (i.e. junction) 
    cluster.counts <- lapply(1:nClust, function(j){
      
      res <- df.exp %>% filter(Classification==j)
      counts <- table(res$Source)
      perc <- as.numeric(counts)/sum(as.numeric(counts))
      perc <- round(perc*100,digits=1)
      
      return(data.frame(Junction = jannot$event_id, Gene = jannot$gene_symbol,
                        Source = factor(names(counts), levels = c("Tumor", "Normal")),
                        Counts = as.numeric(counts),
                        Percent = perc,
                        Cluster = paste0("C", j),
                        stringsAsFactors = FALSE)
      )
      
      
    })
    # rbind clusters for the model
    return(do.call(rbind, cluster.counts))
    
  }, mc.cores = 4)
  
  #rbind for all models
  return(do.call(rbind, main.counts))
}
clust <- list()
clust$composition <- computeClusterComposition(clust_dir, events)

## compute mean PSI
computeClusterMeanPSI <- function(clust_dir = "./clust_output", events){
  
  models <- list.files(clust_dir, pattern = ".RData")
  junction_ids <- gsub(".RData", "", models)
  
  main.counts <- mclapply(seq_along(models), function(i){
    
    jannot <- dplyr::filter(events, event_id == junction_ids[i])
    # Load df.exp data frame
    # Cluster classification for a given model (i.e junction) 
    load(file = file.path(clust_dir, models[i]))
    nClust <- max(as.vector(df.exp$Classification))
    
    #Compute clusters composition in each model (i.e. junction) 
    cluster.counts <- lapply(1:nClust, function(j){
      
      res <- df.exp %>% filter(Classification==j)
      meanPSI <- mean(res$PSI)
      sdPSI <- sd(res$PSI)
      
      return(data.frame(Junction = jannot$event_id, Gene = jannot$gene_symbol,
                        Cluster = j, meanPSI = meanPSI,sdPSI = sdPSI, 
                        stringsAsFactors = FALSE))
      
    })
    # rbind clusters for the model
    return(do.call(rbind, cluster.counts))
    
  }, mc.cores = 4)
  
  #rbind for all models
  return(do.call(rbind, main.counts))
}
clust$meanPSI <- computeClusterMeanPSI(clust_dir, events)

## select clusters 
params <- list(min_size = 5, #number of patients
               tumor_purity = 90, # percent of tumor purity
               deltaPSI = 0.1, #minimum deltaPSI
               pvalue = 0.05) #minimum pvalue  
selectClusters <- function(clust, params, clust_dir){
  res <- clust$composition %>% filter(Source=="Tumor", 
                                      Percent > params$tumor_purity, 
                                      Counts > params$min_size)
  res$Cluster <- as.numeric(gsub("C", "", res$Cluster))

  main.deltaPSI <- mclapply(seq(1, nrow(res)), function(i){
    
    # Load df.exp data frame
    # Cluster classification for a given model (i.e junction) 
    load(file = paste0(clust_dir, res$Junction[i], ".RData"))
    
    tumor.cluster.exp <- df.exp %>% dplyr::filter(Junction==res$Junction[i],
                                                  Classification==res$Cluster[i],
                                                  Source == "Tumor")
    
    normal.exp <- df.exp %>% dplyr::filter(Junction==res$Junction[i],
                                           Source == "Normal")
    
    if (nrow(tumor.cluster.exp) > 2 || nrow(normal.exp) > 2) {
      deltaPSI <- mean(tumor.cluster.exp$PSI)-mean(normal.exp$PSI)
      test <- wilcox.test(tumor.cluster.exp$PSI, normal.exp$PSI, exact = FALSE)
      cluster.test <- data.frame(Junction = res$Junction[i], 
                                 Cluster = res$Cluster[i],
                                 meanPSI_T = mean(tumor.cluster.exp$PSI),
                                 meanPSI_N = mean(normal.exp$PSI),
                                 deltaPSI_TvsN = deltaPSI,
                                 absDeltaPSI_TvsN = abs(deltaPSI),
                                 Wilcox_pval = test$p.value,
                                 stringsAsFactors = FALSE)
      
    }else{
      cluster.test <- data.frame(Junction = res$Junction[i], 
                                 Cluster = res$Cluster[i],
                                 meanPSI_T = NaN,
                                 meanPSI_N = NaN,
                                 deltaPSI_TvsN = NaN,
                                 absDeltaPSI_TvsN = NaN,
                                 Wilcox_pval = NaN,
                                 stringsAsFactors = FALSE)
    }
    return(cluster.test)
    
    
  }, mc.cores = 4)
  
  main.deltaPSI <- do.call(rbind, main.deltaPSI)
  
  res2 <- main.deltaPSI %>% filter(absDeltaPSI_TvsN>params$deltaPSI, 
                                   Wilcox_pval < params$pvalue)
  ts_clusters <- inner_join(res, res2, by = c("Junction", "Cluster") )
  
  return(ts_clusters)
  
}
clust$ts_clusters <- selectClusters(clust, params, clust_dir)
length(unique(clust$ts_clusters$Gene)) #number of unique genes 
length(unique(clust$ts_clusters$Junction)) #number of exon skipping events

## Recover patient assignment and expression for TS clusters 
gatherClusterExp <- function(clust_dir = "./clust_output", events, ts_clusters){
  
  models <- list.files(clust_dir, pattern = ".RData")
  junction_ids <- gsub(".RData", "", models)
  
  idx_keep <- which(junction_ids %in% ts_clusters$Junction)
  models <- models[idx_keep]
  junction_selected <- junction_ids[idx_keep]
  
  list.dfexp <- mclapply(seq_along(models), function(i){
    
    jannot <- dplyr::filter(events, event_id == junction_selected[i])
    # Load df.exp data frame
    # Cluster classification for a given model (i.e junction) 
    load(file = file.path(clust_dir, models[i]))
    
    return(df.exp)
    
  }, mc.cores = 4)
  
  
  #rbind for all models
  return(do.call(rbind, list.dfexp))
}
clust$cluster_exp <- gatherClusterExp(clust_dir, events, clust$ts_clusters)

## save file 
save(clust, file = "6_clustering/output/clust_lr.RData")
save.image("data/6_clustering/output/2_cluster_AS_events.RData")
















