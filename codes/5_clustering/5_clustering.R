#---
#title: "Clustering Isoforms by PSI values"
#author: "You"
#date: "2024-07-29"
#output: html_document
#---

## Setting up =================================    
source("code/0_pkgs&func.R")
PREPARE_MATRIX <- FALSE
CLUSTER_EVENT <- TRUE

## Prepare matrix =============================
if(PREPARE_MATRIX){
  # [1] result from SUPPA2 
  psi_events <- read.table(file="data/5_differential_expression/SUPPA/psi_local.psi", 
                           sep="\t", check.names = FALSE)
  psi_isoforms <- read.table(file="data/5_differential_expression/SUPPA/psi_isoform.psi", 
                             sep="\t", check.names = FALSE)
  
  # [2] tumour and normal sample names
  normal <- colnames(psi_events)[grepl("FN", colnames(psi_events))]
  tumour <- setdiff(colnames(psi_events), normal)
  
  # [3] dealing with missing values 
  #psi_events[is.na(psi_events)] <- 0 # 1540 events 
  #psi_isoforms[is.na(psi_isoforms)] <- 0 # 26625 isoforms 
  
  # [4] remove rows if all of the samples have expression of 0
  #psi_events <- psi_events %>% dplyr::filter(rowSums(psi_events) > 0) # 1037 events 
  #psi_isoforms <- psi_isoforms %>% dplyr::filter(rowSums(psi_isoforms) > 0) # 26625 isoforms
  
  # [5] split the samples into normal and tumour psi 
  psi_events.T <- psi_events[, tumour]
  psi_events.N <- psi_events[, normal]
  psi_isoforms.T <- psi_isoforms[, tumour]
  psi_isoforms.N <- psi_isoforms[, normal]
  
  # [6] save these 
  save(list = c("psi_events.T", "psi_events.N"), file = "data/5_differential_expression/SUPPA/psi_events.rds")
  save(list = c("psi_isoforms.T", "psi_isoforms.N"), file = "data/5_differential_expression/SUPPA/psi_isoforms.rds")
}

## clustering event PSI =======================
if(CLUSTER_EVENT){
  load("data/5_differential_expression/SUPPA/psi_events.rds")
  clust_dir = "data/6_clustering/AS_events"
  
  ## create result directory: 
  if(dir.exists(clust_dir)){unlink(clust_dir, force = T, recursive = T)}
  dir.create(clust_dir, recursive = TRUE)
  
  ## flag events to remove 
  events <- read.table(file = "data/5_differential_expression/SUPPA/sqanti_filter.events.ioe", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  flags <- data.frame(Junction.id = events$event_id,  # column 1 junction ID 
                      Remove = rep(FALSE, nrow(events))) # column 2: T/F to remove the junction? 
  for(i in 1:nrow(events)){
    # [1] get the AS event ID and pacbio ID 
    se.id <- events$event_id[i]
    gene <- events$gene_id[i]
    # [2] check which column is not NA? 
    idx.case <- which(!is.na(psi_events.T[se.id, ]))
    idx.control <- which(!is.na(psi_events.N[se.id, ]))
    # [3] subset expression matrix that is not NA 
    exp.case <- psi_events.T[se.id, idx.case, drop = FALSE]
    exp.control <- psi_events.N[se.id, idx.control, drop = FALSE]
    # [4] check if the subsetted expression matrix has length over 0
    if(length(exp.case) < 10 || length(exp.control) < 3){
      flags$Remove[i] <- TRUE # if not, remove it!
      next}
    if(sd(exp.case)==0 || sd(exp.control)==0){
      flags$Remove[i] <- TRUE
      next}
  }
  
  ## filter out those flagged removal = true 
  id.keep <- flags$Junction.id[!flags$Remove]
  events_filt <- events %>% filter(event_id %in% id.keep)
  cat("INFO: Events to be clustered: ", nrow(events_filt), "\n")
  
  ## Check if the rest of them are normally distributed? 
  distribution <- data.frame(Junction.id = events_filt$event_id,  # column 1 junction ID 
                             p_val.N = rep(FALSE, nrow(events_filt)), 
                             p_val.T = rep(FALSE, nrow(events_filt))) # column 2: T/F to remove the junction? 
  for(i in 1:nrow(events_filt)){
    # [1] get event id 
    se.id <- events_filt$event_id[i]
    gene <- events_filt$gene_id[i]
    # [2] get ID? 
    idx.case <- which(!is.na(psi_events.T[se.id, ]))
    idx.control <- which(!is.na(psi_events.N[se.id, ]))
    # [3] get expression matrix 
    exp.case <- psi_events.T[se.id, idx.case, drop = FALSE]
    exp.control <- psi_events.N[se.id, idx.control, drop = FALSE]
    # [4] perform test 
    distribution$p_val.N[i] <- NA
    distribution$p_val.T[i] <- NA
    tryCatch({
      test_T <- shapiro.test(as.numeric(exp.case))
      distribution$p_val.T[i] <- test_T$p.value
    }, error = function(e) {
      distribution$p_val.T[i] <- NA
    })
    tryCatch({
      test_N <- shapiro.test(as.numeric(exp.control))
      distribution$p_val.N[i] <- test_N$p.value
    }, error = function(e) {
      distribution$p_val.N[i] <- NA
    })
    # [5] FDR? 
    #distribution$p_adj.N <- p.adjust(distribution$p_val.N)
    #distribution$p_adj.T <- p.adjust(distribution$p_val.T) 
    # [6] normal vs not normal?
    table(ifelse(distribution$p_val.N < 0.05, "not normal", "normal"))
    table(ifelse(distribution$p_val.T < 0.05, "not normal", "normal")) 
  }
  
  ## register for parallel computing 
  nmodels <- nrow(events_filt)
  registerDoParallel(cores = 8)
  
  ## GMM modelling 
  foreach(i = 1:nmodels, .combine = rbind) %dopar% {
    # [1] get event id 
    se.id <- events_filt$event_id[i]
    gene <- events_filt$gene_id[i]
    # [2] get ID? 
    idx.case <- which(!is.na(psi_events.T[se.id, ]))
    idx.control <- which(!is.na(psi_events.N[se.id, ]))
    # [3] get expression matrix 
    exp.case <- psi_events.T[se.id, idx.case, drop = FALSE]
    exp.control <- psi_events.T[se.id, idx.control, drop = FALSE]
    # [4] histogram? 
    p <- ggplot(data.frame(value = c(as.numeric(exp.case), as.numeric(exp.control)),
                           group = c(rep("Case", length(exp.case)), rep("Control", length(exp.control)))), 
                aes(x = value, fill = group)) +
      geom_density(alpha = 0.5)
    ggsave(paste0("data/6_clustering/plot/", se.id, ".png"), p, width = 3, height = 3, units = "in", dpi = 600)
    # [5] normal distribution? 
    shapiro.test(as.numeric(exp.case))
    shapiro.test(as.numeric(exp.control))
    
    
    # [5] forming a dataframe 
    df.exp <- data.frame(Junction = se.id, 
                         Gene = gene, 
                         PSI = c(as.numeric(exp.case), as.numeric(exp.control)), 
                         Source = c(rep("Tumour", length(exp.case)), 
                                    rep("Normal", length(exp.control))), 
                         Sample = c(colnames(exp.case), colnames(exp.control)), 
                         stringsAsFactors = FALSE)
    # [6] first fitting, this is to prevent error stopping the loop 
    an.error.occured <- FALSE
    tryCatch( {mod.mixed <- mclust::densityMclust(df.exp$PSI, G=1:3) },
              error = function(e) {an.error.occured <<- TRUE} ) 
    # [7] if an error did NOT occur: 
    if(FALSE){ #!an.error.occured
      # [7.1] get parameters from first fit 
      nClust <- mod.mixed$G
      clust.assigned <- levels(factor(mod.mixed$classification))
      nAssignClust <- length(clust.assigned)
      # [7.2] if there are empty clusters, or else... 
      if(length(clust.assigned)<nClust){ 
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
  }
}



















