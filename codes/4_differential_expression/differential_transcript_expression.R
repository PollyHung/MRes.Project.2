## This script performs differential transcript expression analysis 
## Aim: analysis on differential expression of each isoforms between normal and tumor samples. 
## Methods: 
##        edgeR –– requiring raw counts 
##        NOISeq –– requiring tmm normalized matrix



source("code/pkgs&func.R")
## expression matrix: raw count matrix from long read sequences using flair quantify 
## experiment design: 3 column df with sample name, group, and a new sample name 
## annotation features: a gff3 file recording the annotations for isoforms mentioned in expression matrix 
## classification file: a classification file from SQANTI3 post to filtering 
## gtf: the gtf file from SQANTI3 post filtering 

set.seed(55)


## method 1 - NOISeq -----------------------------------------------------------
expression_matrix <- read.table("data/3_flair/pacbio2/PB2.counts.tsv", 
                                header = T, row.names = 1, check.names = F)
classifications <- read.table("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt", header = T) %>%
  dplyr::filter(filter_result == "Isoform") 
rownames(classifications) <- classifications$isoform
classifications <- classifications[rownames(expression_matrix), ]
experiment_design <- data.frame(sample = colnames(expression_matrix), 
                                group = ifelse(grepl("FN", colnames(expression_matrix)), "normal", "tumour"))

# get small stuffs 
mychroms <- classifications %>% dplyr::select(chrom, CDS_genomic_start, CDS_genomic_end) 
rownames(mychroms) <- classifications$isoform
mybiotypes <- classifications$subcategory
names(mybiotypes) <- classifications$isoform
mylength <- classifications$length
names(mylength) <- classifications$isoform

# check if the isoforms are true isoforms 
if(all(rownames(expression_matrix) %in% classifications$isoform)){
  print("all isoforms in expression matrix are true Isoforms by SQANTI3 filtering")
} else {
  match_id <- match(rownames(expression_matrix), classifications$isoform)
  expression_matrix <- expression_matrix[match_id, ]
}

# building data object 
data_object = readData(data = expression_matrix, 
                       factors = experiment_design, 
                       length = mylength, 
                       chromosome = mychroms,  
                       biotype = mybiotypes)

# normalization 
normalized_expression = NOISeq::tmm(datos = assayData(data_object)$exprs, 
                                    long = length(as.vector(as.matrix(data_object@featureData@data))), 
                                    refColumn = 1, 
                                    logratioTrim = 0.3, 
                                    sumTrim = 0.05, 
                                    k = 0, 
                                    lc = 0)
normalized_expression <- round(normalized_expression, digits = 4) ## round 

# perform differential analysis with NOISeq
myfactors <- (experiment_design[, 2]) %>% 
  as.data.frame() %>% 
  dplyr::rename(group = ".")
rownames(myfactors) <- experiment_design[, 1]
myfactors[,1]=as.factor(myfactors[,1]) 
factor=colnames(myfactors)[1] 
de_noiseq = noiseq(input = data_object, 
                   norm = "n", 
                   k = NULL, 
                   lc = 0 , 
                   factor = factor, 
                   replicates = "no")
rm(list = normalize_files)
save.image("docs/results/DTE_NOISeq.RData")
write.csv(de_noiseq@results[[1]], "docs/results/DTE_NOISeq.csv")


## method 2 - edgeR ------------------------------------------------------------

# normalize data 
tmm_factors = calcNormFactors(round(expression_matrix))

# Create DGEList Object
myedgeR = DGEList(counts = round(expression_matrix),
                  group = as.factor(experiment_design[, 2]),
                  norm.factors = tmm_factors, 
                  remove.zeros = FALSE)  

# Estimate Dispersion
myedgeR = estimateDisp(myedgeR)

# Perform Exact Test for Differential Expression:
res = exactTest(object = myedgeR, 
                pair = sort(levels(as.factor(experiment_design[, 2])), 
                          decreasing = T))
res.sort = topTags(res, n = nrow(expression_matrix))[[1]]

# save files 
rm(list = normalize_files)
save.image("docs/results/DTE_edgeR.RData")
write.csv(x = res.sort, "docs/results/DTE_edgeR.csv")


## plots -----------------------------------------------------------------------
# edgeR volcano plots 
plot <- read.csv("docs/results/DTE_edgeR.csv", row.names = 1)
plot$significance <- with(plot, ifelse(FDR < 0.05, "Significant", "Not Significant"))
plot$label <- ifelse(abs(plot$logFC) > 5 & plot$FDR < 0.05, rownames(plot), "")
edgeR_volcano <- ggplot(plot, aes(x = logFC, y = -log10(PValue), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 10) +
  labs(x = "Log Fold Change", y = "-log10(P-Value)") +
  theme(legend.title = element_blank(), 
        legend.position = "bottom") +
  geom_text_repel(aes(label = label), 
                  size = 3, 
                  box.padding = 0.1, 
                  point.padding = 0.3,
                  segment.color = 'grey50',
                  max.overlaps = 10)
ggsave("plots/DE/01_edgeR_volcano.png", 
       plot = edgeR_volcano, width = 6.3, height = 9.7, units = "in", dpi = 600)


