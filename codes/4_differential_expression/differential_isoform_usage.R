## This script performs differential isoform usage analysis 
## Aim: identify changes in the relative usage of different transcripts (isoforms) 
##      of the same gene between normal and tumour
## Methods: 
##        DRIMSeq –– requiring raw counts 
##        IsoformSwitchAnalyzeR -- for isoform switching analysis 



## expression matrix: raw count matrix from long read sequences using flair quantify 
## experiment design: 3 column df with sample name, group, and a new sample name 
## annotation features: a gff3 file recording the annotations for isoforms mentioned in 
##                      expression matrix 
## classification file: a classification file from SQANTI3 post to filtering 
## gtf: the gtf file from SQANTI3 post filtering 
library(DRIMSeq)
library(dplyr)
library(tidyverse)
library(magrittr)
library(IsoformSwitchAnalyzeR)

set.seed(55)



## Stage 1: Differential isoform usage -----------------------------------------
# loading data...
# sample_list <- experiment_design$sample_id # we will use the modified data name 
# expression_matrix <- as.data.frame(expression_matrix)
# colnames(expression_matrix) <- experiment_design$sample_id # rename expression matrix accordingly 
# expression_matrix$feature_id <- rownames(expression_matrix) # add feature_id column
# expression_matrix$gene_id <- gsub("\\.\\d+$", "", rownames(expression_matrix)) # add gene id column
# expression_matrix <- expression_matrix[, c("gene_id", "feature_id", sample_list)] # reorder dataframe columns 
load("LongReads_TPM.Rd")
load("cor_mtx.Rd")
load("classification.Rd")


condition_mtx <- condition_mtx %>% dplyr::rename(sample_id = sample) %>% 
  dplyr::select(sample_id, group)
long_reads <- long_reads[classification$isoform, ]
long_reads$feature_id <- rownames(long_reads)
long_reads$gene_id <- gsub("\\.\\d+$", "", rownames(long_reads))
long_reads <- long_reads[, c("gene_id", "feature_id", condition_mtx$sample_id)] # reorder dataframe columns 


# build data object 
dtu_object <- dmDSdata(counts = long_reads, samples = condition_mtx)
png("plots/DTU/01_data_summary.png", width = 4, height = 3, units = "in", res = 600)
plotData(dtu_object) ## data summary plot 
dev.off()

# Filtering of lowly expressed transcripts for 
# 1. minimal expression using min_samps_feature_expr and min_feature_expr parameters
# 2. minimal proportion with min_samps_feature_prop and min_feature_prop
dtu_filter <- dmFilter(x = dtu_object, 
                       min_samps_feature_expr = 6, ## situation where transcript is expressed in one condition but not another
                       min_gene_expr = 3, ## at least 3 estimated counts 
                       min_feature_expr = 3) 
png("plots/DTU/02_data_summary_post_filter.png", width = 4, height = 3, units = "in", res = 600)
plotData(dtu_filter) ## data summary plot 
dev.off()

# create design matrix 
design_full <- model.matrix(~ group, data = samples(dtu_filter))
colnames(design_full) <- c("Intercept", "groupC2") # rename design matrix column 

# calculate precision
dtu_filter <- dmPrecision(dtu_filter, design = design_full)
png("plots/DTU/03_precision_estimates.png", width = 5, height = 3, units = "in", res = 600)
plotPrecision(dtu_filter)
dev.off()

# proportion estimates 
dtu_filter <- dmFit(dtu_filter, design = design_full, verbose = 1)

# Testing for differential transcript usage
dtu_filter <- dmTest(dtu_filter, coef = "groupC2", verbose = 1, bb_model = TRUE)
png("plots/DTU/04_p_values.png", width = 4, height = 3, units = "in", res = 600)
plotPValues(dtu_filter)
dev.off()
png("plots/DTU/04_p_values_features.png", width = 4, height = 3, units = "in", res = 600)
plotPValues(dtu_filter, level = "feature")
dev.off()

# results
res <- results(dtu_filter)
res_sig <- res %>% dplyr::filter(adj_pvalue < 0.05) %>% dplyr::arrange(adj_pvalue)
#res_sig <- read.csv("docs/results/DTU_DRIMSeq.csv", row.names = 1)

# individual plot 
gene_id <- res_sig$gene_id
pdf("plots/DTU/05_sig_gene_transcript_change_level.pdf", width = 6.3, height = 9.7)
for(i in seq_along(gene_id)){
  top_gene_id <- gene_id[i]
  p <- plotProportions(dtu_filter, gene_id = top_gene_id, group_variable = "group", plot_type = "ribbonplot")
  print(p)
}
dev.off()

# summary plot 
# [1] volcano plot 
res$log10_adjP <- -log10(res$adj_pvalue)
volcano_plot <- ggplot(res, aes(x = lr, y = log10_adjP)) +
  geom_point(aes(color = adj_pvalue < 0.05), size = 1) +
  labs(x = "Effect Size",
       y = "-log10(adj.p)") +
  theme_bw() +
  theme(legend.position = "right") +
  geom_text_repel(data = subset(res, adj_pvalue < 0.05 & lr > 70), 
                  aes(label = gene_id), size = 3) + 
  guides(color = guide_legend(title = NULL))
ggsave("plots/DTU/volcano_plot.png", volcano_plot, width = 3, height = 2, units = "in", dpi = 600)
# [2] pie chart 
res <- res %>% mutate(significant = adj_pvalue < 0.05)
pie(x = table(res$significant))

# save files 
rm(list = normalize_files)
write.csv(res_sig, "docs/results/DTU_DRIMSeq.csv")
save.image("docs/results/DTU_DRIMSeq.RData")




## Stage 2: Isoform switching analysis =========================================

## [1] Creating a switchAnalyzeRlist --------------
load("surface_proteins.Rd")
isoformCountMatrix <- merge(long_reads, classification[, c("isoform", "associated_gene")], 
                            by.x = "feature_id", by.y = "isoform")
isoformCountMatrix <- isoformCountMatrix %>% 
  dplyr::filter(associated_gene %in% surface_proteins$entrez_symbol) %>% 
  dplyr::rename(isoform_id = feature_id)
isoformCountMatrix <- isoformCountMatrix[, c("isoform_id", condition_mtx$sample_id)]

# process the design matrix
designMatrix <- condition_mtx %>% 
  dplyr::select(sample_id, group) %>% 
  dplyr::rename(sampleID = sample_id, condition = group)

# process the gtf and fasta in terminal 
write.table(rownames(isoformCountMatrix), "~/Desktop/isoform_list.txt",
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("grep -Fwf '$seqid' '$qc_corrected.gtf.cds.gff' > '$filtered_gff'") # filter gff
system("gffread -E '$filtered_gff' -T -o '$filtered_gtf'") # gff to gtf
system("seqkit grep -f '$seqid' '$qc_corrected.fasta' > '$filtered_fasta'") # filter fasta

# creating a switch list 
SwitchList <- importRdata(isoformCountMatrix = isoformCountMatrix, 
                           designMatrix = designMatrix, 
                           isoformExonAnnoation = "~/Desktop/project/docs/data/2_sqanti/filtered_output.gtf", 
                           isoformNtFasta = "~/Desktop/project/docs/data/2_sqanti/filtered_output.fasta", 
                           showProgress = TRUE) 
summary(SwitchList)

## [2] Identifying Isoform Switches -----------------
# pre filtering (not strict, basic)
filterSwitchList <- preFilter(switchAnalyzeRlist = SwitchList,
                        geneExpressionCutoff = 1,
                        isoformExpressionCutoff = 0,
                        removeSingleIsoformGenes = TRUE)


# perform isoform switching analysis with DEXSeq
filterSwitchList <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = filterSwitchList, 
                                            alpha = 0.05, # 95% confident?
                                            dIFcutoff = 0.10, # at least 10% change 
                                            reduceToSwitchingGenes = FALSE) # run on isoform level 
extractSwitchSummary(filterSwitchList)

## [3] Analyzing Open Reading Frames -----------------
# known isoforms are automatically annotated providing fasta and gtf with cds 
# during inital data object construction stage 
"orfAnalysis" %in% names(filterSwitchList) # returns TRUE

# analyse novel isoforms 
# which of those isoforms with isoform switching issues are novel? 
novel_isoforms <- classifications %>% 
  dplyr::filter(isoform %in% filterSwitchList$isoformSwitchAnalysis$isoform_id) %>% # select 
  dplyr::filter(associated_transcript == "novel")

analyzeNovelIsoformORF(switchAnalyzeRlist = filterSwitchList, 
                       analysisAllIsoformsWithoutORF = TRUE, 
                       orfMethod = "longest.AnnotatedWhenPossible")

## [4] Extracting Nucleotide and Amino Acid Sequences ----------------
filterSwitchListAnalyzed <- extractSequence(filterSwitchList, 
                                            pathToOutput = '~/Desktop/',
                                            writeToFile=TRUE) # to avoid output when running this example data
# outputs isoformSwitchAnalyzeR_isoform_nt.fasta and isoformSwitchAnalyzeR_isoform_AA.fasta


## saving files 
rm(list = normalize_files)
write.csv(filterSwitchListAnalyzed$isoformSwitchAnalysis, "docs/results/DTU_isoformSwitch.csv")
save.image("docs/results/DTU_isoformSwitch.RData")


## plots 
plot_df <- filterSwitchList$isoformFeatures
add_gene_name <- classification %>% dplyr::select(isoform, associated_gene) %>% 
  dplyr::rename(isoform_id = isoform)
plot_df <- left_join(plot_df, add_gene_name, by= "isoform_id")


isoform_switch <- ggplot(data=plot_df, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
             size=1) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  scale_color_manual('Signficant\nIsoform Switch', values = c('grey','skyblue')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw() + 
  geom_text_repel(data = subset(plot_df, isoform_switch_q_value < 0.05 & 
                                  abs(dIF) > 0.1 &
                                  abs(gene_log2_fold_change) > 1), 
                  aes(label = isoform_id), size = 2, max.overlaps = 10)

ggsave("plots/DTU/isoform_switch_volcano.png", isoform_switch, units = "in", 
      width = 5, height = 4, dpi = 600)


genes <- subset(plot_df, isoform_switch_q_value < 0.05 & 
                  abs(dIF) > 0.1 &
                  abs(gene_log2_fold_change) > 1)








classification <- read.delim("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt") %>% 
  dplyr::filter(filter_result == "Isoform")
classification$gene_id <-  sub("\\.[^.]*$", "", classification$isoform)
differential_isoform <- read.csv("data/5_differential_expression/DTU_DRIMSeq.csv")
differential_isoform <- merge(differential_isoform, classification, by="gene_id")

load("data/5_differential_expression/DTU_DRIMSeq.RData")
differential_isoform <- merge(res_sig, classification, by="gene_id")


TMHMM <- read.table("data/4_transcoder/tmhmm_regions.txt") %>% 
  dplyr::rename(isoform=V1, location=V3) %>% 
  dplyr::select(isoform, location)
TMHMM$isoform_id <- sub("\\.[^.]*$", "", TMHMM$isoform)

differential_isoform <- merge(differential_isoform, TMHMM, by.x="isoform", by.y="isoform_id")
differential_isoform <- differential_isoform %>% 
  dplyr::filter(location %in% c("outside", "TMhelix"))






