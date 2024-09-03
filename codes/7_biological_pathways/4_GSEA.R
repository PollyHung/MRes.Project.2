source("code/0_pkgs&func.R")


edgeR <- read.csv("data/5_differential_expression/DTE_edgeR.csv") %>% dplyr::rename(isoform_id=X)
gsea_list <- edgeR %>% 
  dplyr::select(isoform_id, logFC) %>% 
  dplyr::arrange(desc(logFC))
classification <- read.delim("data/2_sqanti/filter_rule_RulesFilter_result_classification.txt") %>%  
  dplyr::filter(filter_result == "Isoform") #%>% 
  #dplyr::select(isoform, associated_gene)
gsea_list <- merge(gsea_list, classification, by.x = "isoform_id", by.y = "isoform")
gsea_list <- gsea_list %>% dplyr::select(associated_gene, logFC)


#write.table(gsea_list, "data/7_gsea/differential_exon_expression.txt", 
#            sep = "\t", quote = F, row.names = F, col.names = F)


reactome_neg <- read.delim("data/7_gsea/DTE_reactome/gsea_report_for_na_neg_1721979449212.tsv", check.names = FALSE)
reactome_pos <- read.delim("data/7_gsea/DTE_reactome/gsea_report_for_na_pos_1721979449212.tsv", check.names = FALSE)
gobp_neg <- read.delim("data/7_gsea/DTE_gobp/gsea_report_for_na_neg_1721979533819.tsv", check.names = FALSE)
gobp_pos <- read.delim("data/7_gsea/DTE_gobp/gsea_report_for_na_pos_1721979533819.tsv", check.names = FALSE)

reactome <- rbind(reactome_pos, reactome_neg)
gobp <- rbind(gobp_pos, gobp_neg)

reactome$NAME <- gsub("reactome_", "", tolower(reactome$NAME)) 
reactome$NAME <- gsub("_", " ", reactome$NAME)
reactome$NES <- as.numeric(reactome$NES)
reactome$`FDR q-val` <- as.numeric(reactome$`FDR q-val`)


gsea_df <- reactome %>% mutate(Sign = ifelse(NES >= 0, "Positive", "Negative")) %>% 
  dplyr::filter(`FDR q-val` < FDR) 

## do we want to subset? if top.n is numeric (i.e., specify top x layer) then we would subset, 
if(is.numeric(top.n)){
  top_positive <- gsea_df %>% filter(Sign == "Positive") %>% top_n(top.n, NES)
  top_negative <- gsea_df %>% filter(Sign == "Negative") %>% top_n(-top.n, NES)
  top_gsea <- rbind(top_positive, top_negative)
} else { ## if its NULL, then we will use all 
  top_gsea <- gsea_df
}

## if the initial filtering included pathways with FDR > 0.05, we would tone down the color by increasing transparency 
top_gsea$alpha <- ifelse(top_gsea$`FDR q-val`< 0.05, 1, 0.75)

## make human readable names 
lower_names <- tolower(top_gsea$NAME)
process_name <- function(name) {
  parts <- strsplit(name, "_")[[1]]  # Split each name and select the list element
  if (length(parts) > 1) {
    parts <- parts[-1]  # Remove the first element
  }
  paste(parts, collapse = " ")  # Join the remaining elements with space
}
lower_names <- sapply(lower_names, process_name)
top_gsea$NAME <- lower_names 

height <- nrow(top_gsea)/7+0.5
print(paste0("suggested plot height: ", height, " inches")) 
gsea_height <- make.names("gsea_height")
assign(gsea_height, value = height, envir = .GlobalEnv)

width <- (top_gsea$NAME %>% unique %>% nchar %>% max)/20 + 1
print(paste0("suggested plot width: ", width, " inches")) 
gsea_width <- make.names("gsea_width")
assign(gsea_width, value = width, envir = .GlobalEnv)

## plot! 
p <- ggplot(top_gsea, aes(x = reorder(NAME, as.numeric(NES)), y = as.numeric(NES), fill = Sign, alpha = alpha)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("Positive" = "#FF407D", "Negative" = "#1B3C73"), 
                    labels = c("Positive" = "Up-reg", "Negative" = "Down-reg")) +
  scale_alpha_identity() + 
  labs(x = "Pathway", 
       y = "NES", 
       title = title,
       fill = "Direction") + theme_bw() +
  theme(legend.direction = "vertical", 
        legend.byrow = FALSE, 
        legend.title = element_blank(), 
        legend.position = "top", 
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7, margin(t = 0, unit = "pt")),
        axis.ticks.length = unit(1, "mm"), 
        title = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.key.size = unit(2.5, 'mm'), 
        legend.margin = margin(t = 0, unit = "pt"), 
        strip.text.x = element_text(size = 7, face = "bold"), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"))



















