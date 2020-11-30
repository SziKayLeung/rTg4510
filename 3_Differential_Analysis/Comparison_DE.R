# Szi Kay Leung 
# 23/11/2020: Compare results from Iso_Seq_Aaron_DTU.R and DESeq2_IsoseqFL_TranscriptExpression.R 

library("dplyr")
suppressMessages(library(VennDiagram))
suppressMessages(library(cowplot))
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq3_Tg4510/Rmarkdown_Input.R")

input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN_OLD/SQANTI3"
class_file <- read.table(paste0(input_dir,"/all_samples.chained.rep_classification.filtered_lite_classification.txt"), as.is = T, header = T, sep = "\t")

Aaron_DTU <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq3_Tg4510/Iso_Seq_Aaron_DTU_All.csv") %>% mutate("transcript_name_id" = paste0(associated_transcript,"_", isoform))
DESeq2 <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq3_Tg4510/DESeq2_IsoseqFL_TranscriptExpression.csv")

# variance 
Aaron_DTU <- Aaron_DTU %>% mutate(WT_sd_exp = apply(.[,grepl( "WT",colnames(.))], 1, function(x) sd(x))) %>%
  mutate(TG_sd_exp = apply(.[,grepl( "TG",colnames(.))], 1, function(x) sd(x))) 

# significant expression 
Aaron_DTU_sg <- Aaron_DTU %>% filter(u < 0.05) %>% mutate(Detection = "Aaron_DTU") %>% mutate("transcript_name_id" = paste0(associated_transcript,"_", isoform))
DESeq2_sg <- DESeq2 %>% filter(padj < 0.05) %>% mutate(Detection = "DESeq2") 

# ggplot 
venn_diagram_plot_twocircles <- function(set1, set2, label_set1, label_set2){
  
  p <- venn.diagram(
    x = list(set1, set2),
    category.names = c(label_set1,label_set2),
    filename = NULL,
    output=TRUE,
    
    # Circles
    lwd = 0.2,
    lty = 'blank',
    fill = c("orange","blue"),
    
    # Numbers
    cex = 3,
    fontface = "bold",
    fontfamily = "ArialMT",
    
    # Set names
    cat.cex = 2,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "ArialMT",
    #rotation = 1, 
    main = "\n\n\n\n",
    
    print.mode = "raw"
  )
  
  return(p)
  
}

p1 <- venn_diagram_plot_twocircles(Aaron_DTU_sg$transcript_name_id,DESeq2_sg$X,"NA","NA")
pdf(file="venn.pdf")
grid.draw(venn_diagram_plot_twocircles(Aaron_DTU_sg$transcript_name_id,DESeq2_sg$X,"DTU","DESeq2"))
dev.off()

# Transcripts differentially expressed uniquely in DTU, uniquely in DESeq2, common in DTU and DESeq2
DESeq2_unique <- setdiff(DESeq2_sg$X, Aaron_DTU_sg$transcript_name_id)
DTU_unique <- setdiff(Aaron_DTU_sg$transcript_name_id,DESeq2_sg$X)
common <- intersect(Aaron_DTU_sg$transcript_name_id,DESeq2_sg$X)

# Median variance 
plot_parameter_distribution <- function(list, plot_title, variable){
  if(variable == "median"){
    df <- Aaron_DTU[Aaron_DTU$transcript_name_id %in% list,c("WT_median_exp","TG_median_exp", "transcript_name_id")]
    xvar <- "Median - Iso-Seq Transcript Expression (TPM)"
    yvar <- "Median - Expression (TPM)"
  } else if(variable == "SD"){
    df <- Aaron_DTU[Aaron_DTU$transcript_name_id %in% list,c("WT_sd_exp","TG_sd_exp", "transcript_name_id")]
    xvar <- "Standard deviation - Iso-Seq Transcript Expression (TPM)"
    yvar <- "Sd - Expression (TPM)"
    
  } else {
    print("3rd argument required: median or SD")
  }
  
  p1 <- df %>% reshape2::melt(.) %>%
    ggplot(., aes(x = value, fill = variable)) + 
    geom_density(alpha = 0.2) + 
    labs(x = xvar, title = plot_title) + 
    theme_bw() + mytheme + 
    scale_fill_manual(name = "Phenotype", labels = c("WT","TG"), values = c(label_colour("WT"),label_colour("TG")))
  
  p2 <- df %>% reshape2::melt(.) %>%
    ggplot(., aes(x = transcript_name_id, y = value, group = variable, colour = variable)) + geom_line() + 
    labs(x = "Transcripts", y = yvar, title = plot_title) +
    theme_bw() + mytheme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    scale_colour_manual(name = "Phenotype", labels = c("WT","TG"), values = c(label_colour("WT"),label_colour("TG"))) + ylim(0,5000)
  
  return(list(p1,p2))
}

### grid_layout plots 
p1 <- plot_parameter_distribution(DTU_unique,"\n", "median") # DTU - Median: Transcripts Differentially Expressed
p2 <- plot_parameter_distribution(DTU_unique,"\n", "SD") # DTU - SD: Transcripts Differentially Expressed
p3 <- plot_parameter_distribution(DESeq2_unique,"\n","median") # DESeq2 - Median: Transcripts Differentially Expressed
p4 <- plot_parameter_distribution(DESeq2_unique,"\n","SD") # DESeq2 - SD: Transcripts Differentially Expressed
p5 <- plot_parameter_distribution(common,"\n", "median") #Common - Transcripts Differentially Expressed
p6 <- plot_parameter_distribution(common,"\n", "SD") #Common - Transcripts Differentially Expressed

pdf("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Comparison_Distribution.pdf", width = 11, height = 8.5)
plot_grid(p1[[1]],p2[[1]],p3[[1]],p4[[1]],p5[[1]],p6[[1]], labels = "auto", label_size = 30, label_fontfamily = "ArialMT",ncol = 2)
plot_grid(p1[[2]],p2[[2]],p3[[2]],p4[[2]],p5[[2]],p6[[2]], labels = "auto", label_size = 30, label_fontfamily = "ArialMT",ncol = 2)
dev.off()

# correlation of p-values and expression 
Aaron_DTU %>% mutate(WT_TG_ratio_median_exp = abs(WT_median_exp-TG_median_exp)) %>% .[,c("u","WT_TG_ratio_median_exp")] %>% ggplot(., aes(x = WT_TG_ratio_median_exp, y = u)) + geom_point()
Aaron_DTU %>% filter(transcript_name_id %in% DESeq2_unique) %>% mutate(WT_TG_ratio_median_exp = abs(WT_median_exp-TG_median_exp)) %>% 
  .[,c("u","WT_TG_ratio_median_exp")] %>% ggplot(., aes(x = WT_TG_ratio_median_exp, y = u)) + geom_point()

# significant values only 
DESeq2_sg_corr <- Aaron_DTU %>% mutate(WT_TG_ratio_log10mean_exp = log10(abs(WT_mean_exp-TG_mean_exp))) %>% 
  .[,c("transcript_name_id","WT_TG_ratio_log10mean_exp")] %>% merge(., DESeq2_sg, by.x = "transcript_name_id",by.y = "X", all.y = TRUE) 
# all values
DESeq2_sg_corr_all <- Aaron_DTU %>% mutate(WT_TG_ratio_log10mean_exp = log10(abs(WT_mean_exp-TG_mean_exp))) %>% 
  .[,c("transcript_name_id","WT_TG_ratio_log10mean_exp")] %>% merge(., DESeq2, by.x = "transcript_name_id",by.y = "X", all.y = TRUE) 

ggplot(DESeq2_sg_corr , aes(x = padj, y = WT_TG_ratio_log10mean_exp)) + geom_point() 
ggplot(DESeq2_sg_corr_all , aes(x = padj, y = WT_TG_ratio_log10mean_exp)) + geom_point() 

# 
cor.test(DESeq2_sg_corr_all$padj, DESeq2_sg_corr_all$WT_TG_ratio_log10mean_exp)
cor.test(DESeq2_sg_corr$padj, DESeq2_sg_corr$WT_TG_ratio_log10mean_exp)


