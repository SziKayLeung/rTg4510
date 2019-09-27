library(stringr) # Extract strings
library(dplyr)   # dataframe splitting
library("ggpubr") # correlation plot
library(plotly) # interactive graphs
library(ggplot2)
library(grid) # ggplot add corr value
library(gridExtra) # ggplot add corr value
library(ggthemes) # ggplot minimum theme
library(DT)
library(tidyr) #gather 

# Input Sqanti classification file, tofu's abundance filtered file, and TPM calculated using script processfeaturecounts.R
Merge_Transcript_Input <- function(sqanti_dir, tofu_dir, featurecount_dir,sample_name,ftr_cnt_name){
  sqanti <- read.table(print(paste0(sqanti_dir,"/",sample_name,".collapsed.filtered.rep_classification.txt")),header=T, stringsAsFactors=F, sep="\t")
  datatable(head(sqanti))

  tofu <- read.table(print(paste0(tofu_dir,"/",sample_name,".collapsed.filtered.abundance.txt")),header=T, stringsAsFactors=F, sep="\t")
  datatable(head(tofu))

  TPM <- read.table(print(paste0(featurecount_dir,"/",ftr_cnt_name,"_tpm.txt")),header=T, stringsAsFactors=F, sep="\t")
  colnames(TPM)[7] <- "TPM"
  TPM <- TPM[,c("Geneid","TPM")] #just the columns that matter

  #merge and format
  sqanti_counts <- merge(sqanti,tofu,by.x="isoform", by.y="pbid")
  FSM <- subset(sqanti_counts, sqanti_counts$structural_category=="full-splice_match")
  FSM_TPM <- merge(FSM,TPM,by.x="associated_transcript", by.y="Geneid")
  FSM_TPM <<- FSM_TPM
}

# calculate correlation
corr.value <- cor(FSM_TPM$count_fl,FSM_TPM$TPM)

# plot correlation in density plot 
density_plot <- function(dat,x.var,y.var, x_lab, y_lab){
  
  corr <- grobTree(textGrob(paste("R : ", round(corr.value, 4)), x = 0.05, y = 0.97, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 11, fontface = "italic")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1) +
    theme_tufte() +
    labs(x = print(paste(x_lab)), y = print(paste(y_lab))) + 
    geom_smooth(method=lm, colour = "black")
  print(p)
}