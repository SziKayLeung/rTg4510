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


Merge_Gene_Input <- function(sqanti_dir, tofu_dir, featurecount_dir,sample_name,ftr_cnt_name){
  sqanti <- read.table(print(paste0(sqanti_dir,"/",sample_name,".collapsed.filtered.rep_classification.txt")),header=T, stringsAsFactors=F, sep="\t")
  datatable(head(sqanti))
  
  tofu <- read.table(print(paste0(tofu_dir,"/",sample_name,".collapsed.filtered.abundance.txt")),header=T, stringsAsFactors=F, sep="\t")
  datatable(head(tofu))
  
  TPM <- read.table(print(paste0(featurecount_dir,"/",ftr_cnt_name,"_tpm.txt")),header=T, stringsAsFactors=F, sep="\t")
  colnames(TPM)[7] <- "TPM"
  TPM <- TPM[,c("Geneid","TPM")] #just the columns that matter
  
  # merge and format
  sqanti_counts <- merge(sqanti,tofu,by.x="isoform", by.y="pbid")
  FSM <- subset(sqanti_counts, sqanti_counts$structural_category=="full-splice_match")
  FSM_TPM <- merge(FSM,TPM,by.x="associated_transcript", by.y="Geneid")
  FSM_TPM <<_ FSM_TPM
}


SQANTI2="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/WT8_ISOSEQ/IsoSeq3.1.2/SQANTI2"
TOFU="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/WT8_ISOSEQ/IsoSeq3.1.2/TOFU"
FEATURECOUNTS="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/FeatureCounts/WT8_Transcript"

Merge_Transcript_Input(SQANTI2,TOFU,FEATURECOUNTS,"WT8IsoSeq","WT8MergeBasic.transcript_id.nodup")
