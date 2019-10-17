#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# Szi Kay Leung
# 17/10/2019: Expression Subsampling of RNAseq detected by Isoseq (using chained output from chain_samples.py)

list_of_packages <- c("reshape2","dplyr","ggplot2","wesanderson")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(wesanderson))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))


#********************** Prepare input chained output file (From tofu cupcake chain_sample.py)
chain_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN"
chain_input <- read.table(paste0(chain_input_dir,"/all_samples.chained_count.txt"), header = TRUE)

#### QUESTION: why does the input file of SQANTI2 have NAs and 0s in places where there should be a count??
chain_input_sqanti2_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN/SQANTI2/"
Q21_chain_input <- read.table(paste0(chain_input_sqanti2_dir,"/Q21.collapsed.filtered.rep_classification.txt"), header = TRUE, fill = TRUE)

# check concordance between chain input file to SQANTI2 of FL isoforms and SQANTI2 output 
check_concordance <- function(isoform.ID){
  chain_input <- chain_input[,order(colnames(chain_input))]
  chain_input <- chain_input[,c(17,1:16)]
  print(chain_input[chain_input$superPBID == isoform.ID,])
  print(Q21_chain_input[Q21_chain_input$isoform == isoform.ID,c(1,43:58)])
}
# check_concordance("PB.1.1")
# note some isoforms do not have FL but then observed in SQANTI output
# check_concordance(PB.101.6)


# Calculate Log(TPM) for FL columns
Calculate_LogTPM <- function(column){
  total_fl <- sum(column, na.rm = TRUE)
  TPM <- (column*(10**6))/total_fl
  log10TPM <- log10(TPM)
  log10TPM
}


Log_Q21_chain_input <- data.frame(Q21_chain_input[1:42], apply(Q21_chain_input[43:ncol(Q21_chain_input)],2, Calculate_LogTPM))
colnames(Log_Q21_chain_input)[43:ncol(Q21_chain_input)] <- paste(colnames(Log_Q21_chain_input)[43:ncol(Q21_chain_input)], "LogTPM", sep = "_") 


#********************** Prepare input featurecount file (in TPM)
featurecount_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/FeatureCounts/WT8_OLD"
featurecount_input <- read.table(paste0(featurecount_input_dir,"/WT8Merge.transcript_id_tpm.txt"), header = TRUE)

# simplify column names to only extract sample name
colnames(featurecount_input)[7:dim(featurecount_input)[2]] <- word(colnames(featurecount_input)[7:dim(featurecount_input)[2]],c(11),  sep = fixed ('.'))   


#********************** sample specific: K17

sample_detection_point <- function(FL.K17_LogTPM){
  
  print(FL.K17_LogTPM)
  
  #K17_IsoSeq <- all[all$Sample == "K17",relevant_cols]
  K17_IsoSeq <- Log_Q21_chain_input %>%  select(c(1:42),FL.K17_LogTPM)
  # remove input with NA and 0 counts (log = -inf)
  K17_IsoSeq<- K17_IsoSeq[!is.na(K17_IsoSeq$FL.K17_LogTPM),]
  K17_IsoSeq<- K17_IsoSeq[K17_IsoSeq$FL.K17_LogTPM == "-Inf",]
  
  # Duplicated transcripts if not including isoform into ID
  #K17[duplicated(K17$ID) == "TRUE",] 
  #all[all$associated_transcript == "ENSMUST00000191939.1" & all$Sample == "K17",]
  #all[all$associated_transcript == "ENSMUST00000159814.1" & all$Sample == "K17",]
  
  K17_RNASeq <- featurecount_input[,c("Geneid","Length","K17")]
  # remove input with 0 counts
  K17_RNASeq <- K17_RNASeq[K17_RNASeq$K17 > 1,]
  K17_RNASeq$K17_RNASeq_LogTPM<- log10(K17_RNASeq$K17)
  
  
  #********************** Merge
  K17_Merge <- merge(K17_IsoSeq, K17_RNASeq, by.x = "associated_transcript", by.y = "Geneid", all = "TRUE" )
  #K17_IsoSeq[duplicated(K17_IsoSeq$associated_transcript) == "TRUE",] 
  #K17_RNASeq[duplicated(K17_RNASeq$Geneid) == "TRUE",] 
  
  
  #********************** Separate by expression bins 
  res <- hist(log10(K17_RNASeq$K17), breaks=10)
  res$breaks
  
  # 8 groups
  K17_Merge$Rank <- ifelse(K17_Merge$K17_RNASeq_LogTPM >= 0 & K17_Merge$K17_RNASeq_LogTPM < 1, "0-1",
                           ifelse(K17_Merge$K17_RNASeq_LogTPM >= 1 & K17_Merge$K17_RNASeq_LogTPM < 2, "1-2",
                                  ifelse(K17_Merge$K17_RNASeq_LogTPM >= 2 & K17_Merge$K17_RNASeq_LogTPM < 3, "2-3",
                                         ifelse(K17_Merge$K17_RNASeq_LogTPM >= 3 & K17_Merge$K17_RNASeq_LogTPM < 4, "3-4", "NA"))))
  
  K17_Merge$Rank <- as.factor(K17_Merge$Rank)
  K17_Merge <- K17_Merge[order(-K17_Merge$K17_RNASeq),]
  K17_Merge$Detected <- ifelse(!is.na(K17_Merge$FL.K17_LogTPM), "Yes","No")
  
  
  # plot of the counts of Isoseq detected based on the expression 
  p <- ggplot(K17_Merge, aes(x = Rank, fill = Detected)) +
    geom_bar() + 
    labs(x = "RNASeq Transcript Expression (Log10TPM)", y = "Number of Transcripts", fill = "Detected by IsoSeq") + 
    theme_bw()
  
  p
  
  # proportions of expression detected 
  a <- count(K17_Merge, vars=c("Rank","Detected"))
  a <- a %>% spread(Detected, freq)
  
  b <- count(K17_Merge, vars=c("Rank"))
  
  proportions <- merge(a,b, by = "Rank")
  proportions$perc_covered <- (proportions$Yes/proportions$freq)*100
  
  p <- ggplot(proportions, aes(x = Rank, y = perc_covered)) + 
    geom_point() + 
    labs(x = "RNASeq Transcript Expression (Log10TPM)", y = "Percentage covered in Isoseq", title = "All Isoseq") + 
    theme_bw()
  
  p
}

sample_detection_point("FL.M21_LogTPM")
