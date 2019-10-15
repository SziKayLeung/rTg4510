#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# Szi Kay Leung
# 15/10/2019: Expression Subsampling of RNAseq detected by Isoseq

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


#********************** Define variables
Tg4510_WT_2mos <- c("K17","M21","Q21")
Tg4510_WT_8mos <- c("K23","O23","S23")
Tg4510_TG_2mos <- c("O18","K18","S18")
Tg4510_TG_8mos <- c("L22","Q20","K24")


#********************** Prepare input featurecount file (in TPM)
featurecount_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/FeatureCounts/WT8_OLD"
featurecount_input <- read.table(paste0(featurecount_input_dir,"/WT8Merge.transcript_id_tpm.txt"), header = TRUE)

# simplify column names to only extract sample name
colnames(featurecount_input)[7:dim(featurecount_input)[2]] <- word(colnames(featurecount_input)[7:dim(featurecount_input)[2]],c(11),  sep = fixed ('.'))   

#********************** Prepare input sqanti classification files
sqanti2_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/SQANTI2"
filenames <- list.files(path = sqanti2_input_dir , pattern = "collapsed.filtered.rep_classification.filtered_lite_classification.txt", full.names = TRUE)
sqanti_files <- lapply(filenames, function(x) read.table(paste(x), header = TRUE))

filenames
samples <- c("B21","C20","C21","E18","K17","K18","K23","K24","L22","M21","O18","O23","Q20","Q21","S18","S23")
counter=1
for (i in 1:length(sqanti_files)){
  sqanti_files[[i]]$Sample <- print(samples[counter])
  counter = counter + 1
}

all <- bind_rows(sqanti_files)

# convert SQANTI FL to TPM (based on E.Tseng's SQANTI2.report https://github.com/Magdoll/SQANTI2)
total_fl <- sum(all$FL, na.rm=T)
all$ISOSEQ_TPM <- all$FL*(10**6)/total_fl
all$Log_ISOSEQ_TPM <- log10(all$ISOSEQ_TPM)

all$Sample <- factor(all$Sample)
all$RTS_stage <- factor(all$RTS_stage)
all$Phenotype <-ifelse(all$Sample %in% Tg4510_WT_2mos, "WT_2mos",
                       ifelse(all$Sample %in% Tg4510_WT_8mos, "WT_8mos", 
                              ifelse(all$Sample %in% Tg4510_TG_2mos, "TG_2mos", 
                                     ifelse(all$Sample %in% Tg4510_TG_8mos, "TG_8mos", "J20"))))

# distinguish novel transcripts by 
all$associated_transcript <- ifelse(all$associated_transcript %in% "novel",
                                    paste0("Novel_",all$associated_gene),
                                    all$associated_transcript)

# concatenate at to distinguish FSM_Coding+Multiexonic
all$ID <- as.factor(paste(all$isoform,all$associated_transcript,all$structural_category,all$subcategory,all$coding, 
                          sep = "//"))


#********************** K17 specific
relevant_cols <- c("Sample","isoform","chrom","associated_gene","associated_transcript","FL","ISOSEQ_TPM","Log_ISOSEQ_TPM",
                   "within_cage_peak","ID")

K17_IsoSeq <- all[all$Sample == "K17",relevant_cols]

# Duplicated transcripts if not including isoform into ID
#K17[duplicated(K17$ID) == "TRUE",] 
#all[all$associated_transcript == "ENSMUST00000191939.1" & all$Sample == "K17",]
#all[all$associated_transcript == "ENSMUST00000159814.1" & all$Sample == "K17",]

K17_RNASeq <- featurecount_input[,c("Geneid","Length","K17")]
# remove input with 0 counts
K17_RNASeq <- K17_RNASeq[K17_RNASeq$K17 > 1,]
K17_RNASeq$RNASEQ_TPM <- log10(K17_RNASeq$K17)


#********************** Merge
K17_Merge <- merge(K17_IsoSeq, K17_RNASeq, by.x = "associated_transcript", by.y = "Geneid", all = "TRUE" )

colnames(K17_Merge)[13] <- "K17_RNASeq_LogTPM"
colnames(K17_Merge)[7] <- "K17_IsoSeq_LogTPM"

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
K17_Merge$Detected <- ifelse(!is.na(K17_Merge$FL), "Yes","No")


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

