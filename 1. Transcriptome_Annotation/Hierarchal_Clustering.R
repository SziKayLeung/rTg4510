#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# Szi Kay Leung
# 24/09/2019: Hierachal Clustering on WT samples 

list_of_packages <- c("pheatmap", "RColorBrewer","tidyr", "reshape2")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

suppressMessages(library(pheatmap))# dataframe splitting
suppressMessages(library(RColorBrewer)) # correlation plot
suppressMessages(library(tidyr)) # correlation plot
suppressMessages(library(reshape2))

#********************** Prepare SQANTI2 classification files
sqanti2_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/SQANTI2"
filenames <- list.files(path = sqanti2_input_dir , pattern = "collapsed.filtered.rep_classification.filtered_lite_classification.txt", full.names = TRUE)
sqanti_files <- lapply(filenames, read.table)

simplify <- function(dat){
  dat <- dat[-1,c(1,7,6,22)]
  #print(dim(df))[1]
}

# keep only FSM 
FSM <- function(dat){
  dat <- dat[dat$structural_category == "full-splice_match",]
}

FL_FSM <- function(input_files){
  sqanti_files_FSM <- lapply(input_files,simplify)
  # set column names
  colnames <- c("isoform","associated_gene","structural_category","FL")
  sqanti_files_FSM <- lapply(sqanti_files_FSM, setNames, colnames)
  sqanti_files_FSM <- lapply(sqanti_files_FSM, FSM)
  sqanti_files_FSM <<- sqanti_files_FSM
}

FL_FSM(sqanti_files)

counter=1
samples <- c("B21","C21","K17","K23","M21","O23","Q21", "S23")
for (i in 1:length(sqanti_files_FSM)){
  sqanti_files_FSM[[i]]$Sample <- print(samples[counter])
  counter = counter + 1
}

all <- bind_rows(sqanti_files_FSM)
all$FL <- as.numeric(all$FL)

# sum of transcripts across genes
all <- dcast(all, associated_gene ~ Sample , value.var = "FL", fun.aggregate = sum)
m <- all[,-1]
res <- cor(m)
#cor_t <- 1 - abs(cor(t(m)))
dim(res)
round(res, 2)
res_clus <- hclust(dist(res))
plot(res_clus)


col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)
heatmap(x = cor_t, col = col, symm = TRUE)


