# Szi Kay Leung
# Aim: To run differential transcript usage using IsoformSwitchAnalyze.R
# 23/11/2020: Using output from SQANTI3 (chained data as input) with FL read count for differential transcript usage
# Prequisite: Run sqanti_gtf_modify.py --> all_samples.chained.rep_classification.filtered_lite_mod.gtf

library(dplyr)
library(IsoformSwitchAnalyzeR)

options(scipen=999) # removal of scientific notation 

# read in SQANTI2 classification file of all merged data
input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN_OLD/SQANTI3/"
class <- read.table(paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite_classification.txt"),header=T)

# PacBio isoform ID, keep only columns with isoform FL counts (and convert FL counts to TPM)
counts <- class[, c("isoform","FL.K17","FL.K18","FL.K23","FL.K24","FL.L22","FL.M21","FL.O18","FL.O23","FL.Q20","FL.Q21","FL.S18","FL.S23")] 
abundance <- cbind(counts[,1], as.data.frame(apply(counts[,-1],2, function(x) round(x/sum(x)*1000000,0))))
colnames(counts)[1] <- c("isoform_id")
colnames(abundance)[1] <- c("isoform_id")


# design dataframe input for IsoformSwitchAnalyzeR
myDesign <- data.frame(
  sampleID = colnames(counts)[-1],
  condition = c("WT", "TG", "WT", "TG","TG","WT","TG","WT","TG","WT","TG","WT")
)

# using SQANTI generated gtf and fasta, rather than gencode
#isoformExonAnnoation = paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite.gtf"),
aSwitchList <- importRdata(
  isoformCountMatrix   = counts,
  isoformRepExpression = abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite_mod.gtf"),
  isoformNtFasta       = paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite.fasta"),
)
head(aSwitchList,2)

# isoformExpressionCutoff =, filtering on isoform expression allows removal of non-used isoforms that only appear in the switchAnalyzeRlist because they were in the isoform/gene annotation used. Furthermore, the expression filtering allows removal of lowly expressed isoforms where the expression levels might be untrustworthy. The filtering on gene expression allows for removal of lowly expressed genes which causes the calculations of the Isoform Fractions (IF) to become untrustworthy
# removeSingleIsoformGenes = TRUE: Removal of single isoform genes is the default setting in preFilter() since these genes, per definition, cannot have changes in isoform usage.
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 0,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)
