library(dplyr)
library(IsoformSwitchAnalyzeR)

install.packages("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/IsoformSwitchAnalyzeR_1.13.03.tar.gz", repos = NULL, type="source")
install.packages("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tximeta_1.8.2.tar.gz", repos = NULL, type="source")

# read in SQANTI2 classification file of all merged data
input_dir <- "/gpfs/ts0/home/sl693/files/"
# sqanti classification file filtering removed monoexonic transcripts 
class <- read.table(paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite_classification.txt"),header=T, as.is = T, sep = "\t")

# PacBio isoform ID, keep only columns with isoform FL counts (and convert FL counts to TPM)
counts <- class[, c("isoform","FL.K17","FL.K18","FL.K24","FL.L22","FL.M21","FL.O18","FL.O23","FL.Q20","FL.Q21","FL.S18","FL.S23")] 
abundance <- cbind(counts[,1], as.data.frame(apply(counts[,-1],2, function(x) round(x/sum(x)*1000000,0))))
colnames(counts)[1] <- c("isoform_id")      # column name important for input into IsoformSwitch
colnames(abundance)[1] <- c("isoform_id")


# design dataframe input for IsoformSwitchAnalyzeR
myDesign <- data.frame(
  sampleID = colnames(counts)[-1],
  condition = c("WT", "TG", "TG", "TG","WT","TG","WT","TG","WT","TG","WT")
)

# using SQANTI generated gtf and fasta, rather than gencode
#isoformExonAnnoation = paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite.gtf"),
aSwitchList <- importRdata(
  isoformCountMatrix   = counts,
  isoformRepExpression = abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = paste0(input_dir,"edited.gtf"),
  #isoformNtFasta       = paste0(input_dir,"WT8_all_samples.chained_classification.filtered_lite.fasta"),
)

edited_gtf <- rtracklayer::import("/gpfs/ts0/home/sl693/files/edited.gtf")
original_gtf <- rtracklayer::import("/gpfs/ts0/home/sl693/files/original.gtf")

edited_gtf[edited_gtf$transcript_id == "PB.429.72"]
original_gtf[original_gtf$transcript_id == "PB.429.72"]