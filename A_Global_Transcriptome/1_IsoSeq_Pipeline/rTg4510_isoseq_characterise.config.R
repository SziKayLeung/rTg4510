## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
## 
##   
##
##


## ---------- Directory and input files -----------------

WKD_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/"
METADATA_DIR <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/A_isoseq_whole/"
WHOLESC_DIR <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/Whole_Transcriptome_Paper/"

source(paste0(WHOLESC_DIR, "Output/SQANTI_General.R"))
source(paste0(WHOLESC_DIR, "Alternative_Splicing/AS_Events.R"))

# Phenotypes
twomos <- c("K18","O18","S18", "K17","M21","Q21")
WT <- c("K17","M21","Q21","K23","O23","S23")

# Input Sequencing
sequencing_output <- read.csv(paste0(METADATA_DIR, "Sequenced_mouse_output.csv"))
tg4510_samples <- read.csv(paste0(METADATA_DIR, "Tg4510_fullsample.csv"))[,c("Genotype","Sample.ID","RIN","ng.ul")]
sequenced <- merge(sequencing_output, tg4510_samples, by.x = "Sample", by.y = "Sample.ID", all.x = TRUE)
sequenced$Genotype <- factor(sequenced$Genotype, levels=c("WT", "TG"))

# Input IsoSeq Paths
CCS_input_dir <- paste0(WKD_ROOT, "IsoSeq/CCS")
LiMA_input_dir <- paste0(WKD_ROOT, "IsoSeq/LIMA")
REFINE_input_dir <- paste0(WKD_ROOT, "IsoSeq/REFINE")
CLUSTER_input_dir <- paste0(WKD_ROOT, "IsoSeq/CLUSTER")
CLUSTER_merge <- read.csv(paste0(WKD_ROOT, "IsoSeq/MERGED_CLUSTER/WholeIsoSeq.clustered.cluster_report.csv"))

# Input Lengths Paths
lengths_input_dir <- paste0(WKD_ROOT, "CCS/Lengths")

# sample run 
sample_run <- read.csv(paste0(WHOLESC_DIR, "IsoSeq_Processing/Mouse/All_Tg4510_Demultiplex.csv"))


## ---------- SQANTI classification file -----------------

# Iso-Seq 
class.names.files <- paste0(WKD_ROOT, "Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.txt")
class.files <- SQANTI_class_preparation(class.names.files,"standard")

# ERCC 
ERCC.class.names.files <- paste0(WKD_ROOT, "ERCC/WholeIsoSeq_sqantitamafiltered.classification.txt")
Ercc.class.file <- read.table(ERCC.class.names.files, header=T, as.is=T, sep="\t")

# lncRNA 
lncRNA.class.names.files <- paste0(WKD_ROOT, "Post_IsoSeq/SQANTI_TAMA_FILTER/LNCRNA/WholeIsoSeq_sqantitamafiltered.classification.txt")
lnc.files <- SQANTI_class_preparation(lncRNA.class.names.files,"standard")

# RNASeq Stringtie Transcriptome
rnaseq_transcriptome <- paste0(WKD_ROOT, "RNASeq/SQANTI2/rnaseq_stringtie_merged_final_classification.filtered_lite_classification.txt")
rnaseq.class.files <- SQANTI_class_preparation(rnaseq_transcriptome,"rnaseq") %>% mutate(associated_gene = toupper(.$associated_gene))
rnaseq.class.files$Sample <- "RNA-Seq"


## ---------- Mapping statistics  -----------------
# PAF mapping stats 
map <- read.table(paste0(WKD_ROOT,"Post_IsoSeq/MAP/WholeIsoSeq_reads_with_alignment_statistics.txt"))
colnames(map) <- c("name_of_read","chrom","start","read_length","alignment_length",
                   "alignment_length_perc","length_of_target_seq_alignment","alignment_identity",
                   "strand","num_mismatches","num_insertions","num_deletions")


## ---------- Rarefaction  -----------------

rarefaction_dir <- paste0(WKD_ROOT,"Post_IsoSeq/RAREFACTION/")
read_rarefaction_files(rarefaction_dir, "WholeIsoSeq")


## ---------- Transgene  -----------------

# Human MAPT inputs
# headers of the clustered read names that contain human MAPT sequences
# Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
humanmapt_input_dir <- paste0(WKD_ROOT, "Post_IsoSeq/HUMANMAPT") 
hMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "hMAPT2.header", full.names = T)
mMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "mMAPT1.header", full.names = T)
cluster_reads <- list.files(CLUSTER_input_dir, pattern = "cluster_report.csv", full.names = T)
hMAPT_input_all <- read.table(paste0(humanmapt_input_dir,"/WholeIsoSeq_Ids.txt"))


## ---------- Kallisto alignment  -----------------

RNASeq_Def <- read.table(paste0(WKD_ROOT, "Post_IsoSeq/KALLISTO/WholeIsoSeq.abundance.tsv"), header = T) 
IsoSeq_Def <- read.table(paste0(WKD_ROOT, "Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/KALLISTO/WholeIsoSeq.abundance.tsv"), header = T) %>% 
  mutate(target_id = word(target_id, c(1), sep = fixed("<")))


## ---------- Gffcompare: RNA-Seq vs Iso-Seq  -----------------

# Gffcompare output
cuff_suffix <- "long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf"
cuffrefmap_input <- paste0(WKD_ROOT, "RNASeq/SQANTI2/", cuff_suffix, ".refmap")
cufftmap_input <- paste0(WKD_ROOT, "RNASeq/SQANTI2/", cuff_suffix, ".tmap")
