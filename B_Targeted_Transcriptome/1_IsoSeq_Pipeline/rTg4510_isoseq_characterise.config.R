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

WKD_ROOT="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted"

METADATA_DIR <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/B_isoseq_targeted"
TargetGene <- read.table(paste0(METADATA_DIR, "/TargetGenes.tsv"))[["V1"]]
targetedpheno <- read.csv(paste0(METADATA_DIR, "/Targeted_Sample_Demographics.csv")) 
tg4510_samples <- read.csv(paste0(METADATA_DIR, "/Tg4510_fullsample.csv"))[,c("Genotype","Sample.ID","RIN","ng.ul")]

CCS_input_dir <- paste0(WKD_ROOT,"/1_ccs")
LiMA_input_dir <- paste0(WKD_ROOT,"/2_lima")
REFINE_input_dir <- paste0(WKD_ROOT,"/3_refine")
CLUSTER_input_dir <- paste0(WKD_ROOT,"/4_cluster")
#CLUSTER_Merge <- read.csv(paste0(WKD_ROOT,"/5_merged_cluster/All_Targeted_Merged.clustered.cluster_report.csv"))


## ---------- SQANTI classification files -----------------

# whole transcriptome
whole_sq_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/DiffAnalysis_noRNASEQ/SQANTI3/"
whole.class.names.files = paste0(whole_sq_dir, "WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt")
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard") %>% 
  mutate(ss_category = paste0(structural_category,"_", subcategory)) %>%
  filter(ss_category != "ISM_3prime_fragment") %>% 
  mutate(ADGene = ifelse(associated_gene %in% TargetGene, "Target Genes","Not Target Genes")) %>%
  mutate(within_50_cage = ifelse(abs(dist_to_cage_peak) <= 50, "Within 50bp","Not within 50bp"))

# subsetted target transcriptome
subset_targeted_sq_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/10_characterise/subset/SQANTI3_noRNASEQ/"
subsettargeted.class.names.files <- paste0(subset_targeted_sq_dir, "SubsetAllMouseTargeted.collapsed_classification.filtered_lite_classification.txt") 
subsettargeted.class.files <- SQANTI_class_preparation(subsettargeted.class.names.files,"standard") %>% 
  mutate(ss_category = paste0(structural_category,"_", subcategory)) %>%
  filter(ss_category != "ISM_3prime_fragment") %>%
  mutate(ADGene = ifelse(associated_gene %in% TargetGene, "Target Genes","Not Target Genes"))

# comparison from gffcompare 
cuff_tmap = read.table(paste0(subset_targeted_sq_dir, "Whole_Targeted.SubsetAllMouseTargeted.collapsed_classification.filtered_lite.gtf.tmap"), header = T)


## ---------- Target Rate -----------------

Probes_input <- list.files(path = paste0(WKD_ROOT, "/6b_target_rate"), pattern = "fasta.sam.probe_hit.txt", full.names = T)
Probes_files <- lapply(Probes_input, function(x) read.table(x, header=T, as.is=T, sep="\t"))
names(Probes_files) <- list.files(path = paste0(WKD_ROOT, "/6b_target_rate"), pattern = "fasta.sam.probe_hit.txt")


## ---------- Transgene sequence -----------------

# headers of the clustered read names that contain human MAPT sequences
# Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
humanmapt_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/10_characterise/transgene"
hMAPT_input <- read.table(paste0(humanmapt_input_dir,"/hMAPT_counts.txt")) %>% `colnames<-`(c("Sample", "Counts"))
mMAPT_input <- read.table(paste0(humanmapt_input_dir,"/hMAPT_counts.txt")) %>% `colnames<-`(c("Sample", "Counts"))
cluster_reads <- list.files(CLUSTER_input_dir, pattern = "cluster_report.csv", full.names = T)
hMAPT_input_all <- read.table(paste0(humanmapt_input_dir,"/AllMouseTargeted_Ids.txt"))

