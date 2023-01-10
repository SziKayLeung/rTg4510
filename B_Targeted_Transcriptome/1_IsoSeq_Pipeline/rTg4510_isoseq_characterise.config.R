## ---------- Script -----------------
##
## Purpose: input variables for Iso-Seq targeted mouse transcriptome datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- Directory and input files -----------------

dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/",
  meta = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/B_isoseq_targeted",
  whole_sq = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/0_archive/DiffAnalysis_noRNASEQ/SQANTI3/",
  subset_targeted_sq = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/10_characterise/subset/SQANTI3_noRNASEQ/"
)

# Phenotypes
twomos <- c("K18","O18","S18", "K17","M21","Q21")
control_samples <- c("K19","K23","K21","K17","S19","M21","O23","P19","Q21","S23","Q17","Q23")
case_samples <- c("K18","K20","K24","L22","O18","O22","T20","Q20","S18","Q18","L18","T18")

# Input IsoSeq Paths
spec_dirnames <- list(
  ccs = paste0(dirnames$root, "1_ccs"),
  lima = paste0(dirnames$root, "2_lima"),
  refine = paste0(dirnames$root, "3_refine"),
  cluster = paste0(dirnames$root, "4_cluster")
)


# Miscellaneous input files 
misc_input <- list(
  
  # ccs and lima output
  ccs_output = read.csv(paste0(spec_dirnames$ccs, "/AllMouseTargeted_CCS_output.csv"), header = T),
  lima_output = read.csv(paste0(spec_dirnames$lima, "/AllMouseTargeted_LIMA_summary.csv"), header = T),
  
  # list of target genes
  TargetGene = read.table(paste0(dirnames$meta, "/TargetGenes.tsv"))[["V1"]],
  
  # phenotype
  targetedpheno = read.csv(paste0(dirnames$meta, "/Targeted_Sample_Demographics.csv")),
  
  # list of samples
  tg4510_samples = read.csv(paste0(dirnames$meta, "/Tg4510_fullsample.csv"))[,c("Genotype","Sample.ID","RIN","ng.ul")],
  
  # comparison from gffcompare 
  cuff_tmap = read.table(paste0(dirnames$subset_targeted_sq, "Whole_Targeted.SubsetAllMouseTargeted.collapsed_classification.filtered_lite.gtf.tmap"), header = T)
  
)


## ---------- SQANTI classification files -----------------

input.class.names.files <- list(
  whole_iso = paste0(dirnames$whole_sq, "WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt"),
  sub_targeted_iso = paste0(dirnames$subset_targeted_sq, "SubsetAllMouseTargeted.collapsed_classification.filtered_lite_classification.txt") 
)

input.class.files <- lapply(input.class.names.files, function(x) SQANTI_class_preparation(x,"standard"))
input.class.files <- lapply(input.class.files, function(x) SQANTI_remove3prime(x))
input.class.files <- lapply(input.class.files, function(x) 
  x %>% mutate(TargetGene = ifelse(associated_gene %in% misc_input$TargetGene, "Target Genes","Not Target Genes")))


## ---------- Target Rate -----------------

Probes_files <- read_target_probes(paste0(dirnames$root,"/6b_target_rate"), "fasta.sam.probe_hit.txt")


## ---------- Transgene sequence -----------------

# headers of the clustered read names that contain human MAPT sequences
# Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
humanmapt_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/10_characterise/transgene"
hMAPT_input <- read.table(paste0(humanmapt_input_dir,"/hMAPT_counts.txt")) %>% `colnames<-`(c("Sample", "Counts"))
mMAPT_input <- read.table(paste0(humanmapt_input_dir,"/hMAPT_counts.txt")) %>% `colnames<-`(c("Sample", "Counts"))
cluster_reads <- list.files(CLUSTER_input_dir, pattern = "cluster_report.csv", full.names = T)
hMAPT_input_all <- read.table(paste0(humanmapt_input_dir,"/AllMouseTargeted_Ids.txt"))