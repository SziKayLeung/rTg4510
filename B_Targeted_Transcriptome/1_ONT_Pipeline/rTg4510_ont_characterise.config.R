## ---------- Script -----------------
##
## Purpose: input variables for characterising QC of ONT targeted mouse transcriptome datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- Packages -----------------

suppressMessages(library("dplyr"))


## ---------- LOGEN modules -----------------

LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "transcriptome_stats/read_sq_classification.R"))


## ---------- Directory and input files -----------------

dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/",
  raw = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/1_raw/F_ont_targeted/",
  meta = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/F_ont_targeted/",
  cupcake_sq = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/6b_tofu_sqanti3/"
)


# Miscellaneous input files 
misc_input <- list(
  
  # list of samples
  tg4510_samples = read.csv(paste0(dirnames$meta, "Tg4510_fullsample.csv"))[,c("Genotype","Sample.ID","RIN","ng.ul")],
  
  # list of target genes
  TargetGene = array(read.table(paste0(dirnames$meta, "TargetGenes.tsv"))[["V1"]]),
  
  # list of barcodes for each sample
  BarcodedPhenotype = read.csv(paste0(dirnames$meta, "ONTBarcoded_Phenotype.csv"))
  
)


## ---------- Annotations -----------------

# SQANTI classification files
input.class.names.files <- list(
  
  # SQANTI filtering after cupcake collapse
  cupcake_ont = paste0(dirnames$cupcake_sq, "3_sqanti3/ONT_merged_RulesFilter_result_classification.txt")
)

input.class.files <- lapply(input.class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))


# Abundance 
input.abundance.files <- list(
  
  # demultiplexed abundance after Iso-Seq3 collapse
  cupcake_ont = read.csv(paste0(dirnames$cupcake_sq, "2_collapse/ONT_merged_fl_count.csv"), row.names = "isoform")

)


# Filter sqanti by target genes, and merge with abundance
targeted.class.files <- lapply(input.class.files, function(x)  subset_class_by_targets(x, misc_input$TargetGene))
targeted.class.files$cupcake_ont <- quantify_class_abundance(targeted.class.files$cupcake_ont, input.abundance.files$cupcake_ont)