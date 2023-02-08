## ---------- Script -----------------
##
## Purpose: input variables for differential analysis of Iso-Seq targeted mouse transcriptome datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("readxl"))


# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)

## 
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")

# output directory
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"

## ---------- TappAS input files -----------------

dirnames <- list(
  meta = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/",
  iso = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/thesis_dump/DiffAnalysis_noRNASEQ/",
  ont = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/thesis_dump/TALON/",
  merged = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/"
)


# output files from running tappAS
tappas_dirnames = list(
  iso = paste0(dirnames$iso,"TAPPAS_OUTPUT/IsoSeq_Expression/"),
  ont = paste0(dirnames$ont,"TAPPAS_OUTPUT/")
)

# miscellaneous input
misc_input <- list(
  
  ont_sq_reason = read.table(paste0(dirnames$ont, "All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon_classification.filtered_lite_reasons.txt"), 
                             header = T, sep = ","), 
  
  iso = read.table(paste0(dirnames$meta,"B_isoseq_targeted/TargetedMouse_PhenotypeTAPPAS.txt"), header = T),
  
  ont = read.table(paste0(dirnames$meta,"F_ont_targeted/ONT_phenotype.txt"), header = T)

)


## ---------- SQANTI classification files -----------------

### decide which to use: filtered or unfiltered
# Classification file
class.names.files <- list(
  iso = paste0(dirnames$iso, "SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"),
  ont = paste0(dirnames$ont, "All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon_classification.filtered_lite_classification.txt"),
  ont_unfil = paste0(dirnames$ont, "All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon_classification.txt"),
  merged = paste0(dirnames$merged, "3_sqanti3/IsoSeqONT_final_genename_classification_noISM.txt")
)

input.class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

## ---------- TappAS output -----------------

tappassiggene <- lapply(list("iso" = paste0(dirnames$iso,"/TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx"),
                             "ont" = paste0(dirnames$ont,"TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx")), 
                        function(x) read_dea_files(x))

tappassigtrans <- lapply(list("iso" = paste0(dirnames$iso,"/TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx"),
                              "ont" = paste0(dirnames$ont,"TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx")), 
                         function(x) read_dea_files(x))

# sort by R-sqaured
tappassigtrans$ont$TargetedOnt_Transexp <- tappassigtrans$ont$TargetedOnt_Transexp %>% arrange(-`R-squared`)