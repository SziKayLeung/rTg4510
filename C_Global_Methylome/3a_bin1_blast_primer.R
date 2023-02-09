## ---------- Script -----------------
##
## Purpose of script: Identify Transcript ID associated with Bin1 isoforms for creating a blast library 
##
## ---------- Notes -----------------
## Prerequisite for downstream script (3b_bin1_blast_primer.sh)
## To generate Transcript ID for subsetting and to generate a fasta file to then create a blast library
## Output file: bin1_isoforms_classification.txt 
##
##


## ---------- Source function and config files ---------------

suppressMessages(library("dplyr"))

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/sqanti_general.R")

# Specific analysis related
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/"
source(paste0(SC_ROOT,"B_Targeted_Transcriptome/3_Differential_Analysis/rTg4510_differential.config.R"))

output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/G_Bin1_Integration/BLAST"


## ---------- Subset bin1 isoforms and write output ---------------
bin1_isoforms <- class.files$ont %>% filter(associated_gene == "Bin1") %>% .[,"isoform"]
write.table(bin1_isoforms, paste0(output_dir,"/bin1_isoforms_classification.txt"), quote = F, row.names = F, col.names = F)
