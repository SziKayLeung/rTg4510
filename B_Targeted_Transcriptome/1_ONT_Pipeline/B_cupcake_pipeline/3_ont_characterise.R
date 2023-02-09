## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Purpose: Generate post-sqanti plots for summary of ONT targeted mouse transcriptome datasets
##
## --------------------------------


## ---------- packages -----------------------------------------

suppressMessages(library("dplyr"))
suppressMessages(library("cowplot"))


## ---------- Source function and config files -----------------

# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)

# project related scripts and functions
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/1_ONT_Pipeline/"
source(paste0(SC_ROOT, "02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "rTg4510_ont_characterise.config.R"))

# output directory
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"


## ---------- QC output --------------------------------------

# sensitivity plot of ONT datasets after cupcake collapse to determine threshold
ptargetall <- plot_cupcake_collapse_sensitivity(targeted.class.files$cupcake_ont,"All 20 target genes")
ptargetapp <- plot_cupcake_collapse_sensitivity(targeted.class.files$cupcake_ont %>% filter(associated_gene == "App"), "App")
plot_grid(ptargetall[[1]],ptargetall[[2]],ptargetapp[[1]],ptargetapp[[2]], nrow = 2, ncol = 2, labels = "auto")
