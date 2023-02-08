## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: sources functions for generating downstream plots for SFARI dataset 
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

## ---------- Packages -----------------

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(rjson)) # json files
suppressMessages(library(plyr)) # revalue
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(stringr)) 
suppressMessages(library(viridis)) 
suppressMessages(library(wesanderson)) 
suppressMessages(library(extrafont))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(tibble))
suppressMessages(library(VennDiagram))
suppressMessages(library(directlabels))
suppressMessages(library(cowplot))
suppressMessages(library(readxl))
suppressMessages(library(ggdendro))
suppressMessages(library(pheatmap))
suppressMessages(library(extrafont))
suppressMessages(loadfonts())


## ----------Functions-----------------

# load all the functions
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/plot_aesthetics.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/plot_tappas_analysis.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/sqanti_general.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/comp_characterise.R")

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT,"bin/segregate_tappasresults.R"))


## ----------Theme-----------------

label_colour <- function(genotype){
  if(genotype %in% c("WT","Control")){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "WT_2mos"){colour = alpha(wes_palette("Royal1")[2],0.5)}else{
      if(genotype %in% c("TG","Case")){colour = wes_palette("Royal1")[2]}else{
        if(genotype == "TG_2mos"){colour = alpha(wes_palette("Royal1")[1],0.5)}else{
          if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
            if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
              if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
              }}}}}}}
  return(colour)
}

label_group <- function(genotype){
  if(genotype %in% c("Case","CASE")){group = "TG"}else{
    if(genotype %in% c("Control","CONTROL")){group = "WT"}}
  return(group)
}

## ---------- Load tappAS files -----------------

loaded <- list(
  glob_iso = input_tappasfiles(tappas_dir$glob_iso),
  glob_rna = input_tappasfiles(tappas_dir$glob_rna),
  targ_iso = input_tappasfiles(tappas_dir$targ_iso),
  targ_ont = input_tappasfiles(tappas_dir$targ_ont)
)

## ---------- Annotate tappAS files -----------------

annotated <- list(
  glob_iso = annotate_tappasfiles(class.files$glob_iso,loaded$glob_iso$input_normalized_matrix,phenotype$glob_iso),
  glob_rna = annotate_tappasfiles(class.files$glob_iso,loaded$glob_rna$input_normalized_matrix,phenotype$glob_rna),
  targ_iso = annotate_tappasfiles(class.files$targ_iso,loaded$targ_iso$input_normalized_matrix,phenotype$targ_iso),
  targ_ont = annotate_tappasfiles(class.files$targ_ont,loaded$targ_ont$input_normalized_matrix,phenotype$targ_ont)
)


## ---------- Differential expression models -----------------
# segregate differential gene and transcript expression results by the different models (using beta coefficient as filters)
diff_models <- list(
  glob_iso_siggenes = segregate_tappasresults(tappassiggene$glob$WholeIso_Genexp,"IsoSeq"),
  glob_rna_siggenes = segregate_tappasresults(tappassiggene$glob$WholeRNA_Genexp,"RNASeq")
)

