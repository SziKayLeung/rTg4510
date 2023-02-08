## ----------Script-----------------
##
## Purpose: Compare and characterise Iso-Seq vs ONT targeeted datasets
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##


## ---------- Source function and config files -----------------

# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "aesthetics_basics_plots/draw_density.R"))

sapply(list.files(path = paste0(LOGEN,"compare_datasets"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)

# project related scripts and functions
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/1_ONT_Pipeline/"
source(paste0(SC_ROOT, "02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "rTg4510_ont_characterise.config.R"))

# output directory
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"


## ---------- Transcriptome stats --------------------------------------

# Number of isoforms 
pPreFil_ONT_Num = total_num_iso(input.class.files$ONT_unfiltered,"ONT Pre-Filtered") + theme(legend.position = "none")
pIsoSeq_Num = total_num_iso(input.class.files$iso,"")
pONT_Num = total_num_iso(input.class.files$ont_retained,"")
pMerged_Num = total_num_iso(input.class.files$merged_noISM,"")

# Cage peaks, TTS, TSS
histpeaks = wrapper_hist_peaks_2datasets(input.class.files$iso, input.class.files$ont_retained,"IsoSeq","ONT")

# Iso-Seq vs ONT, and merged transcriptome 
IsoSeq_ONT_p = IsoSeq_vs_ONT_descriptive()
IsovsONT_p = characterise_iso_ont_merge()

# compare number of isoforms to known isoforms
ONT_isonum_diff = compare_isoform(ONT_retained_class,"ONT")
IsoSeq_isonum_diff = compare_isoform(IsoSeq_filtered_class,"Iso-Seq")
Merged_isonum_diff = compare_isoform(Merged_noISM,"Merged")
isonum_diff_corr = cbind(IsoSeq_isonum_diff$Corr,ONT_isonum_diff$Corr,Merged_isonum_diff$Corr)
row.names(isonum_diff_corr) = c("Gencode Isoform Number","Gene Length","Transcript Length","Gene Expression","Exon Number")

