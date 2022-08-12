## ----------Script-----------------
##
## Script name: 
##
## Purpose of script: Functions script for ADBDR dataset 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ----------Notes-----------------
##
## 
##   
##
##


## ---------- Source function and config files -----------------

OUTPUT_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/02_source_characterise_functions.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/rTg4510_isoseq_characterise.config.R")


## ---------- Apply functions ----------------
# QC plots of Iso-Seq pipeline
QC_yield <- QC_yield_plot() 

# Target rate of probes
ontarget <- on_target_plot()

# Whole Transcriptome vs Targeted Transcriptome plots 
wholevstargeted <- whole_vs_targeted_plots()


## ---------- Output ----------------
pdf (paste0(OUTPUT_DIR,"/TargetedTranscriptome.pdf"), width = 10, height = 15)
# pg 1 - QC
top_row <- plot_grid(QC_yield[[1]], labels = c('A'), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
bottom_row <- plot_grid(QC_yield[[3]],QC_yield[[2]],labels = c('B', 'C'), label_size = 30, label_fontfamily = "CM Roman", ncol = 2,scale = 0.9)
plot_grid(top_row, bottom_row, nrow = 2)

# pg 2 - Target rate
plot_grid(ontarget,NULL,NULL, label_size = 30, label_fontfamily = "CM Roman", nrow = 3, scale = 0.9)

# pg 3 - Whole vs Targeted Transcriptome
top_row <- plot_grid(wholevstargeted[[1]],wholevstargeted[[2]],
                     labels = c('A', 'B'), label_size = 30, label_fontfamily = "CM Roman", ncol = 2,scale = 0.9, rel_widths = c(0.6,0.4))
bottom_row <- plot_grid(wholevstargeted[[3]],wholevstargeted[[4]],
                        labels = c('C', 'D'), label_size = 30, label_fontfamily = "CM Roman", ncol = 2,scale = 0.9)
plot_grid(top_row, bottom_row, nrow = 2)
dev.off()