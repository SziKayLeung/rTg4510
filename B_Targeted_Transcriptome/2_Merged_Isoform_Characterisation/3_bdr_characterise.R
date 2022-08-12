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

OUTPUT_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/01_figures_tables"

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR/1_IsoSeq_Pipeline/02_source_characterise_functions.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR/1_IsoSeq_Pipeline/bdr_characterise.config.R")


## ---------- Apply functions ----------------
# QC plots of Iso-Seq pipeline
QC_yield <- QC_yield_plot() 

# Number of isoforms detected per target gene
sq_num <- final_num_iso(targeted.class.files)

# Target rate of probes
ontarget <- on_target_plot()

# Histogram peaks of cage, TSS and TSS peaks
histpeaks <- all_hist_peaks(Control_targeted,AD_targeted)

# Descriptive plots
desc_p <- two_groups_descriptions(Control_targeted,AD_targeted,Rnaseq_gene_counts)


## ---------- Output ----------------
pdf(paste0(OUTPUT_DIR,"/BDRTargetedTranscriptome.pdf"), width = 10, height = 15)
# pg 1 - QC
bottom_row <- plot_grid(QC_yield[[2]],QC_yield[[3]],labels = c('b', 'c'), label_size = 30, label_fontfamily = "CM Roman", ncol = 2,scale = 0.9)
plot_grid(QC_yield[[1]],bottom_row,labels = c('a', ''), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)

# pg 2 - Target rate
plot_grid(ontarget,NULL,NULL, label_size = 30, label_fontfamily = "CM Roman", nrow = 3, scale = 0.9)

# pg 3 - Number of isoforms, rnaseq vs isoseq expression
generate_cowplot(sq_num[[2]], desc_p[[10]],NULL,NULL,num=2,nrow=2,ncol=2)

# pg 4 - Descriptive plots
generate_cowplot(desc_p[[1]],desc_p[[2]],
                 histpeaks[[1]],histpeaks[[2]] + theme(legend.position = "None"),
                 histpeaks[[3]] + theme(legend.position = "None"),desc_p[[7]],num=6,nrow=4,ncol=2)
dev.off()


## ---------- Annotations ----------------
pTREM2 = all_t_plots("TREM2")
pABCA1 = all_t_plots("ABCA1")
dendro_plot("MAPT")
ES_tally_plots("MAPT")
A5A3_tally_plots("MAPT")



