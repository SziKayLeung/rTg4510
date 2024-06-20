## ----------Script-----------------
##
## Purpose: Generate post-sqanti plots for summary of Iso-Seq targeted mouse transcriptome datasets
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##


## ---------- Source function and config files -----------------

# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "aesthetics_basics_plots/draw_density.R"))
source(paste0(LOGEN, "transcriptome_stats/read_sq_classification.R"))

sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"compare_datasets"), pattern="*.R", full = T), source,.GlobalEnv)

# project related scripts and functions
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/"
source(paste0(SC_ROOT, "02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "rTg4510_isoseq_characterise.config.R"))

# output directory
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"


## ---------- QC output --------------------------------------

# Number of Iso-Seq reads
Reads <- number_iso_reads(misc_input$ccs_output,
                          misc_input$lima_output,
                          spec_dirnames$refine,
                          spec_dirnames$cluster)

Reads$Reads$Genotype <- sapply(Reads$Reads$sample, function(x) classify_genotype(x, case_samples, control_samples,"Batch"))
Reads$CCS_values$Genotype <- sapply(Reads$CCS_values$sample, function(x) classify_genotype(x, case_samples, control_samples,"Batch"))
write.csv(Reads$Reads,paste0(output_dir,"/Tg4510_IsoSeqTargetedReadsStats.csv"))

# Plots relating to the number of Iso-Seq reads
QC_yield <- iso_QC_yield_batch(Reads$Reads, misc_input$targetedpheno, misc_input$tg4510_samples)

# Target rate of probes
ontarget <- on_target_plot(Probes_files, misc_input$targetedpheno,"notbatched")

# Transgene sequence (Iso-Seq)
humanMAPT_plots <- find_mapt(dirnames$mapt, misc_input$targeted_tg4510_samples)

## ---------- Compare datasets ------------------------------

# Whole Transcriptome vs Targeted Transcriptome plots 
wholevstargeted <- whole_vs_targeted_plots(misc_input$cuff_tmap, misc_input$TargetGene,input.class.files$whole_isoseq, input.class.files$subsettargeted_isoseq)


## ---------- Pdf Output ----------------

pdf(paste0(output_dir,"/TargetedTranscriptome.pdf"), width = 10, height = 15)
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