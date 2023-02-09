## ----------Script-----------------
##
## Purpose: Generate post-sqanti plots for summary of ONT targeted mouse transcriptome datasets
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##


## ---------- Source function and config files -----------------

# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "aesthetics_basics_plots/draw_density.R"))
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)

# project related scripts and functions
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/1_ONT_Pipeline/"
source(paste0(SC_ROOT, "02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "rTg4510_ont_qc_runs.config.R"))

# output directory
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"


## ---------- QC output --------------------------------------

# Number of Reads across the Batches
ONTDemultiplexReads = number_ont_reads(sequencing_data, misc_input$ReadCounts, samples_removed, num_barcodes, 
                                       misc_input$BarcodedPhenotype,misc_input$tg4510_samples)

# Sequencing summary
AllQCPlots_B2 <- AllQCPlots(B2_sequencingdata)
AllQCPlots_B3 <- AllQCPlots(B3_sequencingdata)

# On Target Rate 
ontarget <- on_target_plot(Probes_input, misc_input$BarcodedPhenotype,"batched")


## ---------- hMAPT --------------------------------------

hmapt_iso = find_mapt_isoseq()
hmapt_ont = find_mapt_ont()


## ---------- Pdf Output ----------------

pdf(paste0(output_dir,"/ONTTargetedTranscriptome.pdf"), width = 10, height = 15)
# QC plots 1
plot_grid(AllQCPlots_B2$pchannel, AllQCPlots_B3$pchannel,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", nrow = 3,scale = 0.9)
# QC plots 2
grobs <- ggplotGrob(AllQCPlots_B2$ptemp + theme(legend.position = c(0.7,0.1),legend.direction = "horizontal"))$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid(NULL,legend,AllQCPlots_B2$ptemp + theme(legend.position = "none"),AllQCPlots_B3$ptemp + theme(legend.position = "none"),AllQCPlots_B2$ptemp_cum + theme(legend.position = "none"),AllQCPlots_B3$ptemp_cum + theme(legend.position = "none"),labels = c("","","A","B","C","D"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
# QC plots 2
grobs <- ggplotGrob(AllQCPlots_B2$plength + theme(legend.position = c(0.7,0.1),legend.direction = "horizontal"))$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid(NULL,legend,AllQCPlots_B2$pscore,AllQCPlots_B3$pscore,AllQCPlots_B2$plength + theme(legend.position = "none"), AllQCPlots_B3$plength + theme(legend.position = "none"),labels = c("","","A","B","C","D"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2)
#plot_grid(AllQCPlots_B2$pspeedtime, AllQCPlots_B3$pspeedtime,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", nrow = 2,scale = 0.9)
#plot_grid(AllQCPlots_B2$pcontour,AllQCPlots_B3$pcontour,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
plot_grid(ONTDemultiplexReads$p1,ONTDemultiplexReads$p2,ONTDemultiplexReads$p3,ONTDemultiplexReads$p4,labels = c("A","B","C","D"), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
plot_grid(ontarget,NULL,NULL, label_size = 30, label_fontfamily = "CM Roman", nrow = 3, scale = 0.9)
dev.off()