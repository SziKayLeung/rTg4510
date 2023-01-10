## ----------Script-----------------
##
## Purpose: Generate post-sqanti plots for summary of Iso-Seq global mouse transcriptome datasets
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
sapply(list.files(path = paste0(LOGEN,"alternative_splicing"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"compare_datasets"), pattern="*.R", full = T), source,.GlobalEnv)

# project related scripts and functions
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome/1_IsoSeq_Pipeline/"
source(paste0(SC_ROOT, "02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "rTg4510_isoseq_characterise.config.R"))

# output directory
output_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Whole_Transcriptome"


## ---------- QC output --------------------------------------

# Number of Iso-Seq reads
Reads <- number_iso_reads(misc_input$sample_run,
                         misc_input$ccs_output,
                         misc_input$lima_output,
                         spec_dirnames$refine,
                         spec_dirnames$cluster,
                         misc_input$CLUSTER_merge)

Reads$Reads$Genotype <- sapply(Reads$Reads$sample, function(x) classify_genotype(x, case_samples, control_samples,"J20"))
Reads$CCS_values$Genotype <- sapply(Reads$CCS_values$sample, function(x) classify_genotype(x, case_samples, control_samples,"J20"))
Reads$Reads <- Reads$Reads %>% filter(Genotype != "J20")
write.csv(Reads$Reads,paste0(output_dir,"/Tg4510_IsoSeqReadsStats.csv"))

# Plots relating to the number of Iso-Seq reads
QC_yield <- iso_QC_yield(Reads$Reads, misc_input$sequenced)

# Lengths table output
all_clustered <- tab_raw_isoseq_length(paste0(spec_dirnames$cluster,"/Lengths"),"clustered.hq.fasta.seqlengths.txt","clustered")
write.csv(all_clustered,paste0(output_dir,"/Tg4510_CLusteredLengthsReadsStats.csv"))

all_ccs <- tab_raw_isoseq_length(paste0(spec_dirnames$ccs,"/Lengths"),"ccs.fasta.seqlengths.txt","CCS")
write.csv(all_ccs,paste0(output_dir,"/Tg4510_CCSLengthsReadsStats.csv"))

CCS_lengths <- plot_raw_isoseq_length(all_ccs,"CCS")

# Rarefaction curve
rarefaction <- rarefaction_distribution(misc_input$rarefaction)


## ---------- Transcriptome summarise stats -----------------

iso_length_plot <- multi_isoform_length(input.class.files$isoseq)
no_of_iso_persample <- no_of_isoforms_sample(input.class.files$isoseq)
novel_anno_plots <- diff_isoform_type_by_features(input.class.files$isoseq)
exon_length_corr <- corr_exon_length_num(input.class.files$isoseq)
lncRNA_plots <- differentiate_lncRNA(input.class.files$isoseq, input.class.files$lncRNA)

## ---------- Compare datasets ------------------------------

ERCC_plots <- number_ERCC_detected(input.class.files$ercc, misc_input$ercc_cal)
rnaseq_isoseq_counts_mouse <- correlate_rnaseq_isoseq_counts(input.class.files$isoseq, misc_input$cuffrefmap_input, misc_input$RNASeq_Def, misc_input$IsoSeq_Def)
rnaseqisoseq <- compare_rna_vs_iso_transcriptome(input.class.files$isoseq, input.class.files$rnaseq, misc_input$cuffrefmap_input, misc_input$cufftmap_input)


## ---------- Alternative splicing  -------------------------

AS_genes <- AS_genes_events(misc_input$AS_IR_genes, misc_input$AS_IR_isoforms, misc_input$AS_IR_knowngenes, "Mouse")
venn_IR_NMD <- associate_IR_NMD(input.class.files$isoseq)
NMDvsnon <- plot_IR_NMD_isoform_expression(input.class.files$isoseq) 


humanMAPT_plots <- find_mapt()


## ---------- Pdf Output ----------------------------------------

pdf(paste0(output_dir,"/IsoSeqWholeTranscriptome.pdf"), width = 10, height = 15)
plot_grid(rarefaction[[1]],iso_length_plot,NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2)
plot_grid(no_of_iso_persample[[1]],no_of_iso_persample[[2]],NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2)
plot_grid(ERCC_plots[[2]],ERCC_plots[[3]],NULL, NULL,NULL,NULL,ncol = 2,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(novel_anno_plots[[1]],novel_anno_plots[[2]],novel_anno_plots[[3]],novel_anno_plots[[4]],novel_anno_plots[[5]],novel_anno_plots[[6]],labels = c("A","B","C","D","E","F"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, nrow = 3)
plot_grid(AS_genes[[1]],AS_genes[[2]],NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2, rel_widths = c(0.62,0.38))
plot_grid(NMDvsnon,grobTree(venn_IR_NMD),NULL,NULL,NULL,NULL,labels = c("A","B"),label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2, rel_widths = c(0.6,0.4))
grobs <- ggplotGrob(lncRNA_plots[[3]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid(lncRNA_plots[[1]],lncRNA_plots[[2]],lncRNA_plots[[5]],lncRNA_plots[[6]],lncRNA_plots[[7]],legend,labels = c("A","B","C","D","E"),label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, nrow = 3)
#plot_grid(humanMAPT_plots[[1]],humanMAPT_plots[[2]],NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
grobs <- ggplotGrob(rnaseqisoseq[[6]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid(rnaseqisoseq[[2]],rnaseqisoseq[[5]],rnaseqisoseq[[6]] + theme(legend.position = "none"),legend,NULL,NULL, labels = c("A","B","C"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
dev.off()


# Chapter 5 plots
pdf(paste0(output_dir,"/rTg4510WholeTranscriptome.pdf"), width = 10, height = 15)
top = plot_grid(QC_yield[[1]],QC_yield[[2]], labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", nrow = 1, scale = 0.9)
mid = plot_grid(QC_yield[[3]],NULL, rel_widths = c(1,0), labels = c("C",""), label_size = 30, label_fontfamily = "CM Roman", nrow = 1, scale = 0.9)
bottom = plot_grid(QC_yield[[4]],CCS_lengths,labels = c("D","E"), label_size = 30, label_fontfamily = "CM Roman", nrow = 1, scale = 0.9)
plot_grid(top, mid, bottom, nrow = 3)
dev.off()
