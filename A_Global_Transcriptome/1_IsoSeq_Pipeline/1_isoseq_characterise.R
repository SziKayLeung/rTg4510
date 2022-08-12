## ----------Script-----------------
##
## Script name: 
##
## Purpose of script: 
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

OUTPUT_DIR <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Whole_Transcriptome"

SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome/1_IsoSeq_Pipeline/"

source(paste0(SC_ROOT, "02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "rTg4510_isoseq_characterise.config.R"))

# results from Whole transcriptome paper
input_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/IsoSeq_Tg4510/Figures_Thesis/Whole_Transcriptome/Figures_IsoSeqWhole_InputVariables.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/IsoSeq_Tg4510/Figures_Thesis/Whole_Transcriptome/Figures_IsoSeqWhole_Functions.R")


## ---------- Apply functions ----------------

QC_yield <- QC_yield_plot() 
Lengths <- lengths_plots()
rarefaction <- rarefaction_distribution()
iso_length_plot <- iso_length(class.files)
no_of_iso_persample <- no_of_isoforms_sample(class.files)
novel_anno_plots <- novel_annotated()
exon_length_corr <- exon_length_isoform_correlation()
ERCC_plots <- run_ERCC(Ercc.class.file)
AS_genes <- AS_genes_events()
NMDvsnon <- NMD_vs_NonNMD() 
venn_IR_NMD <- IR_NMD_run(class.files)
lncRNA_plots <- lncRNA()
humanMAPT_plots <- find_mapt()
rnaseq_isoseq_counts_mouse <- rnaseq_isoseq_counts(class.files)
rnaseqisoseq <- rnaseq_isoseq_transcriptome(cuffrefmap_input,cufftmap_input)


## ---------- Output ----------------

pdf (paste0(output_plot_dir,"/IsoSeqWholeTranscriptome.pdf"), width = 10, height = 15)
plot_grid(rarefaction[[1]],iso_length_plot,NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2)
plot_grid(no_of_iso_persample[[1]],no_of_iso_persample[[2]],NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2)
plot_grid(ERCC_plots[[2]],ERCC_plots[[3]],NULL, NULL,NULL,NULL,ncol = 2,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(novel_anno_plots[[1]],novel_anno_plots[[2]],novel_anno_plots[[3]],novel_anno_plots[[4]],novel_anno_plots[[5]],novel_anno_plots[[6]],labels = c("A","B","C","D","E","F"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, nrow = 3)
plot_grid(AS_genes[[1]],AS_genes[[2]],NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2, rel_widths = c(0.62,0.38))
plot_grid(NMDvsnon,grobTree(venn_IR_NMD),NULL,NULL,NULL,NULL,labels = c("A","B"),label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, ncol = 2, rel_widths = c(0.6,0.4))
grobs <- ggplotGrob(lncRNA_plots[[3]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid(lncRNA_plots[[1]],lncRNA_plots[[2]],lncRNA_plots[[5]],lncRNA_plots[[6]],lncRNA_plots[[7]],legend,labels = c("A","B","C","D","E"),label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, nrow = 3)
plot_grid(humanMAPT_plots[[1]],humanMAPT_plots[[2]],NULL,NULL,NULL,NULL,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
grobs <- ggplotGrob(rnaseqisoseq[[6]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid(rnaseqisoseq[[2]],rnaseqisoseq[[5]],rnaseqisoseq[[6]] + theme(legend.position = "none"),legend,NULL,NULL, labels = c("A","B","C"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
dev.off()


# Chapter 5 plots
top = plot_grid(QC_yield[[1]],QC_yield[[2]], labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", nrow = 1, scale = 0.9)
mid = plot_grid(QC_yield[[3]],NULL, rel_widths = c(1,0), labels = c("C",""), label_size = 30, label_fontfamily = "CM Roman", nrow = 1, scale = 0.9)
bottom = plot_grid(QC_yield[[4]],Lengths,labels = c("D","E"), label_size = 30, label_fontfamily = "CM Roman", nrow = 1, scale = 0.9)
pdf (paste0(output_plot_dir,"/rTg4510WholeTranscriptome.pdf"), width = 10, height = 15)
plot_grid(top, mid, bottom, nrow = 3)
dev.off()
