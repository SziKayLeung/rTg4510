# Szi Kay Leung
# Functions script for Thesis Chapter on Whole Transcriptome IsoSeq

# packages
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
suppressMessages(library(data.table))
suppressMessages(library(readxl))

detach("package:plyr")
library(extrafont)
#font_install('fontcm')
loadfonts()

output_helpfig_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures"
output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Whole_Transcriptome"
# results from Whole transcriptome paper
input_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Whole_Transcriptome/Figures_IsoSeqWhole_InputVariables.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Whole_Transcriptome/Figures_IsoSeqWhole_Functions.R")

rarefaction_files()
ERCC_sqanti_files()
rnaseq_sqanti_files()
lncrna_class_files()

### 
#QC_yield <- QC_yield_plot() 
#Mapping_Stats <- Mapping_stats_plots()
#Lengths <- lengths_plots()
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

#bottom_row <- plot_grid(QC_yield[[4]], Mapping_Stats, labels = c('B', 'C'), label_size = 30, label_fontfamily = "CM Roman", ncol = 2)

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
#plot_grid(QC_yield[[1]],QC_yield[[2]],labels = "auto", label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
#plot_grid(QC_yield[[3]], bottom_row,labels = c('a', ''), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9, nrow = 2)
#plot_grid(Lengths[[1]],Lengths[[2]],Lengths[[3]], labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.8),
dev.off()