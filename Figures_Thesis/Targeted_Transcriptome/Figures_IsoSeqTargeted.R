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

# Adapted from report_detail.R (https://github.com/PacificBiosciences/barcoding)
# 19/01/2021: Lima report output on plot 1, 10, 17 from report_detail.R

#detach("package:plyr")
library(extrafont)
#font_install('fontcm')
loadfonts()

output_helpfig_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures"
output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Targeted_Transcriptome"
# results from Whole transcriptome paper
input_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Targeted_Transcriptome/Figures_IsoSeqTargeted_InputVariables.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Targeted_Transcriptome/Figures_IsoSeqTargeted_Functions.R")

Probe_file()
Filtered_data_staged()
whole_vs_targeted_plots()
sqanti_filter_reason()
find_mapt()

### 
QC_yield <- QC_yield_plot() 
ontarget <- on_target_plot()
postisofilter <- level2filter() 
sq_num <- final_num_iso(targeted.class.files)
wholevstargeted <- whole_vs_targeted_plots()
targeted_filtered <- sqanti_filter_reason()
sqantifilter_plots <- sqantifilter_valdidation()

### Differential Analysis 
tappasiso <- input_tappasfiles(tappasiso_input_dir)
tappasrna <- input_tappasfiles(tappasrna_input_dir)
tappas_removediso(tappasiso$tappAS_Transcripts_InputExpressionMatrix.tsv)
tappas_removediso(tappasrna$tappAS_Transcripts_InputExpressionMatrix.tsv)
targetedtappas_isoexp <- tappas_resultsanno(targeted.class.files,tappasiso$input_normalized_matrix.tsv,tappas_phenotype)
targetedtappas_rnaexp <- tappas_resultsanno(targeted.class.files,tappasrna$input_normalized_matrix.tsv,tappas_phenotype)

# Gene Expression Plots 
# using IsoSeq FL read count as expression input
targeted_isogeneexp_plots <- lapply(lapply(unique(targetedtappas_isoexp$GeneExp$associated_gene), function(gene) plot_mergedexp(gene,"NA",targetedtappas_isoexp$GeneExp,targetedtappas_isoexp$Norm_transcounts)),ggplotGrob)
names(targeted_isogeneexp_plots) <- unique(targetedtappas_isoexp$GeneExp$associated_gene)
# using RNASeq abundance as expression input
targeted_rnageneexp_plots <- lapply(lapply(unique(targetedtappas_rnaexp$GeneExp$associated_gene), function(gene) plot_mergedexp(gene,"NA", targetedtappas_rnaexp$GeneExp, targetedtappas_rnaexp$Norm_transcounts)),ggplotGrob)
names(targeted_rnageneexp_plots) <- unique(targetedtappas_rnaexp$GeneExp$associated_gene)

pdf (paste0(output_plot_dir,"/TargetedTranscriptome.pdf"), width = 10, height = 15)
bottom_row <- plot_grid(QC_yield[[3]],QC_yield[[2]],labels = c('b', 'c'), label_size = 30, label_fontfamily = "CM Roman", ncol = 2,scale = 0.9)
plot_grid(QC_yield[[1]],bottom_row,labels = c('a', ''), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
plot_grid(ontarget,NULL,NULL, label_size = 30, label_fontfamily = "CM Roman", nrow = 3, scale = 0.9)
plot_grid(sq_num[[1]])
plot_grid(sq_num[[2]])
plot_grid(sqantifilter_plots[[2]],sqantifilter_plots[[1]],NULL,labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", nrow = 3, scale = 0.9)
plot_grid(wholevstargeted[[1]],wholevstargeted[[2]],wholevstargeted[[3]],labels = "auto", label_size = 30, label_fontfamily = "CM Roman", nrow = 3, scale = 0.9)
grobs <- ggplotGrob(targeted_filtered[[1]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
top_row <- plot_grid(targeted_filtered[[1]] + theme(legend.position = "none"),legend,labels = c('a', ''), label_size = 30, label_fontfamily = "CM Roman", ncol = 2,scale = 0.9)
plot_grid(top_row,targeted_filtered[[2]],targeted_filtered[[3]],labels = c("","b","c"), label_size = 30, label_fontfamily = "CM Roman", nrow = 3, scale = 0.9)
dev.off()


pdf (paste0(output_plot_dir,"/TargetedDifferentialAnalysis.pdf"), width = 10, height = 15)
# IsoSeq FL as expression input 
group_plots(ADReg_Genes,targeted_isogeneexp_plots)
group_plots(GWAS_Genes,targeted_isogeneexp_plots)
group_plots(FTD_Genes,targeted_isogeneexp_plots)
group_plots(EWAS_Genes,targeted_isogeneexp_plots)
# RNASeq as expression input 
group_plots(ADReg_Genes,targeted_rnageneexp_plots)
group_plots(GWAS_Genes,targeted_rnageneexp_plots)
group_plots(FTD_Genes,targeted_rnageneexp_plots)
group_plots(EWAS_Genes,targeted_rnageneexp_plots)
dev.off()

pdf (paste0(output_plot_dir,"/TargetedDifferentialAnalysis_RNAvsIso.pdf"), width = 10, height = 7)
for(gene in c(ADReg_Genes,GWAS_Genes,FTD_Genes,EWAS_Genes)){print(plot_grid(group_plots_rnavsiso(gene,targeted_isogeneexp_plots,targeted_rnageneexp_plots)))}
dev.off()



