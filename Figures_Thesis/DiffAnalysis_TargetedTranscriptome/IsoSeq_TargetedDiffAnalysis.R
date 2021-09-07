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

# do not output log files for venn diagrams
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

#output_helpfig_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures"
output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_TargetedTranscriptome"
# results from Whole transcriptome paper
input_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_TargetedTranscriptome/IsoSeq_TargetedDiffAnalysis_Variables.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_TargetedTranscriptome/IsoSeq_TargetedDiffAnalysis_Functions.R")
#source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_TargetedTranscriptome/AS_Functions.R")

TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")

##### Load Files ############################# 
tappasiso <- input_tappasfiles(tappasiso_input_dir)
tappasrna <- input_tappasfiles(tappasrna_input_dir)
tappasisocol <- input_tappasfiles(tappasisocol_input_dir)
tappasrnacol <- input_tappasfiles(tappasrnacol_input_dir)
tappasrnacolmerged <- input_tappasfiles(tappasrnamerged_input_dir)

# tappasrna expression file: change the column to include phenotype and age
# tappasrna_phenotype = phenotype file used for input to tappas with all the mouse for rnaseq expresssion input
#for(c in 1:length(colnames(tappasrna$input_normalized_matrix.tsv))){
#  colnames(tappasrna$input_normalized_matrix.tsv)[c] <- tappasrna_phenotype[colnames(tappasrna$input_normalized_matrix.tsv)[c] == tappasrna_phenotype$sample,"variable"]
#}
#colnames(tappasrna$input_normalized_matrix.tsv) <- lapply(colnames(tappasrna$input_normalized_matrix.tsv),
#                                                          function(x) paste0(word(word(x,c(1),sep = fixed("_")),c(2),sep = fixed(".")),".",
#                                                                             word(word(x,c(1),sep = fixed("_")),c(3),sep = fixed("."))))

#tappas_removediso(tappasiso$tappAS_Transcripts_InputExpressionMatrix.tsv)
#tappas_removediso(tappasrna$tappAS_Transcripts_InputExpressionMatrix.tsv)
targetedtappas_isoexp <- tappas_resultsanno(class.files,tappasiso$input_normalized_matrix.tsv,tappasiso_phenotype)
targetedtappas_rnaexp <- tappas_resultsanno(class.files,tappasrna$input_normalized_matrix.tsv,tappasrna_phenotype)
targetedtappas_isocolexp <- tappas_resultsanno(coll.class.files,tappasisocol$input_normalized_matrix.tsv,tappasiso_phenotype)
targetedtappas_rnacolexp <- tappas_resultsanno(coll.class.files,tappasrnacol$input_normalized_matrix.tsv,tappasrna_phenotype)
targetedtappas_rnacolmergedexp <- tappas_resultsanno(merged.class.files,tappasrnacolmerged$input_normalized_matrix.tsv,tappasrna_phenotype)

targetedtappas_rnacolmergedexp$Norm_transcounts <- reannotate_tamamerged_output(targetedtappas_rnacolmergedexp$Norm_transcounts)
#### Annotate TappasMerged files from RNA-Seq expression to original PB.ID 


#####  Different models from Tappas output ############################# 
# segregate differential gene and transcript expression results by the different models (using beta coefficient as filters)
# IsoSeq as Expression Input
gene_sigs_TargetedIso_lst = segregate_tappasresults(tappassiggene$TargetedIso_Genexp,"IsoSeq")
trans_sigs_TargetedIso_lst = segregate_tappasresults(tappassigtrans$TargetedIso_Transexp,"IsoSeq")

gene_sigs_TargetedRNA_lst = segregate_tappasresults(tappassiggene$TargetedRNA_Genexp,"RNASeq")
trans_sig_TargetedRNA_lst = segregate_tappasresults(tappassigtrans$TargetedRNA_Transexp,"RNASeq")
#tappasgenesig_plots <- tappas_genesig()


##### Differential Gene Expression #############################
meanexp_output <- list()   # for output results
generate_plots <- function(genelist,tappasinput,type,name){
  if(type == "Gene"){
    plist <- lapply(lapply(genelist, function(gene) plot_mergedexp(gene,"NA",tappasinput[["GeneExp"]],tappasinput[["Norm_transcounts"]],name)),ggplotGrob)
  }else if(type == "Transcript_Isoseq"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp(gene,tappasinput[["Norm_transcounts"]],"isoseq",name)),ggplotGrob)
  }else if(type == "Transcript_Rnaseq_Targeted"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp(gene,tappasinput[["Norm_transcounts"]],"rnaseq",name)),ggplotGrob)
  }else if(type == "Transcript Trajectory"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],name)),ggplotGrob)
  }else{
    print("Type Required")
  }
  
  names(plist) = genelist
  return(plist)
}

# using IsoSeq FL read count as expression input
targeted_isogeneexp_plots <- generate_plots(TargetGene,targetedtappas_isoexp,"Gene","Iso-Seq Expression")

# using RNASeq abundance as expression input
targeted_rnageneexp_plots <- generate_plots(TargetGene,targetedtappas_rnaexp,"Gene","RNA-Seq Expression")

# using IsoSeq FL read count as expression input collapsed by partial transcripts 
targeted_isocolgeneexpcol_plots <- generate_plots(TargetGene,targetedtappas_isocolexp,"Gene","Iso-Seq Expression")

# using RNA-Seq as expression input, mapped on Iso-Seq transcripts collapsed by partial transcripts 
targeted_rnacolgeneexpcol_plots <- generate_plots(TargetGene,targetedtappas_rnacolexp,"Gene","RNA-Seq Expression")

# using RNA-Seq as expression input, mapped on Targeted and Whole Transcriptome merged 
targeted_rnacolmergedgeneexpcol_plots <- generate_plots(TargetGene,targetedtappas_rnacolmergedexp,"Gene","RNA-Seq Expression")


##### Differential Transcript Expression ############################# 
# using IsoSeq FL read count as expression input
targeted_isotransexp_plots <- generate_plots(TargetGene,targetedtappas_isoexp,"Transcript_Isoseq","Iso-Seq Expression")
  
# trajectory plot 
targeted_isotransexp_plots_traj <- generate_plots(TargetGene,targetedtappas_isoexp,"Transcript Trajectory","Iso-Seq Expression")
targeted_rnatransexp_plots_traj <-  generate_plots(TargetGene,targetedtappas_rnaexp,"Transcript Trajectory","RNA-Seq Expression")
targeted_isocoltransexp_plots_traj <- generate_plots(TargetGene,targetedtappas_isocolexp,"Transcript Trajectory","Iso-Seq Expression")
targeted_rnacoltransexp_plots_traj <- generate_plots(TargetGene,targetedtappas_rnacolexp,"Transcript Trajectory","RNA-Seq Expression")
targeted_rnacolmergedtransexp_plots_traj <- generate_plots(TargetGene,targetedtappas_rnacolmergedexp,"Transcript Trajectory","RNA-Seq Expression")

# Non trajectory plot 
targeted_rnacolmergedtransexp_plots <- generate_plots(TargetGene,targetedtappas_rnacolmergedexp,"Transcript_Rnaseq_Targeted","RNA-Seq Expression")
#targeted_rnacolmergedtransexp_plots$App <- targeted_rnacolmergedtransexp_plots$App %>% gg.gap::gg.gap(ylim = c(0, 20000), segments = list(c(5000, 9000)), tick_width = 2000, c(0.7,0,0.3)) 

targeted_isocolmergedtransexp_plots <- generate_plots(TargetGene,targetedtappas_isocolexp,"Transcript_Isoseq","Iso-Seq Expression")

# Isoforms unique to iso-seq
#DEI_genes_unique_plots <- DEI_genes_unique()

##### Differential Transcript Usage ############################# 
#venndiu = DIU_analysis_output_venn()
#IFiso_Esyt2 = IF_plot("Esyt2",tappasiso$gene_transcripts.tsv, tappasiso$input_normalized_matrix.tsv, "isoseq")
#IFrna_Esyt2 = IF_plot("Esyt2",tappasrna$gene_transcripts.tsv, tappasrna$input_normalized_matrix.tsv, "rnaseq")
#DIU_genes_exp_plots = DIU_genes_exp()
#DEG_DIU = DIU_RNASEQ_results()

#################################### Generate Plots ############
pdf(paste0(output_plot_dir,"/TargetedDifferentialAnalysis.pdf"), width = 10, height = 15)
for(gene in TargetGene){
  print(plot_grid(targeted_isogeneexp_plots[[gene]],targeted_rnageneexp_plots[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
  print(plot_grid(targeted_isotransexp_plots_traj[[gene]],targeted_rnatransexp_plots_traj[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
}
dev.off()

pdf(paste0(output_plot_dir,"/TargetedDifferentialAnalysis_Collapsed.pdf"), width = 10, height = 15)
for(gene in TargetGene){
  print(plot_grid(targeted_isocolgeneexpcol_plots[[gene]],targeted_rnacolgeneexpcol_plots[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9,rel_heights = c(1,2)))
  print(plot_grid(targeted_isocoltransexp_plots_traj[[gene]],targeted_rnacoltransexp_plots_traj[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
  }
dev.off()

pdf(paste0(output_plot_dir,"/TargetedDifferentialAnalysis_CollapsedMerged.pdf"), width = 10, height = 15)
for(gene in TargetGene){
  print(plot_grid(targeted_isocolgeneexpcol_plots[[gene]],targeted_rnacolmergedgeneexpcol_plots[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9,rel_heights = c(1,2)))
  print(plot_grid(targeted_isocoltransexp_plots_traj[[gene]],targeted_rnacolmergedtransexp_plots_traj[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
  print(plot_grid(targeted_isocolmergedtransexp_plots[[gene]],targeted_rnacolmergedtransexp_plots[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
}
dev.off()




