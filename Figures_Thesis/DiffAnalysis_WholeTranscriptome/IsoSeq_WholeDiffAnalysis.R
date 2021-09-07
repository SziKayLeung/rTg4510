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
output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_WholeTranscriptome"
# results from Whole transcriptome paper
input_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_WholeTranscriptome/IsoSeq_WholeDiffAnalysis_Variables.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_WholeTranscriptome/IsoSeq_WholeDiffAnalysis_Functions.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_WholeTranscriptome/AS_Functions.R")

##### Load Files ############################# 
tappasiso <- input_tappasfiles(tappasiso_input_dir)
tappasrna <- input_tappasfiles(tappasrna_input_dir)

# tappasrna expression file: change the column to include phenotype and age
# tappasrna_phenotype = phenotype file used for input to tappas with all the mouse for rnaseq expresssion input
for(c in 1:length(colnames(tappasrna$input_normalized_matrix.tsv))){
  colnames(tappasrna$input_normalized_matrix.tsv)[c] <- tappasrna_phenotype[colnames(tappasrna$input_normalized_matrix.tsv)[c] == tappasrna_phenotype$sample,"variable"]
}

#tappas_removediso(tappasiso$tappAS_Transcripts_InputExpressionMatrix.tsv)
#tappas_removediso(tappasrna$tappAS_Transcripts_InputExpressionMatrix.tsv)
wholetappas_isoexp <- tappas_resultsanno(class.files,tappasiso$input_normalized_matrix.tsv,tappasiso_phenotype)
wholetappas_rnaexp <- tappas_resultsanno(class.files,tappasrna$input_normalized_matrix.tsv,tappasrna_phenotype)

# Not all target genes detected in Whole transcriptome
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")
for(gene in TargetGene){if(!gene %in% wholetappas_isoexp$Norm_transcounts$associated_gene ){print(gene)}}
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Rhbdf2","Trem2","Vgf")

#####  Different models from Tappas output ############################# 
# segregate differential gene and transcript expression results by the different models (using beta coefficient as filters)
# IsoSeq as Expression Input
gene_sigs_WholeIso_lst = segregate_tappasresults(tappassiggene$WholeIso_Genexp,"IsoSeq")
trans_sigs_WholeIso_lst = segregate_tappasresults(tappassigtrans$WholeIso_Transexp,"IsoSeq")

gene_sigs_WholeRNA_lst = segregate_tappasresults(tappassiggene$WholeRNA_Genexp,"RNASeq")
trans_sig_WholeRNA_lst = segregate_tappasresults(tappassigtrans$WholeRNA_Transexp,"RNASeq")
tappasgenesig_plots <- tappas_genesig()

# Find examples of the different models using the 3rd gene from each list 
# 3rd gene rather than top 2 for more variation in thesis
genesigs_model = sapply(gene_sigs_WholeIso_lst$models, function(x) x[3,1]) %>% reshape2::melt(value.name = "Gene")
genesigs_model_plots = lapply(lapply(1:length(genesigs_model$Gene), function(i) plot_mergedexp(genesigs_model$Gene[[i]],"NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,paste0("Model ",i))),ggplotGrob)


##### Differential Gene Expression ############################# 
### Gene Expression Plots of top significant genes
# siggenes = top3 from IsoSeq+IsoSeq except for Mapt
siggenes = c("Gfap","C4b","Tgfbr1","Slc14a1","Unc93b1","Mapt")
meanexp_output <- list()   # for output results

# using IsoSeq FL read count as expression input
whole_isogeneexp_plots <- lapply(lapply(siggenes, function(gene) plot_mergedexp(gene,"NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")),ggplotGrob)
names(whole_isogeneexp_plots) <- siggenes
#meanexp_output <- do.call("rbind",meanexp_output) 

# using RNASeq abundance as expression input
whole_rnageneexp_plots <- lapply(lapply(siggenes, function(gene) plot_mergedexp(gene,"NA", wholetappas_rnaexp$GeneExp, wholetappas_rnaexp$Norm_transcounts, "RNA-Seq Expression")),ggplotGrob)
names(whole_rnageneexp_plots) <- siggenes


### Gene Expression Plots of novelGenes
## common novel genes identified with significant gene expression changes, using IsoSeq and RNA-Seq as expression
signovelGenes = intersect(c(tappassiggene$WholeIso_Genexp[grepl("novelGene",tappassiggene$WholeIso_Genexp$...1),1])$`...1`,
                          c(tappassiggene$WholeRNA_Genexp[grepl("novelGene",tappassiggene$WholeRNA_Genexp$...1),1])$`...1`)

# using IsoSeq as expression input
whole_isonovelgeneexp_plots <- lapply(lapply(signovelGenes, function(gene) plot_mergedexp(gene,"NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")),ggplotGrob)
names(whole_isonovelgeneexp_plots ) <- signovelGenes

# using RNASeq as expression input
whole_rnanovelgeneexp_plots <- lapply(lapply(siggenes, function(gene) plot_mergedexp(gene,"NA", wholetappas_rnaexp$GeneExp, wholetappas_rnaexp$Norm_transcounts, "RNA-Seq Expression")),ggplotGrob)

## Antisense mechanism of novel genes? 
## novelgenes were identified as differentially expressed using RNASeq as expression input
#gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`[grepl("novelGene",gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1),]
# tappassiggene$WholeRNA_Genexp[grepl("novelGene",tappassiggene$WholeRNA_Genexp$...1),] %>% filter(`R-squared` > 0.5)
novelGenesAS = c("novelGene_529","novelGene_Fgfr1op_AS","Fgfr1op","novelGene_Htra1_AS","Htra1")
whole_rnanovelgeneASexp_plots <- lapply(lapply(novelGenesAS, function(gene) plot_mergedexp(gene,"NA", wholetappas_rnaexp$GeneExp, wholetappas_rnaexp$Norm_transcounts, "RNA-Seq Expression")),ggplotGrob)
names(whole_rnanovelgeneASexp_plots) <- novelGenesAS


## Target Genes 
whole_targetgene_genexp_plots <- lapply(lapply(TargetGene, function(gene) plot_mergedexp(gene,"NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")),ggplotGrob)
names(whole_targetgene_genexp_plots) <- TargetGene
whole_targetgene_rnagenexp_plots <- lapply(lapply(TargetGene, function(gene) plot_mergedexp(gene,"NA",wholetappas_rnaexp$GeneExp,wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Expression")),ggplotGrob)
names(whole_targetgene_rnagenexp_plots) <- TargetGene


##### Differential Transcript Expression ############################# 
# Transcript Expression Plots 
sigtrans = c("Gfap","C4b")
sigtrans2 = c("Cd68","Osmr","Cd34","Ubqln1","Gjb2","Adam23","Mapt","Ctsd","H2-D1","Gatm","Padi2")

# using IsoSeq FL read count as expression input
# mean expression per isoform plot  
whole_isotransexp_plots <- lapply(lapply(sigtrans, function(gene) plot_transexp(gene,wholetappas_isoexp$Norm_transcounts,"isoseq","Iso-Seq Expression")),ggplotGrob)
names(whole_isotransexp_plots) <- sigtrans
# trajectory plot 
whole_isotransexp_plots_traj <- lapply(lapply(sigtrans, function(gene) plot_transexp_overtime(gene,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")),ggplotGrob)
names(whole_isotransexp_plots_traj) <- sigtrans

whole_rnatransexp_plots_traj  <- lapply(lapply(sigtrans, function(gene) plot_transexp_overtime(gene,wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Expression")),ggplotGrob)
names(whole_rnatransexp_plots_traj) <- sigtrans

## sigtrans2
whole_isotransexp2_plots_traj <- lapply(lapply(sigtrans2, function(gene) plot_transexp_overtime(gene,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression") + theme(legend.position = "bottom",legend.direction="vertical")),ggplotGrob)
whole_rnatransexp2_plots_traj  <- lapply(lapply(sigtrans2, function(gene) plot_transexp_overtime(gene,wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Expression") + theme(legend.position = "bottom",legend.direction="vertical")),ggplotGrob)
names(whole_rnatransexp2_plots_traj) <- sigtrans2
names(whole_isotransexp2_plots_traj) <- sigtrans2

# Isoforms unique to iso-seq
DEI_genes_unique_plots <- DEI_genes_unique()

# Target Genes
whole_targetgene_isotransexp_plots_traj <- lapply(lapply(TargetGene, function(gene) plot_transexp_overtime(gene,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")),ggplotGrob)
names(whole_targetgene_isotransexp_plots_traj) <- TargetGene

whole_targetgene_rnatransexp_plots_traj <- lapply(lapply(TargetGene, function(gene) plot_transexp_overtime(gene,wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Expression")),ggplotGrob)
names(whole_targetgene_rnatransexp_plots_traj) <- TargetGene


##### Differential Transcript Usage ############################# 
venndiu = DIU_analysis_output_venn()
IFiso_Esyt2 = IF_plot("Esyt2",tappasiso$gene_transcripts.tsv, tappasiso$input_normalized_matrix.tsv, "isoseq")
IFrna_Esyt2 = IF_plot("Esyt2",tappasrna$gene_transcripts.tsv, tappasrna$input_normalized_matrix.tsv, "rnaseq")
DIU_genes_exp_plots = DIU_genes_exp()
DEG_DIU = DIU_RNASEQ_results()

#################################### Generate Plots ############
pdf(paste0(output_plot_dir,"/WholeDifferentialAnalysis.pdf"), width = 10, height = 15)
plot_grid(grobTree(tappasgenesig_plots$p1),grobTree(tappasgenesig_plots$p2), labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 1, scale = 0.9)
plot_grid(whole_isogeneexp_plots$Gfap,whole_rnageneexp_plots$Gfap,whole_isogeneexp_plots$C4b,whole_rnageneexp_plots$C4b,whole_isogeneexp_plots$Tgfbr1,whole_rnageneexp_plots$Tgfbr1, labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 2, nrow = 3, scale = 0.9)
plot_grid(whole_isogeneexp_plots$Slc14a1,whole_rnageneexp_plots$Slc14a1,whole_isogeneexp_plots$Unc93b1,whole_rnageneexp_plots$Unc93b1, NULL,NULL, labels = c("a","b","c","d"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, nrow = 3, scale = 0.9)
plot_grid(whole_isotransexp_plots_traj$Gfap,whole_isotransexp_plots$Gfap, labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", ncol = 1, nrow = 2, scale = 0.9)
plot_grid(gene_sigs_WholeIso_lst$p,NULL,scale = 0.9, nrow = 2)
plot_grid(plotlist = genesigs_model_plots, ncol = 2, labels = "auto",label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(grobTree(venndiu$v1),grobTree(venndiu$v2),grobTree(venndiu$v3),grobTree(venndiu$v4), labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 1, scale = 0.9)
plot_grid(plot_transexp("Esyt2",wholetappas_isoexp$Norm_transcounts,"isoseq","Iso-Seq Expression"),plot_transexp("Esyt2",wholetappas_rnaexp$Norm_transcounts,"rnaseq","RNA-Seq Expression"), labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 1, scale = 0.9)
plot_grid(plot_transexp("Esyt2",wholetappas_isoexp$Norm_transcounts,"isoseq","Iso-Seq Expression")  + theme(legend.position = "bottom"),
          plot_transexp("Esyt2",wholetappas_rnaexp$Norm_transcounts,"rnaseq","RNA-Seq Expression")  + theme(legend.position = "bottom")
          ,IFiso_Esyt2[[3]],IFrna_Esyt2[[3]], labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
plot_grid(DIU_genes_exp_plots[[2]],DIU_genes_exp_plots[[3]],NULL,NULL,labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
plot_grid(whole_rnanovelgeneASexp_plots$novelGene_529,NULL,whole_rnanovelgeneASexp_plots$novelGene_Fgfr1op_AS,whole_rnanovelgeneASexp_plots$Fgfr1op,whole_rnanovelgeneASexp_plots$novelGene_Htra1_AS,whole_rnanovelgeneASexp_plots$Htra1, ncol = 2,  labels = c("a","","b","c","d","e"),label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(whole_isotransexp_plots_traj$C4b,whole_isotransexp_plots$C4b, labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", ncol = 1, nrow = 2, scale = 0.9)
plot_grid(whole_rnatransexp_plots_traj$Gfap,whole_rnatransexp_plots_traj$C4b,labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", ncol = 1, nrow = 2, scale = 0.9)
plot_grid(whole_isotransexp2_plots_traj$Osmr,whole_rnatransexp2_plots_traj$Osmr,whole_isotransexp2_plots_traj$Cd68,whole_rnatransexp2_plots_traj$Cd68,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
plot_grid(whole_isotransexp2_plots_traj$Cd34,whole_rnatransexp2_plots_traj$Cd34,whole_isotransexp2_plots_traj$Ubqln1,whole_rnatransexp2_plots_traj$Ubqln1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
plot_grid(DEI_genes_unique_plots$p1_density,NULL,NULL,NULL,NULL,NULL, ncol = 2, scale = 0.9)
plot_grid(DEI_genes_unique_plots$p2_slc1a3_isoseq,DEI_genes_unique_plots$p3_slc1a3_rnaseq,DEI_genes_unique_plots$p4_Gja1_isoseq,DEI_genes_unique_plots$p5_Gja1_rnaseq,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9 )
plot_grid(whole_isotransexp2_plots_traj$Ctsd,whole_isotransexp2_plots_traj$`H2-D1`,whole_isotransexp2_plots_traj$Gatm,whole_isotransexp2_plots_traj$Padi2,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
plot_grid(whole_rnatransexp2_plots_traj$Ctsd,whole_rnatransexp2_plots_traj$`H2-D1`,whole_rnatransexp2_plots_traj$Gatm,whole_rnatransexp2_plots_traj$Padi2, labels = "auto", label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
dev.off()

pdf(paste0(output_plot_dir,"/DIU_DEG_major.pdf"), width = 10, height = 15)
for(g in 1:10){print(plot_grid(plotlist = diff_across_rnaseq(DEG_DIU$DIU_DEG_maj$gene[[g]]), ncol = 2, labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))}
dev.off()

pdf(paste0(output_plot_dir,"/DIU_DEG_nomajor.pdf"), width = 10, height = 15)
for(g in 1:10){print(plot_grid(plotlist = diff_across_rnaseq(DEG_DIU$DIU_DEG_nomaj$gene[[g]]), ncol = 2, labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))}
dev.off()

pdf(paste0(output_plot_dir,"/DIU_notDEG_nomajor.pdf"), width = 10, height = 15)
for(g in 1:10){print(plot_grid(plotlist = diff_across_rnaseq(DEG_DIU$DIU_notDEG_nomaj$gene[[g]]), ncol = 2, labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))}
dev.off()

pdf(paste0(output_plot_dir,"/DIU_notDEG_major.pdf"), width = 10, height = 15)
for(g in 1:10){print(plot_grid(plotlist = diff_across_rnaseq(DEG_DIU$DIU_notDEG_maj$gene[[g]]), ncol = 2, labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))}
dev.off()


pdf("commonDIU2.pdf", width = 15, height = 10)
commonDIUisorna_prop = intersect(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_rnaseq_prop$gene)
commonDIUisorna_fc = intersect(tappasDIU$DIU_isoseq_fc$gene,tappasDIU$DIU_rnaseq_fc$gene)
commonDIUisorna = intersect(commonDIUisorna_fc,commonDIUisorna_prop)
for(gene in commonDIUisorna){
  print(plot_grid(plot_transexp(gene,wholetappas_isoexp$Norm_transcounts,"isoseq","Iso-Seq Expression"),
                  plot_transexp(gene,wholetappas_rnaexp$Norm_transcounts,"rnaseq","rnaseq Expression")))
  
  IFiso = IF_plot(gene,tappasiso$gene_transcripts.tsv, tappasiso$input_normalized_matrix.tsv, "isoseq")
  IFrna = IF_plot(gene,tappasrna$gene_transcripts.tsv, tappasrna$input_normalized_matrix.tsv, "rnaseq")
  
  print(plot_grid(IFiso[[1]],IFrna[[1]]))
  print(plot_grid(IFiso[[2]],IFrna[[2]]))
  print(plot_grid(IFiso[[3]],IFrna[[3]]))
}
dev.off()

plot_transexp_overtime("Ctse",wholetappas_isoexp$Norm_transcounts,"isoseq")


##### Alternative Splicing Events ############################# 
dASevents_files = read_files_differentialASevents(group_sqanti_dir,suppa_dir)
AS_events_diff_plots = AS_events_diff(dASevents_files$group.mono.class.files, dASevents_files$annotated_sqanti_gtf, dASevents_files$suppa2.output)
IR_ORF_plots = IR_ORF()

grobs <- ggplotGrob(AS_events_diff_plots[[1]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
pdf(paste0(output_plot_dir,"/WholeASDiff.pdf"), width = 10, height = 15)
plot_grid(AS_events_diff_plots[[1]] + theme(legend.position = "none"),AS_events_diff_plots[[2]],AS_events_diff_plots[[3]],legend,labels = c("a","b","c"), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
plot_grid(IR_ORF_plots[[1]],IR_ORF_plots[[2]],labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
plot_grid(IR_ORF_plots[[3]],IR_ORF_plots[[4]],labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
dev.off()

##### Summary Dataset ############################# 
summaryoutput = summary_dataset()


### DEI but not DEG --> differentially expressed isoform but not differentially expressed genes
'%!in%' <- function(x,y)!('%in%'(x,y))
tappassig$WholeIso_Genexp <- tappassig$WholeIso_Genexp %>% rownames_to_column(var = "associated_gene")
DEI_notDEG <- tappassigtrans$WholeIso_Transexp[tappassigtrans$WholeIso_Transexp$associated_gene %!in% tappassig$WholeIso_Genexp$associated_gene,]

# Gpam # isoform switching
dei_ndeg <- c("Gpam","Cacnb3","Gpm6b","Flot2","Qsox2")
for(gene in dei_ndeg){
  print(plot_mergedexp(gene,"NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression"))
  print(plot_transexp_overtime(gene,wholetappas_isoexp$Norm_transcounts,""))
  print(plot_transexp(gene,wholetappas_isoexp$Norm_transcounts,"isoseq"))
}

#### Target Genes 
pdf(paste0(output_plot_dir,"/WholeDifferentialAnalysis_TargetGenes.pdf"), width = 10, height = 15)
for(gene in TargetGene){
  print(plot_grid(whole_targetgene_genexp_plots[[gene]],whole_targetgene_rnagenexp_plots[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
  print(plot_grid(whole_targetgene_isotransexp_plots_traj[[gene]],whole_targetgene_rnatransexp_plots_traj[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
}
dev.off()




