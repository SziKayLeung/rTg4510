## ----------Script-----------------
##
## Purpose: Generate plots for annotation of target genes after running FICLE
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------------------------------


## ---------- Source function and config files -----------------

# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/SQANTI_class_preparation.R"))
sapply(list.files(path = paste0(LOGEN,"aesthetics_basics_plots"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"target_gene_annotation"), pattern="*.R", full = T), source,.GlobalEnv)

# project related scripts and functions
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/1_ONT_Pipeline/"
source(paste0(SC_ROOT, "02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "rTg4510_ont_characterise.config.R"))

# output directory
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"


## ---------- Target Genes overview -----------------

Merged_gene_class_df <- all_summarise_gene_stats(misc_input$Gene_class, input.class.files$merged_noISM,misc_input$cpat,misc_input$noORF,misc_input$TargetGene)

gAS <- plot_summarised_AS_events(Merged_gene_class_df, misc_input$target_anno_dir, misc_input$ref_gencode)
gES <- plot_summarised_ES(misc_input$Gene_class, input.class.files$merged_noISM, misc_input$target_anno_dir, misc_input$ref_gencode, misc_input$ref_altcon)
gIR <- plot_summarised_IR(input.class.files$merged_noISM, misc_input$target_anno_dir, misc_input$ONT_abundance)
gIR1 <- generate_cowplot(gIR_plots[[1]],gIR_plots[[2]],gIR_plots[[3]], num=3,nrow=1,ncol=3)
gIR2 <- generate_cowplot(gIR_plots[[4]],num="4D",nrow=1,ncol=1)


## ---------- Annotation of each target gene -----------------

pES <- lapply(misc_input$TargetGene, function(x) plot_ES_Tgene(misc_input$target_anno_dir, x, input.class.files$merged_noISM))
pIR <- lapply(misc_input$TargetGene, function(x) plot_IR_Tgene(misc_input$target_anno_dir, x, input.class.files$merged_noISM))
pA5A3 <- lapply(misc_input$TargetGene, function(x) plot_A5A3_Tgene(misc_input$target_anno_dir, x))
pDendro <- lapply(misc_input$TargetGene, function(x) plot_dendro_Tgene(misc_input$target_anno_dir, x))
pFirst <- lapply(misc_input$TargetGene, function(x) plot_FirstExon_Tgene(misc_input$target_anno_dir, x))

names(pES) <- misc_input$TargetGene
names(pIR) <- misc_input$TargetGene
names(pA5A3) <- misc_input$TargetGene
names(pDendro) <- misc_input$TargetGene
names(pFirst) <- misc_input$TargetGene


## ---------- Output -----------------

pdf(paste0(output_dir,"/TargetGenes.pdf"), width = 10, height = 15)
# Abca1
generate_cowplot(pDendro$Abca1,NULL,num=1,nrow=2,ncol=1)
# Abca7 
plot_grid(pDendro$Abca7,pES$Abca7[[1]],pIR$Abca7[[2]],NULL,NULL,NULL,rel_widths = c(0.45,0.1,0.45), ncol = 3, nrow = 2)
# Ank1
plot_grid(pDendro$Ank1, pES$Ank1[[1]],nrow=3,ncol=3, rel_widths = c(0.6,0.2,0.2))
# Apoe
generate_cowplot(pDendro$Apoe,NULL,pA5A3$Apoe,NULL,pES$Apoe[[1]],NULL,num="3rw",nrow=3,ncol=2)
# App
generate_cowplot(pDendro$App,NULL,pES$App[[1]],NULL,pES$App[[2]],NULL,num="3rw",nrow=4,ncol=2)
# Bin1
generate_cowplot(pDendro$Bin1,NULL,pES$Bin1[[1]],NULL,pIR$Bin1[[2]],NULL,num="3rw",nrow=4,ncol=2)
# Cd33
plot_grid(pDendro$Cd33,pIR$Cd33[[1]],pIR$Cd33[[2]],NULL,NULL,NULL,rel_widths = c(0.4,0.3,0.3), ncol = 3, nrow = 2)
# Clu
generate_cowplot(pDendro$Clu,NULL,pES$Clu[[1]],NULL,pA5A3$Clu,NULL,pIR$Clu[[2]],num="3rw",nrow=4,ncol=2)
# Fus
generate_cowplot(pDendro$Fus,NULL,pES$Fus[[1]],NULL,pIR$Fus[[2]],num="3rw",nrow=3,ncol=2)
# Fyn
generate_cowplot(pDendro$Fyn,pES$Fyn[[1]],NULL,NULL,NULL,NULL,num="3rw",nrow=3,ncol=2)
# Mapt
generate_cowplot(pDendro$Mapt,NULL, pES$Mapt[[1]],NULL, pES$Mapt[[2]],NULL,pIR$Mapt[[2]],num="3rw",nrow=4,ncol=2)
# Picalm
generate_cowplot(pDendro$Picalm,NULL,pFirst$Picalm,NULL,pES$Picalm[[1]],num="3rw",nrow=3,ncol=2)
# Ptk2b
generate_cowplot(pDendro$Ptk2b,NULL,pES$Ptk2b[[1]],NULL,pA5A3$Ptk2b,num="3rw",nrow=3,ncol=2)
# Rhbdf2
generate_cowplot(pDendro$Rhbdf2,NULL,num=1,nrow=2,ncol=1)
# Snca
generate_cowplot(pDendro$Snca,NULL,pES$Snca[[1]],NULL,pA5A3$Snca,num="3rw",nrow=3,ncol=2)
# Sorl1
plot_grid(pDendro$Sorl1,pA5A3$Sorl1,nrow=3,ncol=2, rel_widths = c(0.6,0.4))
# Tardbp
generate_cowplot(pDendro$Tardbp,NULL,pIR$Tardbp[[2]],NULL,pA5A3$Tardbp,NULL,num="3rw",nrow=3,ncol=2)
# Trpa1
generate_cowplot(pDendro$Trpa1,NULL,num=1,nrow=2,ncol=1)
# Trem2
plot_grid(pDendro$Trem2,NULL,pA5A3$Trem2,NULL,NULL,NULL,rel_widths = c(0.4,0.6), ncol = 2, nrow = 2)
# Vgf
plot_grid(pDendro$Vgf,NULL,pA5A3$Vgf,nrow=2,ncol=2, rel_widths = c(0.6,0.4))
dev.off()