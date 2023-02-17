## ----------Script-----------------
##
## Purpose: Generate plots for annotation of target genes after running FICLE
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------------------------------

suppressMessages(library("cowplot"))

## ---------- Source function and config files -----------------

# source all general scripts related to long-read sequencing
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))
sapply(list.files(path = paste0(LOGEN,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"aesthetics_basics_plots"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)

# project related scripts and functions
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/"
#source(paste0(SC_ROOT, "/1_ONT_Pipeline/02_source_characterise_functions.R"))
source(paste0(SC_ROOT, "/1_ONT_Pipeline/rTg4510_ont_characterise.config.R"))
source(paste0(SC_ROOT, "/2_Merged_Isoform_Characterisation/B_cupcake_pipeline/rTg4510_merged_characterise.config.R"))

# filtering by isoforms
targeted.class.files$cupcake_merged <- targeted.class.files$cupcake_merged %>% filter(filter_result == "Isoform") 

p_threshold <- list(
  both = plot_cupcake_collapse_sensitivity(subset(targeted.class.files$cupcake_merged, dataset == "Both"),"All 20 target genes, both datasets"),
  iso = plot_cupcake_collapse_sensitivity(subset(targeted.class.files$cupcake_merged, dataset == "Iso-Seq"),"All 20 target genes, Iso-Seq"),
  ont = plot_cupcake_collapse_sensitivity(subset(targeted.class.files$cupcake_merged, dataset == "ONT"),"All 20 target genes, ONT")
)

plot_grid(p_threshold$both[[1]],p_threshold$ont[[1]],p_threshold$iso[[1]])
plot_grid(p_threshold$both[[2]],p_threshold$ont[[2]],p_threshold$iso[[2]])
total_num_iso(targeted.class.files$cupcake_merged,"All isoforms (prefiltered)","dataset")
final.class.files <- filter_class_by_counts(targeted.class.files$cupcake_merged,nread_threshold=10,nsample_threshold=5)
total_num_iso(final.class.files,"Filtered isoforms","category")
total_num_iso(final.class.files,"Filtered isoforms","dataset")


## ---------- Target Genes overview -----------------

Merged_gene_class_df <- all_summarise_gene_stats(misc_input$Gene_class, input.class.files$cupcake_merged, misc_input$cpat ,misc_input$noORF, misc_input$TargetGene)

gAS <- plot_summarised_AS_events(Merged_gene_class_df, dirnames$targetgenes, misc_input$ref_gencode)
gES <- plot_summarised_ES(misc_input$Gene_class, input.class.files$merged_noISM, dirnames$targetgenes, misc_input$ref_gencode, misc_input$ref_altcon)
gIR <- plot_summarised_IR(input.class.files$merged_noISM, dirnames$targetgenes, misc_input$ONT_abundance)
gIR1 <- generate_cowplot(gIR_plots[[1]],gIR_plots[[2]],gIR_plots[[3]], num=3,nrow=1,ncol=3)
gIR2 <- generate_cowplot(gIR_plots[[4]],num="4D",nrow=1,ncol=1)


## ---------- Annotation of each target gene -----------------

pES <- lapply(misc_input$TargetGene, function(x) plot_ES_Tgene(dirnames$targetgenes, x, input.class.files$merged_noISM))
pIR <- lapply(misc_input$TargetGene, function(x) plot_IR_Tgene(dirnames$targetgenes, x, input.class.files$merged_noISM))
pA5A3 <- lapply(misc_input$TargetGene, function(x) plot_A5A3_Tgene(dirnames$targetgenes, x))
pDendro <- lapply(misc_input$TargetGene, function(x) plot_dendro_Tgene(dirnames$targetgenes, x))
pFirst <- lapply(misc_input$TargetGene, function(x) plot_FirstExon_Tgene(dirnames$targetgenes, x))

names(pES) <- misc_input$TargetGene
names(pIR) <- misc_input$TargetGene
names(pA5A3) <- misc_input$TargetGene
names(pDendro) <- misc_input$TargetGene
names(pFirst) <- misc_input$TargetGene

