## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##



## ---------- Source function and config files -----------------
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(SC_ROOT, "rTg4510_config.R"))
source(paste0(SC_ROOT,"bin/draw_heatmap_gene_level.R"))
source(paste0(SC_ROOT,"bin/find_mapt.R"))


# heatmap 
#draw_heatmap_gene_level(loaded$glob_iso$results_gene, annotated$glob_iso$GeneExp,"glob_isoseq",diff="yes")
#draw_heatmap_gene_level(loaded$glob_iso$results_gene, annotated$glob_iso$GeneExp,"glob_isoseq",diff="no")


## ----- MAPT transgene -----

pMapt <- list(
  glob_iso = find_mapt_isoseq(paste0(dirnames$glob_root,"/2_post_isoseq3/11_transgene"),phenotype$whole_rTg4510_iso),
  targ_iso = find_mapt_isoseq(paste0(dirnames$targ_iso_root,"/10_characterise/transgene"),phenotype$targeted_rTg4510_iso),
  targ_ont = find_mapt_ont(paste0(dirnames$targ_ont_root,"/0_characterise/transgene"),phenotype$targeted_rTg4510_ont)
)


## ----- compare effect size of gene expression (global isoseq data) ------

pGlobIsoVsRna <- list(
  genotype = density_plot(GlobalDESeq$resGeneComparison$genotype, "WaldStatistic_group_CASE_vs_CONTROL", "WaldStatistic_Genotype_TG_vs_WT", 
                          "Gene Wald statistic (Iso-Seq)", "Gene Wald statistic (RNA-Seq)", ""),
  interaction = density_plot(GlobalDESeq$resGeneComparison$interaction, "LRTStatistic.x", "LRTStatistic.y", 
                             "Gene LRT statistic (Iso-Seq)", "Gene LRT statistic (RNA-Seq)", "")
)

## ----- compare global WT and TG descriptive data -----

pGlobWTvsTG <- list(
  length = plot_iso_length_mdatasets(group_class.files, "length", "Isoform  length (kb)") +
    scale_x_discrete(limits = c("WT","TG")) +
    scale_fill_manual(values = c(label_colour("TG"), label_colour("WT"))) + 
    scale_y_continuous(labels = ks),
  
  exon = plot_iso_length_mdatasets(group_class.files, "exons", "Number of exons") +
    scale_x_discrete(limits = c("WT","TG")) +
    scale_fill_manual(values = c(label_colour("TG"), label_colour("WT"))), 
  
  numiso = group_class.files %>% group_by(Dataset, associated_gene) %>% tally() %>%
    plot_iso_length_mdatasets(., "n", "Number of isoforms per gene") +
    scale_x_discrete(limits = c("WT","TG")) +
    scale_fill_manual(values = c(label_colour("TG"), label_colour("WT")))
)

## ----- global Iso-Seq data C4b: differential gene and expression analysis  -----
pC4b <- list(
  RNAGeneExp = plot_trans_exp_individual_overtime("7022",GlobalDESeq$RresGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "RNA-Seq Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  
  IsoGeneExp = plot_trans_exp_individual_overtime("7022",GlobalDESeq$resGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "Iso-Seq Gene Expression"),
  
  tracks1 = NULL,
  
  # PB.7022.9 = red, PB.7022.8 = blue
  IsoTransExp  = plot_transexp_overtime("C4b",GlobalDESeq$resTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.7022.9")) + 
    scale_colour_manual(values = c("#7CAE00","#F8766D")) + theme(legend.position = "None") + 
    labs(title = "", subtitle = "Iso-Seq Transcript Expression"),
  
  # PB.2973.16 = red
  RNATransExp = plot_transexp_overtime("C4b",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.7022.9")) + 
    theme(legend.position = "None") + scale_y_continuous(labels = ks) + 
    labs(title = "", subtitle = "RNA-Seq Transcript Expression", y = "Normalised counts (K)"),
  
  # PB.7022.9 = red, PB.7022.13 = yellow, PB.7022.1 = black
  RNATransExpRanked = plot_transexp_overtime("C4b",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.7022.9")) +
    scale_colour_manual(values = c("black", wes_palette("Zissou1")[4],"#F8766D")) + 
    theme(legend.position = "None") + scale_y_continuous(labels = ks) +
    labs(title = "", subtitle = "RNA-Seq Transcript Expression", y = "Normalised counts (K)")
)


plot_transexp_overtime("C4b",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.7022.9"))
plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.14646.39341")) 
plot_transexp_overtime("Snca",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.38419.249")) 

pNovelOntTargetedDiff <- list(
  Clu1 = plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.14646.39352")), 
  Fyb1 = plot_transexp_overtime("Fyn",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.3948.4511")),
  Clu2 = plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.14646.35283")),
  Apoe = plot_transexp_overtime("Apoe",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.40586.1023")) 
)

plot_transexp_overtime("Mapt",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.8675.37810"))
plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.22007.224"))
plot_transexp_overtime("Snca",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.38419.87"))
plot_transexp_overtime("App",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.19309.7497"))


## ---------- Isoform Fraction -----------------

pIF <- list(
  ontNorm = lapply(Targeted$Genes, function(x) plotIF(x,
                                                      ExpInput=Exp$targ_ont$normAll,
                                                      pheno=phenotype$targeted_rTg4510_ont,
                                                      cfiles=class.files$targ_all,
                                                      design="time_series",
                                                      majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso))),
  
  isoNorm = lapply(Targeted$Genes, function(x) plotIF(x,
                                                      ExpInput=Exp$targ_iso$normAll,
                                                      pheno=phenotype$targeted_rTg4510_iso,
                                                      cfiles=class.files$targ_all,
                                                      design="time_series",
                                                      majorIso=row.names(TargetedDIU$isoDIUGeno$keptIso)))
)
for(i in 1:length(pIF)){names(pIF[[i]]) <- Targeted$Genes}



## ---------- Target Genes overview -----------------

class.files$targ_filtered %>% select(isoform, contains("ONT")) %>% mutate(annot_transcript_id = isoform)

gES <- plot_summarised_ES(Targeted$Gene_class, class.files$targ_filtered, dirnames$targ_anno, Targeted$ref_gencode, Targeted$ref_altcon)
gIR <- plot_summarised_IR(class.files$targ_filtered, dirnames$targ_anno, Targeted$ont_abundance)
gIR1 <- generate_cowplot(gIR_plots[[1]],gIR_plots[[2]],gIR_plots[[3]], num=3,nrow=1,ncol=3)
gIR2 <- generate_cowplot(gIR_plots[[4]],num="4D",nrow=1,ncol=1)

IR <- input_FICLE_splicing_results(dirnames$targ_anno,"Intron")
IR %>% filter(associated_gene == "Fus")
## ---------- Output -----------------

pdf(paste0(output_dir,"/SuppFigures.pdf"), width = 20, height = 10)
plot_grid(pMapt$glob_iso[[2]],get_legend(p_mapt$glob_iso[[1]]),p_mapt$targ_iso[[2]], p_mapt$targ_ont, labels = c("A","","B","C"))
plot_grid(pGlobIsoVsRna$genotype,pGlobIsoVsRna$interaction, labels = c("A","B"))
plot_grid(plotlist = pGlobWTvsTG , ncol = 2, nrow = 2, labels = c("A","B","C"))
plot_grid(plotlist = pC4b, ncol = 2, nrow = 3, labels = c("A","B","C","D","E","F"))
plot_grid(plotlist =  pNovelOntTargetedDiff, nrow = 1, labels = c("A","B","C","D"))
plot_grid(pIF$ontNorm$Fus[[1]],pIF$ontNorm$Bin1[[1]],pIF$ontNorm$Trpa1[[1]],pIF$ontNorm$Apoe[[1]],pIF$ontNorm$App[[1]], labels = c("A","B","C","D","E"))
plot_grid(pIF$ontNorm$Bin1[[2]],pIF$ontNorm$Vgf[[2]],pIF$ontNorm$Clu[[2]],pIF$ontNorm$Apoe[[2]],pIF$ontNorm$Fus[[2]], labels = c("A","B","C","D","E"))
dev.off()

pdf(paste0(dirnames$targ_output,"/IsoTargetedIFUpdated.pdf"), width = 14, height = 8)
for(i in Targeted$Genes){
  print(i)
  grid.arrange(grobs =  pIF$isoNorm[[i]], nrow = 1)
}
dev.off()

pdf(paste0(dirnames$targ_output,"/ONTTargetedIFUpdated.pdf"), width = 14, height = 8)
for(i in Targeted$Genes){
  print(i)
  grid.arrange(grobs =  pIF$ontNorm[[i]], nrow = 1)
}
dev.off()

pdf(paste0(dirnames$targ_output,"/ComparisonTargetedIF.pdf"), width = 14, height = 8)
for(i in Targeted$Genes){
  print(i)
  print(plot_grid(pIF$ontNorm[[i]][[2]],pIF$isoNorm[[i]][[2]],nrow = 1))
}
dev.off()