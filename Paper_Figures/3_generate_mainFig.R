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
source(paste0(SC_ROOT, "rTg4510_config.R"))
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/plot_basic_stats.R"))
source(paste0(LOGEN_ROOT,"longread_QC/plot_cupcake_collapse_sensitivty.R"))
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Mouse_Isoseq/"


## ---------- Figure 1: Gfap ----------

Gfap_p <- list(
  RNAGeneExp = plot_trans_exp_individual_overtime("2973",GlobalDESeq$RresGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "RNA-Seq Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  
  IsoGeneExp = plot_trans_exp_individual_overtime("2973",GlobalDESeq$resGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "Iso-Seq Gene Expression"),
  
  tracks1 = NULL,
  
  # PB.2973.16 = red, PB.2973.28 = green, PB.2973.35 = blue
  IsoTransExp  = plot_transexp_overtime("Gfap",GlobalDESeq$resTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.2973.16")) + 
    scale_colour_manual(values = c("#F8766D","#7CAE00","#00BFC4")) + theme(legend.position = "None") + 
    labs(title = "", subtitle = "Iso-Seq Transcript Expression"),
  
  # PB.2973.16 = red
  RNATransExp = plot_transexp_overtime("Gfap",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.2973.16")) + 
    theme(legend.position = "None") + 
    labs(title = "", subtitle = "RNA-Seq Transcript Expression"),
  
  # PB.2973.16 = red, PB.2973.23 = grey, PB.2973.33 = black, PB.2973.27 = yellow
  RNATransExpRanked  = plot_transexp_overtime("Gfap",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.2973.16")) +
    scale_colour_manual(values = c("#F8766D","gray",wes_palette("Zissou1")[4],"black")) + 
    theme(legend.position = "None") + scale_y_continuous(labels = ks) +
    labs(title = "", subtitle = "RNA-Seq Transcript Expression", y = "Normalised counts (K)")
)


## ---------- Figure 2: Targeted ----------

Targeted_p <- list(
  comp = NULL,
  venn = grobTree(vennONTvsIso(class.files$targ_filtered)),
  cumulative = pSensitivity(class.files$targ_all),
  num = total_num_iso(class.files$targ_filtered %>% mutate(dataset=Dataset),""),
  splicing = plot_summarised_AS_events(Merged_gene_class_df, dirnames$targ_anno, Targeted$ref_gencode)[[1]],
  distribution = plotIFAll(Exp=Exp$targ_ont$normAll %>% select(-associated_gene),
                           classf=class.files$targ_all,
                           pheno=phenotype$targeted_rTg4510_ont,
                           majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso)) + theme(legend.position = "None")
)


## ---------- Figure 3: Trem2 ----------

Trem2_p <- list(
  pheat = draw_heatmap_gene(gene="Trem2",cf=class.files$targ_all,normCounts=TargetedDESeq$ontResTranAnno$wald$norm_counts,type="targeted"),
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Trem2"),
  ONTTransExp = plot_transexp_overtime("Trem2",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=3,
                                       isoSpecific=c("PB.20818.54"),setorder=c("CONTROL","CASE")) +
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  ONTGeneExp = plot_trans_exp_individual_overtime("20818",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("Trem2",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL,isoSpecific=c("PB.20818.54"))
)
Trem2_p$IF[[1]] <- Trem2_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
Trem2_p$IF[[2]] <- Trem2_p$IF[[2]] + labs(title = "")


## ---------- Figure 4: Bin1 ---------

Bin1_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Bin1"),
  ES = plot_ES_Tgene(dirnames$targ_anno,"Bin1",class.files$targ_filtered)[[1]],
  IR = plot_IR_Tgene(dirnames$targ_anno,"Bin1",class.files$targ_filtered)[[2]],
  tracks = NULL,
  diff_tracks = NULL,
  ONTGeneExp  = plot_trans_exp_individual_overtime("22007",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  ONTTransExp = plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=2,isoSpecific=c("PB.22007.224"),
                                       setorder=c("CONTROL","CASE")) + 
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("Bin1",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",majorIso=NULL,isoSpecific=c("PB.22007.224","PB.22007.99","PB.22007.176"))
)
Bin1_p$IF[[1]] <- Bin1_p$IF[[1]] + labs(title = "") #+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
Bin1_p$IF[[2]] <- Bin1_p$IF[[2]] + labs(title = "")


## ---------- Figure 5: Clu ---------

Clu_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Clu"),
  ES = plot_ES_Tgene(dirnames$targ_anno,"Clu",class.files$targ_filtered)[[1]],
  IR = plot_IR_Tgene(dirnames$targ_anno,"Clu",class.files$targ_filtered)[[2]],
  ONTTransExp = plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=1,isoSpecific=c("PB.14646.39341"),setorder=c("CONTROL","CASE")) +
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  ONTGeneExp = plot_trans_exp_individual_overtime("14646",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("Clu",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=NULL,majorIso=NULL,isoSpecific=c("PB.14646.139","PB.14646.60837","PB.14646.39341"))
)
Clu_p$IF[[1]] <- Clu_p$IF[[1]] + labs(title = "") #+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
Clu_p$IF[[2]] <- Clu_p$IF[[2]] + labs(title = "")


## ---------- Figure 6: App ---------

App_p <- list(
  ES = plot_ES_Tgene(dirnames$targ_anno,"App",class.files$targ_filtered)[[1]],
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "App"),
  ONTTransExp = plot_transexp_overtime("App",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=3,setorder=c("CONTROL","CASE")) +
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  ONTGeneExp = plot_trans_exp_individual_overtime("19309",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("App",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL)
)
App_p$IF[[1]] <- App_p$IF[[1]] + labs(title = "") #+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
App_p$IF[[2]] <- App_p$IF[[2]] + labs(title = "")


## ---------- Figure 7: Apoe ---------

Apoe_p <- list(
  #dendro = plot_dendro_Tgene(dirnames$targ_anno, "Apoe"),
  ONTTransExp = plot_transexp_overtime("Apoe",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=3,setorder=c("CONTROL","CASE")) +
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  ONTGeneExp = plot_trans_exp_individual_overtime("40586",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("Apoe",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL)
)


## ---------- Output Pdf -----------------

pdf(paste0(output_dir,"/MainFigures1.pdf"), width = 20, height = 10)
plot_grid(plotlist = Gfap_p, ncol = 2, nrow = 3, labels = c("A","B","C","D","E","F"))
dev.off()

pdf(paste0(output_dir,"/MainFigures2.pdf"), width = 15, height = 15)
plot_grid(plot_grid(Targeted_p$comp,Targeted_p$cumulative,Targeted_p$venn,nrow=1,labels=c("A","B","C"),rel_widths = c(0.4,0.3,0.3)),
          plot_grid(Targeted_p$n,Targeted_p$distribution,nrow=1,labels=c("D","E"),rel_widths = c(0.4,0.6)),
          plot_grid(Targeted_p$splicing,nrow=1,labels=c("F")),nrow=3, rel_heights = c(0.3,0.4,0.3))
dev.off()

pdf(paste0(output_dir,"/MainFigures3.pdf"), width = 20, height = 15)
plot_grid(NULL,
          plot_grid(
            plot_grid(Trem2_p$dendro,Trem2_p$pheat$gtable,nrow=1, rel_widths = c(0.6,0.4)),
            plot_grid(Trem2_p$ONTTransExp,Trem2_p$ONTGeneExp,nrow=1, rel_widths = c(0.6,0.4)),
            plot_grid(Trem2_p$IF[[1]],Trem2_p$IF[[2]],nrow=1, rel_widths = c(0.6,0.4)), 
            ncol = 1),
          ncol = 2, rel_widths = c(0.4,0.6)
)
dev.off()

pdf(paste0(output_dir,"/MainFigures4.pdf"), width = 20, height = 15)
plot_grid(plot_grid(
  plot_grid(Bin1_p$dendro,Bin1_p$ES,Bin1_p$IR, nrow = 1),
  plot_grid(Bin1_p$ONTTransExp, Bin1_p$ONTGeneExp),
  plot_grid(Bin1_p$IF[[1]], Bin1_p$IF[[2]]),
  ncol = 1), NULL,
  ncol = 2, rel_widths = c(0.6,0.4)
)
dev.off()

pdf(paste0(output_dir,"/OtherMainFigures.pdf"), width = 15, height = 10)
plot_grid(Apoe_p$ONTTransExp + theme(legend.position = "right"),Apoe_p$ONTGeneExp,Apoe_p$IF[[1]],Apoe_p$IF[[2]])
plot_grid(App_p$ONTTransExp + theme(legend.position = "right"),App_p$ONTGeneExp,App_p$IF[[1]],App_p$IF[[2]])
plot_grid(Clu_p$dendro, Clu_p$ONTTransExp, Clu_p$ONTGeneExp, Clu_p$IF[[1]],Clu_p$IF[[2]])
dev.off()