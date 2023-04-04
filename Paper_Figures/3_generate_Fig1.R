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


# Figure 1: Gfap
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

# Figure 2: Targeted 
# venn diagram of ONT vs Iso-Seq
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



# Figure 4: Trem2
Trem2_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Trem2"),
  tracks = NULL,
  ONTGeneExp  = generate_plots("Trem2",annotated$targ_ont,"Gene_time","")$Trem2,
  ONTTransExp = generate_plots("Trem2",annotated$targ_ont,"Transcript_Iso_Trajectory"," ",tappassigtrans$targ_ont$TargetedOnt_Transexp)$Trem2,
  IF = plot_grid(generate_plots("Trem2",loaded$targ_ont,"IF_time_series",phenotype$targ_ont,"ONT_Targeted")[[1]]),
  pheat = draw_heatmap_gene("Trem2", class.files$targ_ont, annotated$targ_ont$Norm_transcounts)$gtable
)

# Figure 5: Bin1
Bin1_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Bin1"),
  tracks = NULL,
  diff_tracks = NULL,
  ONTGeneExp  = generate_plots("Bin1",annotated$targ_ont,"Gene_time","")$Bin1,
  ONTTransExp = generate_plots("Bin1",annotated$targ_ont,"Transcript_Iso_Trajectory"," ",tappassigtrans$targ_ont$TargetedOnt_Transexp)$Bin1,
  IF = plot_grid(generate_plots("Bin1",loaded$targ_ont,"IF_time_series",phenotype$targ_ont,"ONT_Targeted")[[1]])
  #pheat = draw_heatmap_gene("Bin1", class.files$targ_ont, annotated$targ_ont$Norm_transcounts)$gtable
)



pdf(paste0(output_dir,"/MainFigures2.pdf"), width = 20, height = 10)
plot_grid(plotlist = Gfap_p, ncol = 2, nrow = 3, labels = c("A","B","C","D","E","F"))
dev.off()

pdf(paste0(output_dir,"/MainFigures2.pdf"), width = 15, height = 15)
plot_grid(plot_grid(Targeted_p$comp,Targeted_p$cumulative,Targeted_p$venn,nrow=1,labels=c("A","B","C"),rel_widths = c(0.4,0.3,0.3)),
          plot_grid(Targeted_p$n,Targeted_p$distribution,nrow=1,labels=c("D","E"),rel_widths = c(0.4,0.6)),
          plot_grid(Targeted_p$splicing,nrow=1,labels=c("F")),nrow=3, rel_heights = c(0.3,0.4,0.3))
dev.off()
#plot_grid(plot_grid(Trem2_p$dendro,Trem2_p$tracks, ncol = 2, nrow = 1, labels = c("A","B")),
#          plot_grid(Trem2_p$ONTGeneExp,Trem2_p$ONTTransExp,Trem2_p$IF,Trem2_p$pheat, ncol = 4, nrow = 1, labels = c("C","D","E","F")),
#          nrow=2)
#plot_grid(plotlist = Bin1_p, ncol = 3, nrow = 2, labels = c("A","B","","C","D","E"))
dev.off()