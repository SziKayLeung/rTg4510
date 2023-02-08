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

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "rTg4510_config.R"))
source(paste0(SC_ROOT, "0_source_functions.R"))
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Mouse_Isoseq/"


# Figure 1: Gfap
Gfap_p <- list(
  RNAGeneExp = generate_plots("Gfap",annotated$glob_rna,"Gene_time","")$Gfap,
  IsoGeneExp = generate_plots("Gfap",annotated$glob_iso,"Gene_time","")$Gfap,
  tracks1 = NULL,
  tracks2 = NULL,
  RNATransExp = generate_plots("Gfap",annotated$glob_rna,"Transcript_Rna_Trajectory"," ",tappassigtrans$glob$WholeIso_Transexp)$Gfap,
  IsoTransExp = generate_plots("Gfap",annotated$glob_iso,"Transcript_Iso_Trajectory"," ",tappassigtrans$glob$WholeIso_Transexp)$Gfap
)


# Figure 4: Trem2
Trem2_p <- list(
  dendro = dendro_plot(dirnames$targ_anno, "Trem2"),
  tracks = NULL,
  ONTGeneExp  = generate_plots("Trem2",annotated$targ_ont,"Gene_time","")$Trem2,
  ONTTransExp = generate_plots("Trem2",annotated$targ_ont,"Transcript_Iso_Trajectory"," ",tappassigtrans$targ_ont$TargetedOnt_Transexp)$Trem2,
  IF = plot_grid(generate_plots("Trem2",loaded$targ_ont,"IF_time_series",phenotype$targ_ont,"ONT_Targeted")[[1]]),
  pheat = draw_heatmap_gene("Trem2", class.files$targ_ont, annotated$targ_ont$Norm_transcounts)$gtable
)

# Figure 5: Bin1
Bin1_p <- list(
  dendro = dendro_plot(dirnames$targ_anno, "Bin1"),
  tracks = NULL,
  diff_tracks = NULL,
  ONTGeneExp  = generate_plots("Bin1",annotated$targ_ont,"Gene_time","")$Bin1,
  ONTTransExp = generate_plots("Bin1",annotated$targ_ont,"Transcript_Iso_Trajectory"," ",tappassigtrans$targ_ont$TargetedOnt_Transexp)$Bin1,
  IF = plot_grid(generate_plots("Bin1",loaded$targ_ont,"IF_time_series",phenotype$targ_ont,"ONT_Targeted")[[1]])
  #pheat = draw_heatmap_gene("Bin1", class.files$targ_ont, annotated$targ_ont$Norm_transcounts)$gtable
)



pdf(paste0(output_dir,"/MainFigures.pdf"), width = 20, height = 10)
plot_grid(plotlist = Gfap_p, ncol = 2, nrow = 3, labels = c("A","B","C","","D","E"))
plot_grid(plot_grid(Trem2_p$dendro,Trem2_p$tracks, ncol = 2, nrow = 1, labels = c("A","B")),
          plot_grid(Trem2_p$ONTGeneExp,Trem2_p$ONTTransExp,Trem2_p$IF,Trem2_p$pheat, ncol = 4, nrow = 1, labels = c("C","D","E","F")),
          nrow=2)
plot_grid(plotlist = Bin1_p, ncol = 3, nrow = 2, labels = c("A","B","","C","D","E"))
dev.off()