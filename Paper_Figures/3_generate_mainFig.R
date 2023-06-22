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
source(paste0(LOGEN_ROOT, "transcriptome_stats/plot_basic_stats.R"))
source(paste0(LOGEN_ROOT,"longread_QC/plot_cupcake_collapse_sensitivty.R"))
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Mouse_Isoseq/"


## ---------- Figure 1: Gfap ----------

Gfap_p <- list(
  RNAGeneExp = plot_trans_exp_individual_overtime("2973",GlobalDESeq$RresGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "RNA-Seq Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  
  IsoGeneExp = plot_trans_exp_individual_overtime("2973",GlobalDESeq$resGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "Iso-Seq Gene Expression"),
  
  tracks1 = ggTranPlots(gtf$glob_iso, class.files$glob_iso,
                        isoList = c("PB.2973.16","PB.2973.15","PB.2973.35","PB.2973.23","PB.2973.33","PB.2973.27"),
                        colours = c("#F8766D","#7CAE00","#00BFC4","gray","black",wes_palette("Zissou1")[4]),
                        lines = c("#F8766D","#7CAE00","#00BFC4","gray","black",wes_palette("Zissou1")[4])) ,
  
  # PB.2973.16 = red (#F8766D), PB.2973.15 = green (#7CAE00), PB.2973.35 = blue (#00BFC4)
  IsoTransExp  = plot_transexp_overtime("Gfap",GlobalDESeq$resTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.2973.16")) + 
    scale_colour_manual(values = c("#00BFC4","#F8766D","#7CAE00")) + theme(legend.position = "None") + 
    labs(title = "", subtitle = "Iso-Seq Transcript Expression") +
    facet_grid(cols = vars(group), labeller = labeller(group = as_labeller(c("CONTROL" = "WT", "CASE" = "TG")))),
  
  # PB.2973.16 = red
  #RNATransExp = plot_transexp_overtime("Gfap",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.2973.16")) + 
  #  theme(legend.position = "None") + 
  #  labs(title = "", subtitle = "RNA-Seq Transcript Expression"),
  
  # PB.2973.16 = red, PB.2973.23 = grey, PB.2973.33 = black, PB.2973.27 = yellow
  RNATransExpRanked  = plot_transexp_overtime("Gfap",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.2973.16")) +
    scale_colour_manual(values = c("#F8766D","gray",wes_palette("Zissou1")[4],"black")) + 
    theme(legend.position = "None") + scale_y_continuous(labels = ks) +
    labs(title = "", subtitle = "RNA-Seq Transcript Expression", y = "Normalised counts (K)") +
    facet_grid(cols = vars(group), labeller = labeller(group = as_labeller(c("CONTROL" = "WT", "CASE" = "TG"))))
)


## ---------- Figure 2: Targeted ----------

Targeted_p <- list(
  comp = whole_vs_targeted_plots(class.files$iso_match,paste0("FL.WholeIso", wholesamples), paste0("FL.TargetedIso", wholesamples), TargetGene)[[1]],
  venn = grobTree(vennONTvsIso(class.files$targ_filtered)),
  cumulative = pSensitivity(class.files$targ_all),
  num = total_num_iso(class.files$targ_filtered,""),
  splicing = plot_summarised_AS_events(Merged_gene_class_df, dirnames$targ_anno, Targeted$ref_gencode)[[1]],
  distribution = plotIFAll(Exp=Exp$targ_ont$normAll %>% select(-associated_gene),
                           classf=class.files$targ_all,
                           pheno=phenotype$targeted_rTg4510_ont,
                           majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso)) + theme(legend.position = "None")
)


## ---------- Figure 3: Trem2 ----------

Trem2ES <- ES %>% filter(associated_gene == "Trem2") %>% filter(ES %in% c("Gencode_4","Gencode_5")) 
Trem2A5A3 <- A5A3 %>% filter(associated_gene == "Trem2") %>% filter(gencode_exon %in% c("Gencode_2")) 
Trem2NE <- NovelExons %>% filter(associated_gene == "Trem2") %>% filter(novelexons %in% c("Beyond_First","Internal_NovelExon"))

Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    All = c("PB.20818.261"),
    A5A3 = setdiff(as.character(unique(Trem2A5A3$transcript_id)), as.character(unique(Trem2ES$transcript_id)))[7:12],
    ES = as.character(unique(Trem2ES$transcript_id)[1:5]),
    `Novel Exons` = as.character(unique(Trem2NE$transcript_id)[1:5]),
    DTE = c("PB.20818.54","PB.20818.55","PB.20818.62")
  )),
  Category = rep(names(Trem2Iso), lengths(Trem2Iso))
)
Trem2Iso$colour <- c(rep(NA,length(Trem2Iso$Category[Trem2Iso$Category != "DTE"])),wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")

Trem2_p <- list(
  pheat = draw_heatmap_gene(gene="Trem2",cf=class.files$targ_all,normCounts=TargetedDESeq$ontResTranAnno$wald$norm_counts,type="targeted"),
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Trem2"),
  
  # PB.20818.54 = wes_palette("Darjeeling1")[3], PB.20818.55 = blue (#00BFC4), PB.201818.62 = green (#7CAE00)
  ONTTransExp = plot_transexp_overtime("Trem2",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=3,
                                       isoSpecific=c("PB.20818.54"),setorder=c("CONTROL","CASE")) + 
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")) +
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks)+
    theme(legend.position = "None"),
  ONTGeneExp = plot_trans_exp_individual_overtime("20818",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("Trem2",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL,isoSpecific=c("PB.20818.54")),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(Trem2Iso$Isoform)),
                       selfDf = Trem2Iso, gene = "Trem2")
  
)
Trem2_p$IF[[1]] <- Trem2_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
Trem2_p$IF[[2]] <- Trem2_p$IF[[2]] + labs(title = "")


## ---------- Figure 4: Bin1 ---------

Bin1ES <- ES %>% filter(associated_gene == "Bin1") 

Bin1Iso <- data.frame(
  Isoform = unlist(Bin1Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Bin1" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    ES = as.character(unique(Bin1ES$transcript_id)[6:10]),
    IR = c("PB.22007.45261","PB.22007.46397","PB.22007.45271","PB.22007.924","PB.22007.803"),
    AP = c(as.character(Bin1IRFirstExon$V4))[46:56],
    DTE = c("PB.22007.101","PB.22007.224","PB.22007.99")
  )),
  Category = rep(names(Bin1Iso), lengths(Bin1Iso))
)
Bin1Iso$colour <- c(rep(NA,length(Bin1Iso$Category[Bin1Iso$Category != "DTE"])),wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")

Bin1_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Bin1"),
  ES = plot_ES_Tgene(dirnames$targ_anno,"Bin1",class.files$targ_filtered)[[1]],
  IR = plot_IR_Tgene(dirnames$targ_anno,"Bin1",class.files$targ_filtered)[[2]] + theme(legend.position = c(0.25,0.9)),
  ONTGeneExp  = plot_trans_exp_individual_overtime("22007",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  ONTTransExp = plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=2,isoSpecific=c("PB.22007.224"),
                                       setorder=c("CONTROL","CASE")) + 
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")) +
    theme(legend.position = "None"),
  IF = plotIF("Bin1",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",majorIso=NULL,isoSpecific=c("PB.22007.224","PB.22007.99","PB.22007.176")),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(Bin1Iso$Isoform)),
                       selfDf = Bin1Iso, gene = "Bin1")
)
Bin1_p$IF[[1]] <- Bin1_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
Bin1_p$IF[[2]] <- Bin1_p$IF[[2]] + labs(title = "")


## ---------- Figure 5: Clu ---------

CluES <- ES %>% filter(associated_gene == "Clu") 

CluIso <- data.frame(
  Isoform = unlist(CluIso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Clu" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    ES = c("PB.14646.1019","PB.14646.4411","PB.14646.1094","PB.14646.4585","PB.14646.11450"),
    AP = c("PB.14646.134",as.character(CluAF$V4[1:10])),
    IR = c("PB.14646.4741","PB.14646.4624","PB.14646.1317","PB.14646.2515"),
    DTE = c("PB.14646.139","PB.14646.39341")
  )),
  Category = rep(names(CluIso), lengths(CluIso))
)
CluIso$colour <- c(rep(NA,length(CluIso$Category[CluIso$Category != "DTE"])),wes_palette("Darjeeling1")[3],"#00BFC4")


Clu_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Clu"),
  ES = plot_ES_Tgene(dirnames$targ_anno,"Clu",class.files$targ_filtered)[[1]],
  IR = plot_IR_Tgene(dirnames$targ_anno,"Clu",class.files$targ_filtered)[[2]],
  ONTTransExp = plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=1,isoSpecific=c("PB.14646.39341"),setorder=c("CONTROL","CASE"))+
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks)  +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4")) +
    theme(legend.position = "None"),
  ONTGeneExp = plot_trans_exp_individual_overtime("14646",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("Clu",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=NULL,majorIso=NULL,isoSpecific=c("PB.14646.139","PB.14646.60837","PB.14646.39341")),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(CluIso$Isoform)),
                       selfDf = CluIso, gene = "Clu")
)
Clu_p$IF[[1]] <- Clu_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
Clu_p$IF[[2]] <- Clu_p$IF[[2]] + labs(title = "")


## ---------- Figure 6: App ---------

AppES <- ES %>% filter(associated_gene == "App") %>% filter(ES %in% paste0("Gencode_",c(7,8,14,15,16)))

AppIso <- data.frame(
  Isoform = unlist(AppIso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "App" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    ES = c("PB.19309.10098","PB.19309.10232","PB.19309.12318","PB.19309.28614","	PB.19309.11446","PB.19309.12434","PB.19309.31219","PB.19309.31372"),
    DTE = c("PB.19309.7368","PB.19309.7374","PB.19309.7389")
  )),
  Category = rep(names(AppIso), lengths(AppIso))
)
AppIso$colour <- c(rep(NA,length(AppIso$Category[AppIso$Category != "DTE"])),wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")


App_p <- list(
  ES = plot_ES_Tgene(dirnames$targ_anno,"App",class.files$targ_filtered)[[1]],
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "App"),
  ONTTransExp = plot_transexp_overtime("App",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=3,setorder=c("CONTROL","CASE")) +
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")) +
    theme(legend.position = "None"),
  ONTGeneExp = plot_trans_exp_individual_overtime("19309",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("App",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(AppIso$Isoform)),
                       selfDf = AppIso, gene="App")
)
App_p$IF[[1]] <- App_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
App_p$IF[[2]] <- App_p$IF[[2]] + labs(title = "")


## ---------- Figure 7: Apoe ---------

ApoeIso <- data.frame(
  Isoform = unlist(ApoeIso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Apoe" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    A5A3 = as.character(unique(ApoeA5A3[ApoeA5A3$gencode_exon == "Gencode_8","transcript_id"])[3:25]),
    ES = setdiff(as.character(unique(ES[ES$associated_gene == "Apoe" & ES$ES %in% c(paste0("Gencode_",c(5,6))),"transcript_id"])),
                              as.character(unique(ApoeA5A3[ApoeA5A3$gencode_exon == "Gencode_8","transcript_id"])[1:30]))[5:15],
    DTE = c("PB.40586.871","PB.40586.872","PB.40586.875")
  )),
  Category = rep(names(ApoeIso), lengths(ApoeIso))
)
ApoeIso$colour <- c(rep(NA,length(ApoeIso$Category[ApoeIso$Category != "DTE"])),wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")


Apoe_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Apoe"),
  A5A3 = plot_A5A3_Tgene(dirnames$targ_anno,"Apoe") + theme(legend.position = c(0.3,0.8)),
  ONTTransExp = plot_transexp_overtime("Apoe",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=3,setorder=c("CONTROL","CASE")) +
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")) +
    theme(legend.position = "None"),
  ONTGeneExp = plot_trans_exp_individual_overtime("40586",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalised counts (K)") + scale_y_continuous(labels = ks),
  IF = plotIF("Apoe",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(ApoeIso$Isoform)),
                       selfDf = ApoeIso, gene="Apoe")
)
Apoe_p$IF[[1]] <- Apoe_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
Apoe_p$IF[[2]] <- Apoe_p$IF[[2]] + labs(title = "")

## ---------- Output Pdf -----------------

pdf(paste0(output_dir,"/MainFigures1.pdf"), width = 10, height = 15)
plot_grid(
plot_grid(NULL,labels=c("A")),
plot_grid(plotlist = Gfap_p[1:2], labels = c("B","C")),
plot_grid(Gfap_p$tracks1,labels=c("D")),
plot_grid(Gfap_p$RNATransExpRanked,Gfap_p$IsoTransExp, labels = c("E","F")), ncol = 1
)
dev.off()

pdf(paste0(output_dir,"/MainFigures2.pdf"), width = 15, height = 15)
plot_grid(plot_grid(Targeted_p$comp,Targeted_p$cumulative,Targeted_p$venn,nrow=1,labels=c("A","B","C"),rel_widths = c(0.4,0.3,0.3)),
          plot_grid(Targeted_p$distribution,Targeted_p$n,nrow=1,labels=c("D","E"),rel_widths = c(0.6,0.4)),
          plot_grid(Targeted_p$splicing,nrow=1,labels=c("F")),nrow=3, rel_heights = c(0.3,0.4,0.3))
dev.off()

pdf(paste0(output_dir,"/MainFigures3.pdf"), width = 18, height = 12)
plot_grid(plot_grid(
            plot_grid(Trem2_p$dendro,Trem2_p$pheat$gtable,nrow=1, rel_widths = c(0.5,0.5), labels = c("A","B")),
            plot_grid(Trem2_p$ONTTransExp,Trem2_p$ONTGeneExp,nrow=1, rel_widths = c(0.6,0.4), labels = c("D","E")),
            plot_grid(Trem2_p$IF[[1]],Trem2_p$IF[[2]],nrow=1, rel_widths = c(0.6,0.4),labels = c("F","G")), ncol = 1),
          Trem2_p$tracks,
          ncol = 2, rel_widths = c(0.5,0.5), labels = c("","C")
)
dev.off()

pdf(paste0(output_dir,"/MainFigures4.pdf"), width = 18, height = 14)
plot_grid(plot_grid(
  plot_grid(Bin1_p$dendro,plot_grid(Bin1_p$ES,Bin1_p$IR,ncol=1,labels = c("B","C")),nrow = 1,labels = c("A",""), rel_widths = c(0.4,0.6)),
  plot_grid(Bin1_p$ONTTransExp, Bin1_p$ONTGeneExp,rel_widths = c(0.6,0.4),labels = c("E","F")),
  plot_grid(Bin1_p$IF[[1]], Bin1_p$IF[[2]],rel_widths = c(0.6,0.4),labels = c("G","H")),
  ncol = 1, rel_heights = c(0.45,0.275,0.275)), Bin1_p$tracks,
  ncol = 2, rel_widths = c(0.5,0.5), labels = c("","D")
)
dev.off()

pdf(paste0(output_dir,"/MainFigures5.pdf"), width = 18, height = 12)
plot_grid(plot_grid(
    plot_grid(Clu_p$dendro,Clu_p$ES, nrow = 1, labels = c("A","B")),
    plot_grid(Clu_p$ONTTransExp,Clu_p$ONTGeneExp, nrow = 1, rel_widths = c(0.6,0.4), labels = c("D","E")),
    plot_grid(Clu_p$IF[[1]],Clu_p$IF[[2]], nrow = 1, rel_widths = c(0.6,0.4), labels = c("F","G")), ncol = 1), Clu_p$tracks, 
  ncol = 2, rel_widths = c(0.5,0.5), labels = c("","C"))
dev.off()

pdf(paste0(output_dir,"/MainFigures6.pdf"), width = 18, height = 12)
plot_grid(plot_grid(
    plot_grid(App_p$dendro,App_p$ES, nrow = 1, labels = c("A","B")),
    plot_grid(App_p$ONTTransExp,App_p$ONTGeneExp, nrow = 1, rel_widths = c(0.6,0.4), labels = c("D","E")),
    plot_grid(App_p$IF[[1]],App_p$IF[[2]], nrow = 1, rel_widths = c(0.6,0.4), labels = c("F","G")),ncol = 1),
  App_p$tracks, rel_widths = c(0.5,0.4), ncol = 2, labels = c("","C"))
dev.off()

pdf(paste0(output_dir,"/MainFigures7.pdf"), width = 18, height = 12)
plot_grid(plot_grid(
  plot_grid(Apoe_p$dendro, Apoe_p$A5A3,nrow=1,labels=c("A","B")),
  plot_grid(Apoe_p$ONTTransExp,Apoe_p$ONTGeneExp,rel_widths = c(0.6,0.4),labels=c("D","E")),
  plot_grid(Apoe_p$IF[[1]],Apoe_p$IF[[2]],rel_widths = c(0.6,0.4),labels=c("F","G")),
  ncol=1),
  Apoe_p$tracks, rel_widths = c(0.5,0.4), ncol = 2, labels = c("","C")
)
dev.off()

