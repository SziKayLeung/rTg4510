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

LOGEN = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(SC_ROOT, "rTg4510_config.R"))
output_dir = "/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/01_figures_tables/Mouse_Isoseq/"


## ---------- Figure 2: Gfap ----------

Gfap_p <- list(
  RNAGeneExp = plot_trans_exp_individual_overtime("2973",GlobalDESeq$RresGeneAnno$lrt$norm_counts,type="gene") + 
    theme(legend.position = c(0.15,0.9)) +
    theme(legend.title=element_blank()) + 
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG")) +
    labs(title = "", subtitle = "RNA-Seq Gene Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks),
  
  IsoGeneExp = plot_trans_exp_individual_overtime("2973",GlobalDESeq$resGeneAnno$lrt$norm_counts,type="gene") + 
    theme(legend.position = c(0.15,0.9)) +
    theme(legend.title=element_blank()) + 
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG")) +
    labs(title = "", subtitle = "Iso-Seq Gene Expression"),
  
  tracks1 = ggTranPlots(inputgtf=gtf$glob_iso, classfiles=class.files$glob_iso,
                        isoList = c("Gfap-201","PB.2973.16","PB.2973.15","PB.2973.35","PB.2973.23","PB.2973.33","PB.2973.27"),
                        colours = c("black", wes_palette("Darjeeling2")[2],"#F8766D","#7CAE00", wes_palette("Zissou1")[4], wes_palette("GrandBudapest2")[2], "#00BFC4"),
                        lines = c("black", wes_palette("Darjeeling2")[2],"#F8766D","#7CAE00", wes_palette("Zissou1")[4], wes_palette("GrandBudapest2")[2], "#00BFC4"),
                        gene = "Gfap")
  ,
  
  # PB.2973.16 = red (#F8766D), PB.2973.15 = wes_palette("Darjeeling2")[2], PB.2973.35 = blue (#00BFC4)
  IsoTransExp  = plot_transexp_overtime("Gfap",GlobalDESeq$resTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.2973.16"), setorder = c("CONTROL","CASE")) + 
    scale_colour_manual(values = c("#00BFC4","#F8766D",wes_palette("Darjeeling2")[2])) + theme(legend.title=element_blank()) + 
    theme(legend.position = c(0.2,0.8)) +
    labs(title = "", subtitle = "Iso-Seq Transcript Expression") +
    facet_grid(cols = vars(group)),
  
  # PB.2973.16 = red
  RNATransExp = plot_transexp_overtime("Gfap",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.2973.16"), setorder = c("CONTROL","CASE")) + 
    theme(legend.position = "None") + 
    labs(title = "", subtitle = "RNA-Seq Transcript Expression"),
  
  # PB.2973.16 = red, PB.2973.23 = green (#7CAE00), PB.2973.33 = purple, PB.2973.27 = yellow
  RNATransExpRanked  = plot_transexp_overtime("Gfap",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.2973.16"), setorder = c("CONTROL","CASE")) +
    scale_colour_manual(values = c("#F8766D","#7CAE00",wes_palette("Zissou1")[4],wes_palette("GrandBudapest2")[2])) + 
    theme(legend.position = c(0.2,0.75)) +
    theme(legend.title=element_blank()) + scale_y_continuous(labels = ks) +
    labs(title = "", subtitle = "RNA-Seq Transcript Expression", y = "Normalized counts (K)") +
    facet_grid(cols = vars(group)),
  
  Sorted = plot_boxplot_SCN(Exp$targ_sorted_all,"PB.39126.482") + theme(legend.title=element_blank()) + labs(title = "", subtitle = "Gfap-201 (LR.Gfap.16) Transcript Expression") + theme(legend.position = c(0.15,0.8))
)
Gfap_p$tracks1 <- Gfap_p$tracks1 + theme(legend.position = "None", 
                       axis.line.x = element_line(colour = "grey80"),
                       panel.background = element_rect(fill = "white", colour = "grey50"),
                       panel.border = element_rect(fill = NA, color = "grey50", linetype = "dotted"),
                       axis.text.y= element_text(size=12),
                       strip.text.y = element_text(size = 12, color = "black"),
                       strip.background = element_rect(fill = "white", colour = "grey50"))

## ---------- Figure 3: Targeted ----------


ApoeIso <- data.frame(
  Isoform = unlist(ApoeIso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Apoe" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    #Reference = c("ENSMUST00000167646.8","ENSMUST00000003066.15"),
    A5A3 = as.character(unique(ApoeA5A3[ApoeA5A3$gencode_exon == "Gencode_8","transcript_id"])[3:23]))),
  Category = rep(names(ApoeIso), lengths(ApoeIso))
)
ApoeIso$colour <- c(rep(NA,length(ApoeIso$Category[ApoeIso$Category != "DTE"])))
Apoe_p = ggTranPlots(inputgtf=gtf$targ_merged, classfiles=class.files$targ_filtered,
                   isoList = c(as.character(ApoeIso$Isoform)),
                   selfDf = ApoeIso, gene="Apoe", inputPfam = Targeted$pfam, 
                   rnaseqDir=dirnames$rna_aligned, rnaseqTransID="ENSMUST00000174064.8")

Targeted_p <- list(
  comp = whole_vs_targeted_plots(class.files$iso_match,paste0("FL.WholeIso", wholesamples), paste0("FL.TargetedIso", wholesamples), TargetGene)[[1]],
  venn = grobTree(vennONTvsIso(class.files$targ_all)),
  cumulative = pSensitivity(class.files$targ_all),
  num = total_num_iso(class.files$targ_filtered,"") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "Gene"),
  splicing = plot_summarised_AS_events(Merged_gene_class_df, dirnames$targ_anno, Targeted$ref_gencode)[[1]],
  distribution = plotIFAll(Exp=Exp$targ_ont$normAll %>% select(-associated_gene),
                           classf=class.files$targ_all,
                           pheno=phenotype$targeted_rTg4510_ont,
                           majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso))[1],
  apoe = plot_grid(Apoe_p[[1]] + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
                             plot.margin = unit(c(0.5, 0, 0, 0), "cm")),
                   Apoe_p[[2]] + theme(plot.margin = unit(c(0.1, 0, 0, 0.22), "cm")),ncol=1, rel_heights = c(0.85,0.15))
)


## ---------- Figure 4: Differential expression ----------

OntTargetedDiff <- list()
for(num in 1:length(TargetedDESeq$ontResTranAnno$wald$anno_res$isoform)){
  i = TargetedDESeq$ontResTranAnno$wald$anno_res$isoform[num]
  g = class.files$targ_all[class.files$targ_all$isoform == i,"associated_gene"][[1]]
  OntTargetedDiff[[num]] = plot_transexp_overtime(g,TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",
                                                  rank=0,isoSpecific=i,setorder=c("CONTROL","CASE"),
                                                  classfiles=class.files$targ_filtered) + labs(x = "", y = "")
}
OntTargetedDiff[[13]] <- plot_transexp_overtime("Mapt",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",
                                                rank=0,isoSpecific="PB.8675.37810",setorder=c("CONTROL","CASE"),
                                                classfiles=class.files$targ_filtered) + labs(x = "", y = "")
OntTargetedDiff[[14]] <- plot_transexp_overtime("Mapt",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",
                                                rank=0,isoSpecific="PB.8675.41059",setorder=c("CONTROL","CASE"),
                                                classfiles=class.files$targ_filtered) + labs(x = "", y = "")


RefIsoforms <- lapply(c("Clu","Fyn","Apoe","Trem2","Cd33","App","Mapt"), function(x) unique(gtf$ref_target[gtf$ref_target$gene_name == x & !is.na(gtf$ref_target$transcript_id), "transcript_id"]))
names(RefIsoforms ) <- c("Clu","Fyn","Apoe","Trem2","Cd33","App","Mapt")

OntTargetedTracks <- list(
  Trem2 = ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
                      isoList = c("PB.20818.54","PB.20818.62","PB.20818.55",RefIsoforms$Trem2[1:2]),
                      colours = c(wes_palette("Cavalcanti1")[5],wes_palette("Royal1")[2],alpha(wes_palette("Royal1")[2],0.3),rep("#0C0C78",2)), 
                      lines = c(wes_palette("Cavalcanti1")[5],wes_palette("Royal1")[2],alpha(wes_palette("Royal1")[2],0.3),rep("#0C0C78",2)),
                      #lines = c(wes_palette("Rushmore1")[3],wes_palette("Royal2")[5],alpha(wes_palette("Rushmore1")[3],0.3),rep("#0C0C78",2)), 
                      gene = "Trem2",simple=TRUE),
  Clu = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
                    isoList = c("PB.14646.39341","PB.14646.139","PB.14646.39352","PB.14646.35283",RefIsoforms$Clu[1]),
                    colours = c(wes_palette("Darjeeling2")[2],wes_palette("Zissou1")[1],
                                wes_palette("Darjeeling1")[5],wes_palette("Darjeeling2")[4],rep("#0C0C78",2)), 
                    lines = c(wes_palette("Rushmore1")[3],wes_palette("Zissou1")[1],
                              wes_palette("Darjeeling1")[5],wes_palette("Darjeeling2")[4],rep("#0C0C78",2)), 
                    gene = "Clu",simple=TRUE),
  Fyn = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
                    isoList = c("PB.3948.4511",RefIsoforms$Fyn[2]),
                    colours = c(wes_palette("Royal2")[5],rep("#0C0C78",2)), lines = c(wes_palette("Royal2")[5],rep("#0C0C78",2)), gene = "Fyn", simple=TRUE), 
  Apoe = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
                     isoList = c("PB.40586.1023","PB.40586.875",RefIsoforms$Apoe[1]),
                     colours = c(wes_palette("Royal1")[4],alpha(wes_palette("Royal1")[4],0.4),rep("#0C0C78",2)), 
                     lines = c(wes_palette("Royal1")[4],alpha(wes_palette("Royal1")[4],0.4),rep("#0C0C78",2)), 
                     gene = "Apoe", simple=TRUE),
  Cd33 = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
                     isoList = c("PB.41115.1365",RefIsoforms$Cd33[1:2]),
                     colours = c(wes_palette("GrandBudapest2")[1],rep("#0C0C78",2)), lines = c(wes_palette("GrandBudapest2")[1],rep("#0C0C78",2)), gene = "Cd33", simple=TRUE),
  App = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
                    isoList = c("PB.19309.7564",RefIsoforms$App[1]),
                    colours = c(wes_palette("GrandBudapest2")[2],rep("#0C0C78",2)), lines = c(wes_palette("GrandBudapest2")[2],rep("#0C0C78",2)), gene = "App", simple=TRUE),
  Mapt = ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
                     isoList = c("PB.8675.37810","PB.8675.41059","ENSMUST00000106992.9","ENSMUST00000106993.9" ),
                     colours = c(wes_palette("Darjeeling2")[5],alpha(wes_palette("Darjeeling2")[5],0.3),rep("#0C0C78",2)), 
                     lines = c(wes_palette("Darjeeling2")[5],alpha(wes_palette("Darjeeling2")[5],0.3),rep("#0C0C78",2)),
                     #lines = c(wes_palette("Rushmore1")[3],wes_palette("Royal2")[5],alpha(wes_palette("Rushmore1")[3],0.3),rep("#0C0C78",2)), 
                     gene = "Mapt",simple=TRUE)
)

OntTargetedDiffColours <- data.frame(
  iso = c(
    TargetedDESeq$ontResTranAnno$wald$anno_res$isoform, 
    "PB.8675.37810",
    "PB.8675.41059"
  ),
  col = c(
    wes_palette("Darjeeling2")[2],
    wes_palette("Zissou1")[1],
    wes_palette("Cavalcanti1")[5],
    wes_palette("Royal1")[2],
    wes_palette("Darjeeling1")[5],
    wes_palette("Royal2")[5],
    wes_palette("Darjeeling2")[4],
    wes_palette("Royal1")[4],
    wes_palette("GrandBudapest2")[1],
    wes_palette("GrandBudapest2")[2],
    alpha(wes_palette("Royal1")[4], 0.4),
    alpha(wes_palette("Royal1")[2], 0.3),
    wes_palette("Darjeeling2")[5],
    alpha(wes_palette("Darjeeling2")[5], 0.3)
  )
)

for(i in 1:nrow(OntTargetedDiffColours)){
  #OntTargetedDiff[[i]] = OntTargetedDiff[[i]] + theme(plot.title = element_text(colour = as.character(OntTargetedDiffColours$col[[i]])))
  OntTargetedDiff[[i]] = OntTargetedDiff[[i]] + facet_grid(. ~  LRID_struc) + theme(strip.background = element_rect(fill=as.character(OntTargetedDiffColours$col[[i]]))) + 
    labs(title = NULL, x=NULL, y=NULL,subtitle="")
  if(i %in% c(1,2,3,4,13)){
    OntTargetedDiff[[i]] = OntTargetedDiff[[i]] + theme(strip.text = element_text(size=15, colour="white"))
  }else{
    OntTargetedDiff[[i]] = OntTargetedDiff[[i]] + theme(strip.text = element_text(size=15, colour="black"))
  }
}

tag_ggplot_seq <- function(p, id){
  p1 <- p + labs(tag = as.character(id)) + theme(text = element_text(size = 12))
  return(p1)
}
OntTargetedDiff <- lapply(seq_len(14), function(x) tag_ggplot_seq(OntTargetedDiff[[x]],x))
OntTargetedDiff[[13]] <- OntTargetedDiff[[13]] + labs(tag = as.character(1)) + theme(text = element_text(size = 12))
OntTargetedDiff[[14]] <- OntTargetedDiff[[14]] + labs(tag = as.character(2)) + theme(text = element_text(size = 12))


Diffa <- arrangeGrob(grobs=lapply(OntTargetedDiff, function(p) p + guides(colour=FALSE)), ncol=2, 
              bottom=textGrob("Age (months)", gp=gpar(fontsize=15)), 
              left=textGrob("Normalized counts", gp=gpar(fontsize=15), rot=90),
              labels = c("a"))

Diffb <- plot_grid(plotlist = OntTargetedTracks,ncol=1,rel_heights = c(0.3,0.3,0.1,0.15,0.15,0.1,0.25))

Diff <- plot_grid(Diffa,Diffb)
plot_grid(Diff)

## ---------- Figure 4: Trem2 ----------

Trem2ES <- ES %>% filter(associated_gene == "Trem2") %>% filter(ES %in% c("Gencode_4","Gencode_5")) 

Trem2ES5 <- ES %>% filter(associated_gene == "Trem2") %>% filter(ES %in% c("Gencode_5")) 
Trem2A5A3 <- A5A3 %>% filter(associated_gene == "Trem2") %>% filter(gencode_exon %in% c("Gencode_2")) 
Trem2NE <- NovelExons %>% filter(associated_gene == "Trem2") %>% filter(novelexons %in% c("Beyond_First","Internal_NovelExon"))

Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    All = c("PB.20818.261"),
    A5A3 = setdiff(as.character(unique(Trem2A5A3$transcript_id)), as.character(unique(Trem2ES$transcript_id)))[7:12],
    ES = as.character(unique(Trem2ES$transcript_id)[1:5]),
    NE = c("PB.20818.16","PB.20818.32"),
    CE = paste0("PB.20818.", c("80","192","573","1074","493","362","1096")),
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
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks)+
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  ONTGeneExp = plot_trans_exp_individual_overtime("20818",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG"), name = "") +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  IF = plotIF("Trem2",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL,isoSpecific=c("PB.20818.54","PB.20818.55","PB.20818.62")),
  tracks = ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
                       isoList = c(as.character(Trem2Iso$Isoform)),
                       selfDf = Trem2Iso, gene = "Trem2", inputPfam=Targeted$pfam), 
  
  sorted1 = plot_boxplot_SCN(Exp$targ_sorted_all,iso=c("PB.95419.5","PB.95419.7"), ageDiv = FALSE) + 
    labs(title = "", subtitle = "ONT Transcript Expression") + theme(legend.position = c(0.9,0.8)) +  
    facet_grid(~isoform, labeller = labeller(isoform = as_labeller(c("PB.95419.5" = "Trem2-201", "PB.95419.7" = "Trem2-202"))))
  
)
Trem2_p$IF[[1]] <- Trem2_p$IF[[1]] + labs(title = "", y = "IF (%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) + guides(fill = FALSE) 
Trem2_p$IF[[2]] <- Trem2_p$IF[[2]] + labs(title = "", y = "IF (%)") + guides(colour = FALSE) 

# stats
plot_boxplot_SCN(Exp$targ_sorted_all,iso=c("PB.95419.5"), ageDiv = FALSE)
plot_boxplot_SCN(Exp$targ_sorted_all,iso=c("PB.95419.7"), ageDiv = FALSE)

## ---------- Figure 5: Bin1 ---------

Bin1ES <- ES %>% filter(associated_gene == "Bin1") %>% filter(ES %in% c("Gencode_13","Gencode_15","Gencode_14","Gencode_16")) 

Bin1Iso <- data.frame(
  Isoform = unlist(Bin1Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Bin1" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    ES = c("PB.22007.100","PB.22007.1005","PB.22007.423","PB.22007.10118","PB.22007.10147","PB.22007.42921"),
    AP = paste0("PB.22007.", c(50820, 50175, 50143, 50079, 49307, 48852, 48704, 48337, 47598)),
    CE = paste0("PB.22007.",c("925","1554","1033","4967","1470","14222","1014")),
    DTE = c("PB.22007.224","PB.22007.99")
  )),
  Category = rep(names(Bin1Iso), lengths(Bin1Iso))
)
Bin1Iso$colour <- c(rep(NA,length(Bin1Iso$Category[Bin1Iso$Category != "DTE"])),wes_palette("Darjeeling1")[3],"#00BFC4")

Bin1_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Bin1"),
  ES = plot_ES_Tgene(dirnames$targ_anno,"Bin1",class.files$targ_filtered)[[1]],
  ONTGeneExp  = plot_trans_exp_individual_overtime("22007",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG"), name = "") +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  ONTTransExp = plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=NULL,
                                       isoSpecific=c("PB.22007.224","PB.22007.99"),
                                       setorder=c("CONTROL","CASE")) + 
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4")) +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  IF = plotIF("Bin1",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",majorIso=NULL,isoSpecific=c("PB.22007.224"),rank=0),
  tracks = ggTranPlots(inputgtf=gtf$targ_merged, classfiles=class.files$targ_filtered,
                       isoList = c(as.character(Bin1Iso$Isoform)),
                       selfDf = Bin1Iso, gene = "Bin1", inputPfam=Targeted$pfam)
)


## ---------- Figure 6: Clu ---------

CluES <- ES %>% filter(associated_gene == "Clu") 

CluIso <- data.frame(
  Isoform = unlist(CluIso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Clu" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    ES = c("PB.14646.1019","PB.14646.4411","PB.14646.1094","PB.14646.4585","PB.14646.11450"),
    AP = paste0("PB.14646.",c("39554","39402","39345","38679","34808","38406","38351","34798","329","1380","134")),
    IR = c("PB.14646.4741","PB.14646.4624","PB.14646.1317","PB.14646.2515"),
    CE = paste0("PB.14646.",c("6992")),
    DTE = c("PB.14646.139","PB.14646.39341","PB.14646.39352","PB.14646.35283")
  )),
  Category = rep(names(CluIso), lengths(CluIso))
)
CluIso$colour <- c(rep(NA,length(CluIso$Category[CluIso$Category != "DTE"])),wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00",wes_palette("IsleofDogs1")[1])


Clu_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "Clu"),
  ES = plot_ES_Tgene(dirnames$targ_anno,"Clu",class.files$targ_filtered)[[1]],
  IR = plot_IR_Tgene(dirnames$targ_anno,"Clu",class.files$targ_filtered)[[2]],
  ONTTransExp = plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="toprank",rank=1,isoSpecific=c("PB.14646.139","PB.14646.39341","PB.14646.39352","PB.14646.35283"),setorder=c("CONTROL","CASE"))+
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks)  +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00",wes_palette("IsleofDogs1")[1])) +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  ONTGeneExp = plot_trans_exp_individual_overtime("14646",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG"), name = "") +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  IF = plotIF("Clu",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=0,majorIso=NULL,isoSpecific=c("PB.14646.139")),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(CluIso$Isoform)),
                       selfDf = CluIso, gene = "Clu", inputPfam=Targeted$pfam)
)
Clu_p$IF[[1]] <- Clu_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), legend.position = "None")   
Clu_p$IF[[2]] <- Clu_p$IF[[2]] + labs(title = "")


## ---------- Figure 7: App ---------

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
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")) +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  ONTGeneExp = plot_trans_exp_individual_overtime("19309",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG"), name = "") +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  IF = plotIF("App",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(AppIso$Isoform)),
                       selfDf = AppIso, gene="App")
)
App_p$IF[[1]] <- App_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) + guides(fill = FALSE) 
App_p$IF[[2]] <- App_p$IF[[2]] + labs(title = "")


## ---------- Figure 8: Apoe ---------

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
    labs(title = "", subtitle = "ONT Transcript Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00")) +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  ONTGeneExp = plot_trans_exp_individual_overtime("40586",TargetedDESeq$ontResGeneAnno$wald$norm_counts,type="gene") + 
    labs(title = "", subtitle = "ONT Gene Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG"), name = "") +
    theme(legend.title=element_blank(), legend.justification = c(0, 1), legend.position = "top"),
  IF = plotIF("Apoe",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
              cfiles=class.files$targ_all,design="time_series",rank=3,majorIso=NULL),
  tracks = ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                       isoList = c(as.character(ApoeIso$Isoform)),
                       selfDf = ApoeIso, gene="Apoe")
)
Apoe_p$IF[[1]] <- Apoe_p$IF[[1]] + labs(title = "") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) + guides(fill = FALSE) 
Apoe_p$IF[[2]] <- Apoe_p$IF[[2]] + labs(title = "")

## ---------- Output Pdf -----------------

pdf(paste0(output_dir,"/MainFigures2.pdf"), width = 10, height = 12)
plot_grid(
plot_grid(plotlist = Gfap_p[2:1], labels = c("A","B")),
plot_grid(Gfap_p$tracks1,Gfap_p$IsoTransExp, labels=c("C","D")),
plot_grid(Gfap_p$RNATransExpRanked, Gfap_p$Sorted, labels = c("E","F")), 
ncol = 1
)
dev.off()

pdf(paste0(output_dir,"/MainFigures3.pdf"), width = 15, height = 20)
plot_grid(plot_grid(Targeted_p$comp,Targeted_p$cumulative,Targeted_p$venn,nrow=1,labels=c("A","B","C"),rel_widths = c(0.4,0.3,0.3)),
          plot_grid(Targeted_p$distribution[[1]],Targeted_p$n,nrow=1,labels=c("D","E"),rel_widths = c(0.55,0.45)),
          plot_grid(Targeted_p$apoe, labels = c("F")),
          plot_grid(Targeted_p$splicing,nrow=1,labels=c("G")),
          nrow=4, rel_heights = c(0.225,0.225,0.4,0.15))
dev.off()

pdf(paste0(output_dir,"/MainFigures3b.pdf"), width = 15, height = 20)
plot_grid(Diffa,Diffb, labels = c("A","B"),label_size = 20)
dev.off()

pdf(paste0(output_dir,"/MainFigures4.pdf"), width = 18, height = 16)
plot_grid(plot_grid(
            plot_grid(Trem2_p$dendro,Trem2_p$pheat$gtable,nrow=1, rel_widths = c(0.5,0.5), labels = c("A","B"),  scale = 0.95),
            plot_grid(Trem2_p$sorted, NULL,nrow=1, rel_widths = c(0.5,0.5), labels = c("D","E")),
            plot_grid(Trem2_p$ONTTransExp,Trem2_p$ONTGeneExp,nrow=1, rel_widths = c(0.55,0.45), labels = c("F","G")),
            plot_grid(Trem2_p$IF[[1]],Trem2_p$IF[[2]],nrow=1, rel_widths = c(0.55,0.45),labels = c("H","I")), ncol = 1, 
            rel_heights = c(0.25,0.3,0.25,0.25)),
          Trem2_p$tracks,
          ncol = 2, rel_widths = c(0.55,0.45), labels = c("","C")
)
dev.off()

pdf(paste0(output_dir,"/MainFigures5.pdf"), width = 18, height = 12)
plot_grid(plot_grid(
  plot_grid(Clu_p$dendro, nrow = 1, labels = c("A")),
  plot_grid(Clu_p$ONTTransExp,Clu_p$ONTGeneExp, nrow = 1, rel_widths = c(0.6,0.4), labels = c("C","D")),
  plot_grid(Clu_p$IF[[1]],Clu_p$IF[[2]], nrow = 1, rel_widths = c(0.6,0.4), labels = c("E","F")), ncol = 1), Clu_p$tracks, 
  ncol = 2, rel_widths = c(0.6,0.4), labels = c("","B"))
dev.off()

pdf(paste0(output_dir,"/MainFigures6.pdf"), width = 18, height = 14)
plot_grid(plot_grid(
  plot_grid(Bin1_p$dendro,Bin1_p$ES,ncol=1,labels = c("A","B"),rel_heights = c(0.5,0.5)),
  plot_grid(Bin1_p$ONTTransExp, Bin1_p$ONTGeneExp,rel_widths = c(0.55,0.45),labels = c("D","E")),
  plot_grid(Bin1_p$IF[[1]] + theme(legend.position = "None"), Bin1_p$IF[[2]],rel_widths = c(0.55,0.45),labels = c("F","G")),
  ncol = 1, rel_heights = c(0.45,0.275,0.275)), Bin1_p$tracks,
  ncol = 2, rel_widths = c(0.5,0.5), labels = c("","C")
)
dev.off()

pdf(paste0(output_dir,"/MainFigures7.pdf"), width = 20, height = 12)
plot_grid(plot_grid(
    plot_grid(App_p$dendro,App_p$ES, nrow = 1, labels = c("A","B")),
    plot_grid(App_p$ONTTransExp,App_p$ONTGeneExp, nrow = 1, rel_widths = c(0.6,0.4), labels = c("D","E")),
    plot_grid(App_p$IF[[1]],App_p$IF[[2]], nrow = 1, rel_widths = c(0.6,0.4), labels = c("F","G")),ncol = 1),
  App_p$tracks, rel_widths = c(0.5,0.4), ncol = 2, labels = c("","C"))
dev.off()

pdf(paste0(output_dir,"/MainFigures8.pdf"), width = 18, height = 12)
plot_grid(plot_grid(
  plot_grid(Apoe_p$dendro, Apoe_p$A5A3,nrow=1,labels=c("A","B")),
  plot_grid(Apoe_p$ONTTransExp,Apoe_p$ONTGeneExp,rel_widths = c(0.6,0.4),labels=c("D","E")),
  plot_grid(Apoe_p$IF[[1]],Apoe_p$IF[[2]],rel_widths = c(0.6,0.4),labels=c("F","G")),
  ncol=1),
  Apoe_p$tracks, rel_widths = c(0.5,0.4), ncol = 2, labels = c("","C")
)
dev.off()

