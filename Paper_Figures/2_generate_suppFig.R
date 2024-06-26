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
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Mouse_Isoseq/"

# heatmap 
#draw_heatmap_gene_level(loaded$glob_iso$results_gene, annotated$glob_iso$GeneExp,"glob_isoseq",diff="yes")
#draw_heatmap_gene_level(loaded$glob_iso$results_gene, annotated$glob_iso$GeneExp,"glob_isoseq",diff="no")


## ----- MAPT transgene -----

pMapt <- list(
  glob_iso = find_mapt_isoseq(paste0(dirnames$glob_root,"/2_post_isoseq3/11_transgene"),phenotype$whole_rTg4510_iso),
  targ_iso = find_mapt_isoseq(paste0(dirnames$targ_iso_root,"/10_characterise/transgene"),phenotype$targeted_rTg4510_iso),
  targ_ont = find_mapt_ont(paste0(dirnames$targ_ont_root,"/0_characterise/transgene"),phenotype$targeted_rTg4510_ont)
)

## ----- compare global WT and TG descriptive data -----

pGlobWTvsTG <- list(
  length = plot_iso_length_mdatasets(group_class.files.diff, "length", "Isoform length (kb)") +
    scale_x_discrete(limits = c("Both","WT","TG")) +
    scale_fill_manual(values = c(label_colour("novel"), label_colour("TG"), label_colour("WT"))) + 
    scale_y_continuous(labels = ks),
  
  exon = plot_iso_length_mdatasets(group_class.files.diff, "exons", "Number of exons") +
    scale_x_discrete(limits = c("Both","WT","TG")) +
    scale_fill_manual(values = c(label_colour("novel"), label_colour("TG"), label_colour("WT")))#, 
  
 # numiso = group_class.files.diff %>% group_by(Dataset, associated_gene) %>% tally() %>%
#    plot_iso_length_mdatasets(., "n", "Number of isoforms per gene") +
#    scale_x_discrete(limits = c("Both","WT","TG")) +
#    scale_fill_manual(values = c(label_colour("novel"), label_colour("TG"), label_colour("WT")))
)

## ----- compare effect size of gene expression (global isoseq data) ------

pGlobIsoVsRna <- recapitulate_gene_level()


## ----- sensitivity -------
ptargetall <- plot_cupcake_collapse_sensitivity(class.files$targ_all,"All 20 target genes")


## ----- genotype differential expression ----

top10Genotype <- TargetedMergedDESeq$waldGenotype %>% arrange(padj_ont) %>% .[1:10,"isoform"]
ptop10Genotype <- lapply(top10Genotype, function(x) time_case_boxplot(TargetedDESeq$ontResTranAnno$wald$norm_counts,x))

## ----- global Iso-Seq data C4b: differential gene and expression analysis  -----

pC4b <- list(
  IsoGeneExp = plot_trans_exp_individual_overtime("7022",GlobalDESeq$resGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "Iso-Seq Gene Expression") +
    theme(legend.position = c(0.15,0.9)) +
    theme(legend.title=element_blank()) + 
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG")),
  
  RNAGeneExp = plot_trans_exp_individual_overtime("7022",GlobalDESeq$RresGeneAnno$lrt$norm_counts,type="gene") + 
    labs(title = "", subtitle = "RNA-Seq Gene Expression", y = "Normalized counts (K)") + scale_y_continuous(labels = ks) +
    theme(legend.position = c(0.15,0.9)) +
    theme(legend.title=element_blank()) + 
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG")),
  
  tracks1 = ggTranPlots(gtf$glob_iso,class.files$glob_iso,
                        isoList = c("PB.7022.9","PB.7022.8","PB.7022.13","PB.7022.1","PB.7022.7","PB.7022.30","PB.7022.28"),
                        colours = c("#F8766D","#7CAE00","black",wes_palette("Zissou1")[4],"gray","gray","gray"), gene = "C4b"),
  
  # PB.7022.9 = red, PB.7022.8 = blue
  IsoTransExp  = plot_transexp_overtime("C4b",GlobalDESeq$resTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.7022.9")) + 
    scale_colour_manual(values = c("#7CAE00","#F8766D")) + theme(legend.position = c(0.2,0.75)) + theme(legend.title=element_blank()) +  
    labs(title = "", subtitle = "Iso-Seq Transcript Expression", y = "Normalized counts") +
    facet_grid(cols = vars(group), labeller = labeller(group = as_labeller(c("CONTROL" = "WT", "CASE" = "TG")))),
  
  # PB.2973.16 = red
  #RNATransExp = plot_transexp_overtime("C4b",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.7022.9")) + 
  #  theme(legend.position = "None") + scale_y_continuous(labels = ks) + 
  #  labs(title = "", subtitle = "RNA-Seq Transcript Expression", y = "Normalised counts (K)"),
  
  # PB.7022.9 = red, PB.7022.13 = black, PB.7022.1 = yellow
  RNATransExpRanked = plot_transexp_overtime("C4b",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.7022.9")) +
    scale_colour_manual(values = c(wes_palette("Zissou1")[4],"black","#F8766D")) + 
    theme(legend.position = c(0.2,0.75)) + theme(legend.title=element_blank()) + scale_y_continuous(labels = ks) +
    labs(title = "", subtitle = "RNA-Seq Transcript Expression", y = "Normalized counts (K)") +
    facet_grid(cols = vars(group), labeller = labeller(group = as_labeller(c("CONTROL" = "WT", "CASE" = "TG"))))
)


plot_transexp_overtime("C4b",GlobalDESeq$RresTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.7022.9"))
plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="toprank",rank=3,isoSpecific=c("PB.14646.39341")) 
plot_transexp_overtime("Snca",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.38419.249")) 

pNovelOntTargetedDiff <- list(
  Clu1 = plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.14646.39352")), 
  Fyn1 = plot_transexp_overtime("Fyn",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.3948.4511")),
  Clu2 = plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.14646.35283")),
  Apoe = plot_transexp_overtime("Apoe",TargetedDESeq$ontResTranAnno$lrt$norm_counts,show="specific",isoSpecific=c("PB.40586.1023")) 
)
pNovelOntTargetedDiff <- lapply(pNovelOntTargetedDiff, 
                                function(x) x + facet_grid(cols = vars(group),labeller=labeller(group=as_labeller(c("CONTROL"="WT","CASE"="TG")))))

FSMIsoforms <- lapply(c("Clu","Fyn","Apoe"), 
                      function(x) class.files$targ_filtered[class.files$targ_filtered$associated_gene == x & class.files$targ_filtered$structural_category == "FSM","isoform"])
names(FSMIsoforms) <- c("Clu","Fyn","Apoe")

pNovelOntTargetedTracks <- list(
  Clu = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
            isoList = c("PB.14646.39352","PB.14646.35283",FSMIsoforms$Clu[1:2]),
            colours = c("#F8766D","#F8766D",rep("#7CAE00",2)), lines = c("#F8766D","#F8766D",rep("#7CAE00",2)), gene = "Clu"),
  Fyn = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
            isoList = c("PB.3948.4511",FSMIsoforms$Fyn[1:2]),
            colours = c("#F8766D",rep("#7CAE00",2)), lines = c("#F8766D",rep("#7CAE00",2)), gene = "Fyn"), 
  Apoe = ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
            isoList = c("PB.40586.1023",FSMIsoforms$Apoe[1:2]),
            colours = c("#F8766D",rep("#7CAE00",2)), lines = c("#F8766D",rep("#7CAE00",2)), gene = "Apoe")
)

plot_transexp_overtime("Mapt",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.8675.37810"))
plot_transexp_overtime("Bin1",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.22007.224"))
plot_transexp_overtime("Snca",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.38419.87"))
plot_transexp_overtime("App",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",isoSpecific=c("PB.19309.7497"))


plot_transexp_overtime("Clu",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",
                       isoSpecific=TargetedDESeq$ontResTranAnno$wald$anno_res[TargetedDESeq$ontResTranAnno$wald$anno_res$associated_gene == "Clu","isoform"])

ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
            isoList = TargetedDESeq$ontResTranAnno$wald$anno_res[TargetedDESeq$ontResTranAnno$wald$anno_res$associated_gene == "Clu","isoform"],
            colours = c("#F8766D",wes_palette("Darjeeling2")[2],"#00BFC4","#7CAE00",wes_palette("GrandBudapest2")[2],wes_palette("Zissou1")[4]),
            lines = c("#F8766D",wes_palette("Darjeeling2")[2],"#00BFC4","gray",wes_palette("GrandBudapest2")[2],wes_palette("Zissou1")[4]),
            gene = "Clu")

## ----- relationship between number of isoforms and other features -----
numIsoRel <- numIso_relationship(classfiles=class.files$targ_filtered,
                    geneExp=TargetedDESeq$ontResGeneAnno$wald$norm_counts_all,
                    transExp=TargetedDESeq$ontResTranAnno$wald$norm_counts_all,
                    refGencode=Targeted$ref_gencode)


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
for(i in 1:length(pIF)){pIF[[2]][[i]]}

pIFClu <- plotIF("Clu",
       ExpInput=Exp$targ_ont$normAll,
       pheno=phenotype$targeted_rTg4510_ont,
       cfiles=class.files$targ_all,
       design="time_series", rank = 5,
       majorIso=NULL)

pIFBin1 <- plotIF("Bin1",
       ExpInput=Exp$targ_ont$normAll,
       pheno=phenotype$targeted_rTg4510_ont,
       cfiles=class.files$targ_all,
       design="time_series", rank = 5,
       majorIso=NULL)


## ---------- Target Genes overview -----------------

class.files$targ_filtered %>% select(isoform, contains("ONT")) %>% mutate(annot_transcript_id = isoform)

gA5A3 <- plot_summarised_AS_events(Merged_gene_class_df, dirnames$targ_anno, Targeted$ref_gencode)[[2]]
gES <- plot_summarised_ES(Targeted$Gene_class, class.files$targ_filtered, dirnames$targ_anno, Targeted$ref_gencode, Targeted$ref_altcon)
gIR <- plot_summarised_IR(class.files$targ_filtered, dirnames$targ_anno, TargetedDESeq$ontResTranAnno$wald$norm_counts)

AppESMax <- ES %>% filter(associated_gene == "App") %>% group_by(transcript_id) %>% tally() %>% filter(n > 15)
Bin1ESMax <- ES %>% filter(associated_gene == "Bin1") %>% group_by(transcript_id) %>% tally() %>% filter(n > 15)

AppESTrack <- ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
            isoList = c(as.character(AppESMax$transcript_id[1:5]),as.character(FSM$App$isoform[6:10])),
            colours = c(rep("#F8766D",5),rep("#7CAE00",5)),
            lines = c(rep("#F8766D",5),rep("#7CAE00",5)), gene = "App") 

Bin1ESTrack <- ggTranPlots(gtf$targ_merged,class.files$targ_filtered,
                          isoList = c(as.character(Bin1ESMax$transcript_id[1:5]),as.character(FSM$Bin1$isoform[1:5])),
                          colours = c(rep("#F8766D",5),rep("#7CAE00",5)),
                          lines = c(rep("#F8766D",5),rep("#7CAE00",5)), gene = "Bin1") 


IR_tab_exonoverlap <- input_FICLE_splicing_results(dirnames$targ_anno, "IntronRetentionExonOverlap")
TardbpIso <- data.frame(
  Isoform = unlist(TardbpIso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Tardbp" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    IR = as.character(IR_tab_exonoverlap[IR_tab_exonoverlap$IRNumExonsOverlaps > 7,"transcript_id"])
  )),
  Category = rep(names(TardbpIso), lengths(TardbpIso))
)
TardbpIso$colour <- c(rep(NA,length(TardbpIso$Category[TardbpIso$Category != "DTE"])))
TardbpIRTrack <- ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
            isoList = c(as.character(TardbpIso$Isoform)),
            selfDf = TardbpIso, gene = "Tardbp")

Cd33Iso <- data.frame(
  Isoform = unlist(Cd33Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Cd33" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    IR = c("PB.41115.1383")
  )),
  Category = rep(names(Cd33Iso), lengths(Cd33Iso))
)
Cd33Iso$colour <- c(rep(NA,length(Cd33Iso$Category[Cd33Iso$Category != "DTE"])))
Cd33IRTrack <- ggTranPlots(gtf$targ_merged, class.files$targ_filtered,
                             isoList = c(as.character(Cd33Iso$Isoform)),
                             selfDf = Cd33Iso, gene = "Cd33")

## ---------- Revision -----------------------------------------------------------------


## ---------- Crytpic exons ----------
InternalNovelExons <- lapply(FICLENE, function(x) x[x$classification == "Internal_NovelExon",])
InternalNovelExons <- lapply(InternalNovelExons, function(x) merge(x, 
                                                                   mouseProtein$cpat[,c("seq_ID","Coding_prob")], 
                                                                   by.x = "transcriptID", by.y = "seq_ID"))
InternalNovelExons <- lapply(InternalNovelExons, function(x) x %>% mutate(coding = ifelse(Coding_prob < 0.44, "noncoding","coding")))
InternalNovelExons <- lapply(InternalNovelExons, function(x) merge(x, 
                                                                   nmd[,c("pb","is_nmd","has_stop_codon")], 
                                                                   by.x = "transcriptID", by.y = "pb",all.x=TRUE)) 

plotCE <- function(gene, type){
  representative <- as.data.frame(gtf$ref_target %>% filter(type == "transcript" & transcript_type == "protein_coding") %>% group_by(gene_name) %>% top_n(1, width))
  allRep <- as.data.frame(gtf$ref_target)
  # manual drop gene
  print(gene)
  df <- InternalNovelExons[[gene]]
  if(gene == "Apoe"){
    # PB.40586.2026 = reference match ("PB.40586.2026")
    df <- df[df$transcriptID %in% c("PB.40586.26305"),]
    ref <-data.frame("PB.40586.2026",NA,NA,NA,NA,NA,"coding","False")
    colnames(ref) <- colnames(df)
    df <- rbind(df, ref)
  }else if(gene == "Bin1"){
    df <- df[df$transcriptID %in% paste0("PB.22007.",c("925","1554","1033","4967","1470","14222","1014")),]
  }else if(gene == "Clu"){
    df <-  df[df$transcriptID %in% paste0("PB.14646.",c("6992")),]
    # PB.40586.2026 = reference match ("PB.14646.139", "PB.14646.483")
    ref <- data.frame("PB.14646.139",NA,NA,NA,NA,NA,"coding","False")
    colnames(ref) <- colnames(df)
    df <- rbind(df, ref)
  }else if(gene == "Snca"){
    df <-  df[df$transcriptID %in% paste0("PB.38419.",c("3145")),]
  }else if(gene == "Trem2"){
    df <- df[df$transcript %in% paste0("PB.20818.", c("80","192","573","1074","493","362","1096")),]
  }else{
    return(NULL)
  }
  transcriptList <- list(
    #Reference = c("ENSMUST00000024791.14","ENSMUST00000113237.3"),
    Reference =  unique(allRep[allRep$gene_name == gene, "transcript_id"]),
    #Reference =  representative[ representative$gene_name == gene, "transcript_id"],
    `not NMD` = as.character(unique(df %>% filter(coding == "coding" & is_nmd == "False") %>% .[,c("transcriptID")])),
    `NMD` = as.character(unique(df %>% filter(coding == "coding" & is_nmd == "True") %>% .[,c("transcriptID")])),
    `non-coding` = as.character(unique(df %>% filter(coding == "noncoding") %>% .[,c("transcriptID")]))
  )
  proteinList <- list(
    #Reference =  unique(allRep[allRep$gene_name == gene, "transcript_id"]),
    `not NMD` =  unique(gtf$ptarg_merged %>% filter(transcript %in% transcriptList$`not NMD`) %>% .[["gene_id"]]),
    `NMD` = unique(gtf$ptarg_merged %>% filter(transcript %in% transcriptList$`NMD`) %>% .[["gene_id"]]),
    `non-coding` = unique(gtf$ptarg_merged %>% filter(transcript %in% transcriptList$`non-coding`) %>% .[["gene_id"]]))
  
  if(type == "transcript"){
    p <-  generalggTranPlots(transcriptList, gtf$targ_merged, class.files$targ_filtered, gene) 
  }else if(type == "protein"){
    p <- generalggTranPlots(proteinList, gtf$targ_merged, class.files$targ_filtered, gene) 
  }else{
    bothList <- c(transcriptList, proteinList)
    p <- generalggTranPlots(bothList, gtf$targ_merged, class.files$targ_filtered, gene) 
  }
  
  return(p)
}


pInternalNovelExonsTranscript <- lapply(names(InternalNovelExons), function(x) plotCE(x,"transcript"))
names(pInternalNovelExonsTranscript) <- names(InternalNovelExons)

pInternalNovelExonsProtein <- lapply(names(InternalNovelExons), function(x) plotCE(x,"protein"))
names(pInternalNovelExonsProtein) <- names(InternalNovelExons)

pInternalNovelExonsProteinBoth <- lapply(names(InternalNovelExons), function(x) plotCE(x,"proteinTranscript"))
names(pInternalNovelExonsProteinBoth) <- names(InternalNovelExons)

CE <- c(paste0("PB.22007.",c("925","1554","1033","4967","1470","14222","1014")),
        "PB.40586.26305",
        paste0("PB.14646.",c("6992")),
        paste0("PB.38419.",c("3145")),
        paste0("PB.20818.", c("80","192","573","1074","493","362","1096")))
plot_expression_summed(CE, "Cryptic exons", AgeDiv = TRUE)
plot_expression_summed(CE)

plot_grid(pInternalNovelExonsProteinBoth$Apoe,pInternalNovelExonsProteinBoth$Bin1, pInternalNovelExonsProteinBoth$Clu, 
          pInternalNovelExonsProteinBoth$Trem2, pInternalNovelExonsProtein$Snca, ncol = 2)
  
## ---------- IR ----------

IRList <- lapply(TargetGene,function(x) IR[IR$associated_gene == x,"transcript_id"])
names(IRList) <- TargetGene
IRp <- lapply(IRList, function(x) plot_expression_summed(x, AgeDiv = TRUE))
IRGenop <- lapply(IRList, function(x) plot_expression_summed(x, AgeDiv = FALSE))


IR <- list(
  Trem2Exp = IRp$Trem2 + labs(subtitle = "Trem2: IR transcripts"),
  Cd33Exp = IRp$Cd33 + labs(subtitle = "Cd33: IR transcripts"),
  Trem2Tracks = ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
                             isoList = c(as.character(IRList$Trem2),RefIsoforms$Trem2[1:2]),
                             colours = c(rep(wes_palette("Royal1")[2],5),rep("#0C0C78",10)), 
                             lines = c(rep(wes_palette("Royal1")[2],5),rep("#0C0C78",10)),
                             gene = "Trem2",simple=TRUE),
  Cd33Tracks = ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
                           isoList = c(as.character(IRList$Cd33),RefIsoforms$Cd33[1:2]),
                           colours = c(rep(wes_palette("Royal1")[2],9),rep("#0C0C78",10)), 
                           lines = c(rep(wes_palette("Royal1")[2],9),rep("#0C0C78",10)),
                           gene = "Cd33",simple=TRUE)
  
)


## ---------- NMD ----------

# list of NMD transcripts
NMD <- lapply(TargetGene,function(x) class.files$protein_filtered_final[class.files$protein_filtered_final$is_nmd == "TRUE" & class.files$protein_filtered_final$associated_gene == x,"corrected_acc"])
names(NMD) <- TargetGene

# plot summed expression of NMD by gene (genotype * age)
NMDp <- lapply(NMD, function(x) plot_expression_summed(x, AgeDiv = TRUE))
for(i in 1:length(TargetGene)){NMDp[[i]] <- NMDp[[i]] + labs(subtitle = TargetGene[[i]])}
# plot summed expression of NMD by gene (genotype )
NMDGenop <- lapply(NMD, function(x) plot_expression_summed(x, AgeDiv = FALSE))
for(i in 1:length(TargetGene)){NMDGenop[[i]] <- NMDGenop[[i]] + labs(subtitle = TargetGene[[i]])}
names(NMDp) <- TargetGene
names(NMDGenop) <- TargetGene

NMD_plots <- list(
  Trem2Exp = NMDp$Trem2 + labs(subtitle = "Trem2: NMD transcripts"),
  Cd33Exp = NMDp$Cd33 + labs(subtitle = "Cd33: NMD transcripts"),
  Trem2Tracks = visualise_ORFs_NMD(NMD$Trem2, "Trem2"),
  Cd33Tracks = visualise_ORFs_NMD(NMD$Cd33, "Cd33")
)


## ---------- Proteogenomics ----------

pGeneralProtein <- plot_protein_general(class.files$ptarg_filtered, class.files$protein_filtered_final, ORFCorrectedCollapsedID)
pTrem2SameORF <- visualise_ORFs(refgtf=gtf$ref_target, tgtf=gtf$targ_merged,pgtf=gtf$ptarg_merged,
               tclassfiles=class.files$targ_filtered, pclassfiles=class.files$ptarg_filtered,gene="Trem2", transcript="PB.20818.54", cpat=mouseProtein$cpat, species="mouse")

pTrem2UniqueORF <- visualise_ORFs(refgtf=gtf$ref_target, tgtf=gtf$targ_merged,pgtf=gtf$ptarg_merged,
               tclassfiles=class.files$targ_filtered, pclassfiles=class.files$ptarg_filtered,gene="Trem2",cpat=mouseProtein$cpat, species="mouse")


## ---------- Output -----------------

pdf(paste0(output_dir,"/SuppFigures.pdf"), width = 10, height = 12)
plot_grid(pMapt$glob_iso[[2]],get_legend(pMapt$glob_iso[[1]]),pMapt$targ_iso[[2]],pMapt$targ_ont, labels = c("A","","B","C"))
plot_grid(plotlist = pGlobWTvsTG, ncol = 2, labels = c("A","B"))
plot_grid(plot_grid(pC4b$IsoGeneExp,pC4b$RNAGeneExp,nrow=1,labels=c("A","B")),
          pC4b$tracks1,
          plot_grid(pC4b$IsoTransExp,pC4b$RNATransExpRanked,labels=c("D","E")),
          ncol=1, labels = c("","C",""))
plot_grid(pGlobIsoVsRna)
plot_grid(ptargetall[[2]])
plot_grid(plotlist = numIsoRel, labels = c("A","B","C","D"))
plot_grid(gA5A3,gES[[4]],AppESTrack,Bin1ESTrack,labels=c("A","B","C","D"))
plot_grid(Cd33IRTrack,plot_grid(gIR[[1]],gIR[[2]]+theme(legend.position = c(0.8,0.7)),gIR[[3]],nrow=1,ncol=3, labels = c("B","C","D")),TardbpIRTrack,ncol=1,labels = c("A","","E"), rel_heights = c(0.2,0.4,0.5))
top = plot_grid(plot_grid(pNovelOntTargetedDiff[[3]], pNovelOntTargetedDiff[[1]], ncol = 1),pNovelOntTargetedTracks$Clu, labels = c("A","B"))
mid = plot_grid(pNovelOntTargetedDiff[[2]], pNovelOntTargetedTracks$Fyn, labels = c("C","D"))
bot = plot_grid(pNovelOntTargetedDiff[[4]], pNovelOntTargetedTracks$Apoe, labels = c("E","F"))
plot_grid(top,mid,bot,ncol=1,rel_heights = c(0.5,0.25,0.25))
legend <- get_legend(ptop10Genotype[[1]] + theme(legend.position = "top"))
ptop10Genotype <- lapply(ptop10Genotype, function(p) p + theme(legend.position = "none") + labs(x=NULL,y=NULL))
ptop10Genotype <- plot_grid(plotlist = ptop10Genotype,legend = legend,labels=c("","A","B","C","D","E","F","G","H","I","J"))
y.grob <- textGrob("Normalized counts", gp=gpar(fontsize=18), rot=90)
x.grob <- textGrob("Genotype", gp=gpar(fontsize=18))
grid.arrange(arrangeGrob(ptop10Genotype, left = y.grob, bottom = x.grob))
plot_grid(pIF$ontNorm$Fus[[1]],pIF$ontNorm$Bin1[[1]],pIF$ontNorm$Trpa1[[1]],
          pIF$ontNorm$Bin1[[2]] + theme(legend.position = "None"),
          pIF$ontNorm$Vgf[[2]] + theme(legend.position = "None"),
          pIF$ontNorm$Clu[[2]] + theme(legend.position = "None"), labels = c("A","B","C","D","E","F"))
dev.off()

pdf(paste0(output_dir,"/SuppFigures2.pdf"), width = 20, height = 12)
plot_grid(pIFClu[[1]],pIFClu[[2]],pIFBin1[[1]],pIFBin1[[2]])
dev.off()

# redundant 
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

pdf(paste0(dirnames$targ_output,"/SuppProteomics.pdf"), width = 22, height = 13)
A <- plot_grid(pGeneralProtein[[1]], pGeneralProtein[[2]], pTrem2SameORF, ncol = 1, labels = c("A","B","C"), scale = 0.95, 
          rel_heights = c(0.3,0.45,0.35))
B <- plot_grid(pTrem2UniqueORF, labels = c("D"), scale = 0.95)
plot_grid(A,B, rel_widths = c(0.35,0.65))
dev.off()

pdf(paste0(dirnames$targ_output,"/IR_NMD.pdf"), width = 22, height = 13)
plot_grid(IR$Trem2Tracks, IR$Trem2Exp, IR$Cd33Tracks, IR$Cd33Exp, labels = c("A","B","C","D"))
plot_grid(NMD_plots$Trem2Tracks, NMD_plots$Trem2Exp, NMD_plots$Cd33Tracks, NMD_plots$Cd33Exp, labels = c("A","B","C","D"))
dev.off()
