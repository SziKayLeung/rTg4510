class.files$ptarg_filtered <- SQANTI_class_preparation(paste0(dirnames$targ_root, "/2_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered_pCollapsed.txt"),"nstandard")
annopResTran <- readRDS(paste0(dirnames$targ_output, "/DESeq2ProteinLevel.RDS"))

protein$t2p.collapse <- protein$t2p.collapse %>% filter(!orf_calling_confidence == "Low Quality ORF") %>% 
  filter(gene%in%TargetGene) %>% 
  mutate(numtxCollapsed = count.fields(textConnection(as.character(pb_accs)), sep = "|"))

## same protein sequence as LR.Trem2.54
Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    RNA_Transcript = class.files$ptarg_filtered %>% filter(corrected_acc == "PB.20818.54") %>% .[["isoform"]],
    RNA_Protein = paste0(class.files$ptarg_filtered %>% filter(corrected_acc == "PB.20818.54") %>% .[["isoform"]],"_ORF_1")
  )),
  Category = rep(names(Trem2Iso), lengths(Trem2Iso))
)
Trem2Iso$colour <- c(rep(NA,length(Trem2Iso$Category[Trem2Iso$Category != "DTE"])))

p1 <- ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(Trem2Iso$Isoform)),
            selfDf = Trem2Iso, gene = "Trem2")

## different protein sequence annotated to Trem2
Trem2pID <- unique(class.files$ptarg_filtered[class.files$ptarg_filtered$associated_gene == "Trem2","corrected_acc"])
unique(gtf$ptarg_merged %>% filter(transcript %in% Trem2pID) %>% .[["gene_id"]])


/lustre/projects/Research_Project-MRC148213/sl693/rTg4510_FICLE/FICLE/TargetGenes/Trem2


proteinclass.files <- read.table("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510/G_Merged_Targeted/4_proteogenomics/7_classified_protein/all_iso_ont.sqanti_protein_classification.tsv", sep = "\t", header = T)
nmd <- read.table("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510/G_Merged_Targeted/4_proteogenomics/7_classified_protein/all_iso_ont.classification_filtered.tsv", sep = "\t", header = T)
nmd %>% filter(pr_gene == "Trem2")
idx <- match(nmd$pb,class.files$ptarg_filtered$base_acc)
nmd <- transform(nmd, corrected_acc = ifelse(!is.na(idx), as.character(class.files$ptarg_filtered$corrected_acc[idx]), pb))
proteinclass.files 

nmdTrem2 <- nmd[nmd$is_nmd == "True" & nmd$pr_gene == "Trem2","corrected_acc"]
nonnmDTrem2 <-  nmd[nmd$is_nmd == "False" & nmd$pr_gene == "Trem2","corrected_acc"]
Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    RNA_TranscriptNMD = as.character(nmdTrem2),
    RNA_Protein = unique(gtf$ptarg_merged %>% filter(transcript %in% nmdTrem2) %>% .[["gene_id"]])
  )),
  Category = rep(names(Trem2Iso), lengths(Trem2Iso))
)
Trem2Iso$colour <- c(rep(NA,length(Trem2Iso$Category[Trem2Iso$Category != "DTE"])))

NMD <- ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(Trem2Iso$Isoform)),
            selfDf = Trem2Iso, gene = "Trem2")


Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    RNA_TranscriptNMD = as.character(nonnmDTrem2),
    RNA_Protein = unique(gtf$ptarg_merged %>% filter(transcript %in% nonnmDTrem2) %>% .[["gene_id"]])
  )),
  Category = rep(names(Trem2Iso), lengths(Trem2Iso))
)
Trem2Iso$colour <- c(rep(NA,length(Trem2Iso$Category[Trem2Iso$Category != "DTE"])))

nonNMD <- ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
                   isoList = c(as.character(Trem2Iso$Isoform)),
                   selfDf = Trem2Iso, gene = "Trem2")


Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    NMD = as.character(nmdTrem2),
    not_NMD = as.character(nonnmDTrem2)
  )),
  Category = rep(names(Trem2Iso), lengths(Trem2Iso))
)
Trem2Iso$colour <- c(rep(NA,length(Trem2Iso$Category[Trem2Iso$Category != "DTE"])))
ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(Trem2Iso$Isoform)),
            selfDf = Trem2Iso, gene = "Trem2")


plot_trem2 <- function(transcript){
  p <- plot_transexp_overtime("Trem2",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=3,
                         isoSpecific=c(transcript),setorder=c("CONTROL","CASE")) + 
    scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00"))
  return(p)
}
plot_trem2("PB.20818.260")
plot_trem2("PB.20818.274")
plot_trem2("PB.20818.379")
plot_trem2("PB.20818.382")
plot_trem2("PB.20818.547")

## WT vs TG 
TargetedDESeq$ontResTranAnno$wald$norm_counts <- TargetedDESeq$ontResTranAnno$wald$norm_counts %>% mutate(sampleID = word(sample,c(2),sep=fixed("_")))

WTTGTargetedCounts <- merge(TargetedDESeq$ontResTranAnno$wald$norm_counts %>% filter(sampleID %in% targetedTG) %>% group_by(isoform) %>% summarise(TGSum = sum(normalised_counts)),
      TargetedDESeq$ontResTranAnno$wald$norm_counts %>% filter(sampleID %in% targetedWT) %>% group_by(isoform) %>% summarise(WTSum = sum(normalised_counts))
)
WTTGTargetedCounts[WTTGTargetedCounts$TGSum == 0,]
WTTGTargetedCounts[WTTGTargetedCounts$WTSum == 0,]


### 
ES <- read.csv("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510_FICLE/FICLE/TargetGenes/Trem2/Stats/Trem2_Exonskipping_generaltab.csv") 
ES4 <- ES %>% filter(Gencode_4 == "Yes" & Gencode_5 == "No")
ES5 <- ES %>% filter(Gencode_4 == "No" & Gencode_5 == "Yes")
ES4and5 <- ES %>% filter(Gencode_4 == "Yes" & Gencode_5 == "Yes")
Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    ES_exon4 = as.character(unique(ES4$X)),
    ES_exon4 = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(ES4$X),) %>% .[["gene_id"]]),
    ES_exon5 = as.character(unique(ES5$X)),
    ES_exon5 = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(ES5$X),) %>% .[["gene_id"]]),
    ES_exon4and5 = as.character(unique(ES4and5$X)),
    ES_exon4and5 = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(ES4and5$X),) %>% .[["gene_id"]])
  )),
  Category = rep(names(Trem2Iso), lengths(Trem2Iso))
)
Trem2Iso$colour <- c(rep(NA,length(Trem2Iso$Category[Trem2Iso$Category != "DTE"])))
ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(Trem2Iso$Isoform)),
            selfDf = Trem2Iso, gene = "Trem2")

## cryptic exon
# in frame or nmd and low-quality ORFs
NE <- read.csv("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510_FICLE/FICLE/TargetGenes/Trem2/Stats/Trem2_NE.csv") 
Trem2Iso <- data.frame(
  Isoform = unlist(Trem2Iso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    Reference = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(unique(paste0("PB.20818.", c("54","62")))),) %>% .[["gene_id"]]),
    NE = as.character(paste0("PB.20818.", c("80","192","573","1074"))),
    NE = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(unique(paste0("PB.20818.", c("80","192","573","1074")))),) %>% .[["gene_id"]]),
    NE2 = as.character(paste0("PB.20818.", c("493","362","1096"))),
    NE2 = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(unique(paste0("PB.20818.", c("493","362","1096")))),) %>% .[["gene_id"]])
  )),
  Category = rep(names(Trem2Iso), lengths(Trem2Iso))
)
Trem2Iso$colour <- c(rep(NA,length(Trem2Iso$Category[Trem2Iso$Category != "DTE"])))
ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(Trem2Iso$Isoform)),
            selfDf = Trem2Iso, gene = "Trem2")


### Mapt
# 2N4R vs 2N3R

Maptprotein <- unique(class.files$ptarg_filtered[class.files$ptarg_filtered$associated_gene == "Mapt","corrected_acc"])
MaptES <- read.csv("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510_FICLE/Mapt/Stats/Mapt_general_exon_level.csv")
Mapt0N4R <- MaptES %>% filter(Gencode_3 == "Yes" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "No" & Gencode_14 == "No" & Gencode_15 == "No")
Mapt1N4R <- MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "No" & Gencode_14 == "No" & Gencode_15 == "No")
Mapt2N4R <- MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "No" & Gencode_11 == "No" & Gencode_12 == "No" & Gencode_14 == "No" & Gencode_15 == "No")
Mapt2N3R <- MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "No" & Gencode_11 == "No" & Gencode_12 == "Yes" & Gencode_14 == "No" & Gencode_15 == "No")
Mapt1N3R <- MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "Yes" & Gencode_14 == "No" & Gencode_15 == "No")
Mapt0N3R <- MaptES %>% filter(Gencode_3 == "Yes" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "Yes" & Gencode_14 == "No" & Gencode_15 == "No")

MaptIso <- data.frame(
  Isoform = unlist(MaptIso <- list(
    Reference = c("ENSMUST00000106992.9","ENSMUST00000100347.10"),
    Mapt0N4R = intersect(Maptprotein,Mapt0N4R$isoform),
    Mapt1N4R = intersect(Maptprotein,Mapt1N4R$isoform),
    Mapt2N4R = intersect(Maptprotein,Mapt2N4R$isoform),
    Mapt0N3R = intersect(Maptprotein,Mapt0N3R$isoform),
    Mapt1N3R = intersect(Maptprotein,Mapt1N3R$isoform),
    Mapt2N3R = intersect(Maptprotein,Mapt2N3R$isoform)
    #MaptProtein = unique(gtf$ptarg_merged %>% filter(transcript %in% intersect(Maptprotein,Mapt2N4R$isoform),) %>% .[["gene_id"]])
  )),
  Category = rep(names(MaptIso), lengths(MaptIso))
)
MaptIso$colour <- c(rep(NA,length(MaptIso$Category[MaptIso$Category != "DTE"])))
ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(MaptIso$Isoform)),
            selfDf = MaptIso, gene = "Mapt")

MaptIso <- data.frame(
  Isoform = unlist(MaptIso <- list(
    Reference = c("ENSMUST00000106992.9","ENSMUST00000100347.10"),
    Mapt = intersect(Maptprotein,Mapt0N4R$isoform),
    MaptProtein = unique(gtf$ptarg_merged %>% filter(transcript %in% intersect(Maptprotein,Mapt0N4R$isoform),) %>% .[["gene_id"]])
  )),
  Category = rep(names(MaptIso), lengths(MaptIso))
)
MaptIso$colour <- c(rep(NA,length(MaptIso$Category[MaptIso$Category != "DTE"])))
ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(MaptIso$Isoform)),
            selfDf = MaptIso, gene = "Mapt")

MaptIso <- data.frame(
  Isoform = unlist(MaptIso <- list(
    Reference = "ENSMUST00000106992.9",
    Mapt = intersect(Maptprotein,Mapt2N3R$isoform),
    MaptProtein = unique(gtf$ptarg_merged %>% filter(transcript %in% intersect(Maptprotein,Mapt2N3R$isoform),) %>% .[["gene_id"]]),
    TGONly = "PB.8675.38019"
  )),
  Category = rep(names(MaptIso), lengths(MaptIso))
)
MaptIso$colour <- c(rep(NA,length(MaptIso$Category[MaptIso$Category != "DTE"])))
ggTranPlots(inputgtf=gtf$targ_merged,classfiles=class.files$targ_filtered,
            isoList = c(as.character(MaptIso$Isoform)),
            selfDf = MaptIso, gene = "Mapt")


plot_transexp_overtime("Mapt",TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=3,
                       isoSpecific=c("PB.8675.557"),setorder=c("CONTROL","CASE")) + 
  scale_colour_manual(values = c(wes_palette("Darjeeling1")[3],"#00BFC4","#7CAE00"))

MaptSpec <- rbind(data.frame(Mapt = "2N4R", isoform = Mapt2N4R$isoform),
                  data.frame(Mapt = "0N4R", isoform = Mapt0N4R$isoform),
                  data.frame(Mapt = "1N4R", isoform = Mapt1N4R$isoform))
Mapt3RSpec <- rbind(
                  data.frame(Mapt = "0N3R", isoform = Mapt0N3R$isoform),
                  data.frame(Mapt = "1N3R", isoform = Mapt1N3R$isoform))

TargetedDESeq$ontResTranAnno$wald$norm_counts %>% 
  filter(isoform %in% MaptSpec$isoform) %>%
  left_join(., MaptSpec) %>% 
  group_by(isoform, group, time, Mapt) %>% 
  summarise(meanCounts = mean(normalised_counts)) %>% 
  ungroup()%>% 
  ggplot(., aes(x = group, y = log10(meanCounts), colour = as.factor(time))) + geom_boxplot() + 
  #geom_point(aes(fill = as.factor(time)), size = 1, shape = 21, position = position_jitterdodge())  + 
  facet_grid(~Mapt)
t.test(meanCounts ~ group, dat)

TargetedDESeq$isoResTranAnno$wald$norm_counts %>% filter(isoform %in% Mapt2N4R$isoform) %>% 
  group_by(isoform, group, time) %>% 
  summarise(meanCounts = median(normalised_counts)) %>% 
  ungroup() %>%
  ggplot(., aes(x = group, y = log10(meanCounts))) + geom_boxplot()




class.files$targ_filtered[class.files$targ_filtered$isoform %in% WTTGTargetedCounts[WTTGTargetedCounts$WTSum == 0,"isoform"],]


ratioplots <- function(){
  dat1 <- TargetedDESeq$ontResTranAnno$wald$norm_counts_all %>% filter(isoform %in% MaptSpec$isoform) %>% 
    mutate(normalised_counts = normalised_counts + 1) %>%
    group_by(sample) %>% 
    summarise(meanCounts = mean(normalised_counts))
  
  dat2 <- TargetedDESeq$ontResTranAnno$wald$norm_counts_all %>% filter(isoform %in% Mapt3RSpec$isoform) %>% 
    mutate(normalised_counts = normalised_counts + 1) %>%
    group_by(sample) %>% 
    summarise(meanCountsnot = mean(normalised_counts))
  
  dat <- merge(dat1,dat2) %>% mutate(Ratio = meanCounts/meanCountsnot) %>% mutate(sampleID = word(sample,c(2),sep=fixed("_"))) %>% 
    merge(., phenotype$targ_ont, by.x = "sampleID", by.y = "sample") 
  
  p <- ggplot(dat, aes(x = group, y = Ratio)) + geom_boxplot()
  
  print(dat)
  
  return(p)
}

ratioplots()

