#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Mapt 4R vs 3R
## Trem2 proteogenomics characterisation
## --------------------------------

suppressMessages(library("cowplot"))

## ---------- config file -----------------

source("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/Reviews/rTg4510_config.R")

generalggTranPlots <- function(isolist, inputgtf, classfiles, gene, cpat = NULL, species = NULL){
  IsoDf <- data.frame(
    Isoform = unlist(IsoDf <- isolist),
    Category = rep(names(IsoDf), lengths(IsoDf))
  )
  IsoDf$colour <- c(rep(NA,length(IsoDf$Category[IsoDf$Category != "DTE"])))
  p <- ggTranPlots(inputgtf=inputgtf,classfiles=classfiles,
                   isoList = c(as.character(IsoDf$Isoform)),
                   selfDf = IsoDf, gene = gene, inputCpat = cpat, cpatSpecies = species)
  return(p)
}

label_group <- function(genotype){
  if(genotype %in% c("TG","Case","CASE")){group = "TG"}else{
    if(genotype %in% c("WT","Control","CONTROL")){group = "WT"}}
  return(group)
}

## ---------- input from 2_proteogenomics.R -----------------

# original transcript classification file with additional column of the number of collapsed transcripts by ORF 
# col: corrected_acc = representative isoform after collapsing by ORF 
class.files$ptarg_filtered <- SQANTI_class_preparation(paste0(dirnames$targ_root, "/2_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered_pCollapsed.txt"),"nstandard")
protein$t2p.collapse.refined <- read.table(paste0(dirnames$protein,"6_refined_database/all_iso_ont_orf_refined_collapsed.tsv"))
annopResTran <- readRDS(paste0(dirnames$targ_output, "/DESeq2ProteinLevel.RDS"))

dirnames$protein <- paste0(dirnames$targ_root,"/4_proteogenomics/7_classified_protein/")

class.files$protein <- read.table(paste0(dirnames$protein,"all_iso_ont.sqanti_protein_classification.tsv"), sep = "\t", header = T)
nmd <- read.table(paste0(dirnames$protein, "all_iso_ont.classification_filtered.tsv"), sep = "\t", header = T)
idx <- match(nmd$pb,class.files$ptarg_filtered$base_acc)
nmd <- transform(nmd, corrected_acc = ifelse(!is.na(idx), as.character(class.files$ptarg_filtered$corrected_acc[idx]), NA))

# basic summary stats
message("Number of RNA Transcripts:", nrow(class.files$targ_filtered))
message("Number of coding RNA isoforms:", length(unique(protein$t2p.collapse.refined$pb_acc)))
message("Number of coding RNA isoforms:", length(unique(protein$t2p.collapse.refined$corrected_acc)))
length(unique(class.files$ptarg_filtered$corrected_acc)) # -1 to not include "NA" 

NumTranscriptIsoform <- merge(class.files$ptarg_filtered %>% group_by(associated_gene) %>% summarize(RNATranscript = n()),
      class.files$ptarg_filtered %>% select(corrected_acc, associated_gene) %>% filter(!is.na(corrected_acc)) %>% 
        distinct() %>% group_by(associated_gene) %>% summarize(RNAIsoform = n()), by = "associated_gene") 

ggplot(NumTranscriptIsoform, aes(x = RNAIsoform, y = RNATranscript, label = associated_gene)) + geom_point() + geom_text_repel() +
  geom_abline(intercept = 0, linetype = "dashed") + mytheme

dat <- class.files$ptarg_filtered %>% filter(!is.na(corrected_acc)) %>% filter(!isoform %in% class.files$ptarg_filtered$corrected_acc) %>% 
  group_by(structural_category, subcategory) %>% tally() 

  ggplot(., aes(y = n, x = subcategory)) + geom_bar(stat = "identity") + facet_grid(rows = vars(structural_category), scales = "free") + coord_flip()

dat <- class.files$ptarg_filtered %>% filter(!is.na(corrected_acc)) %>% filter(isoform %in% class.files$ptarg_filtered$corrected_acc) %>% 
  group_by(structural_category, subcategory) %>% tally() 

class.files$ptarg_filtered %>% group_by(corrected_acc) %>% tally() %>% group_by(n) %>% tally() 
ggplot(aes(x = n, y = nn)) + geom_bar(stat = "identity")

table(class.files$protein$pr_splice_cat)
nrow(class.files$protein)
3496/nrow(class.files$protein)


## ---------- WT vs TG -----------------

## identify transcripts that are unique to WT and TG mice using normalized expression counts (ONT)
# summarised the normalised counts by group and isoform
TargetedDESeq$ontResTranAnno$wald$norm_counts <- TargetedDESeq$ontResTranAnno$wald$norm_counts %>% mutate(sampleID = word(sample,c(2),sep=fixed("_")))
WTTGTargetedCounts <- merge(TargetedDESeq$ontResTranAnno$wald$norm_counts %>% filter(sampleID %in% targetedTG) %>% group_by(isoform) %>% summarise(TGSum = sum(normalised_counts)),
                            TargetedDESeq$ontResTranAnno$wald$norm_counts %>% filter(sampleID %in% targetedWT) %>% group_by(isoform) %>% summarise(WTSum = sum(normalised_counts))
)

# subset classification file based on the unique transcript
class.files$targ_filtered_WTUnique <- class.files$targ_filtered[class.files$targ_filtered$isoform %in% WTTGTargetedCounts[WTTGTargetedCounts$WTSum == 0,"isoform"],]
class.files$targ_filtered_TGUnique <- class.files$targ_filtered[class.files$targ_filtered$isoform %in% WTTGTargetedCounts[WTTGTargetedCounts$TGSum == 0,"isoform"],]


## ---------- Mapt -----------------

# identification of Mapt Transcripts by exon skipping of original exon 2 and 3, 9, 10, 11, 12 
# FICLE Gencode_13 = alternative first exon
# FICLE Gencode_3, Gencode_4 = origianl exon 2, 3; Gencode_11,_12,_14,15 = original exons 9 - 12
# Yes = Exon skipped; No = Exon present and not skipped
MaptRTranscripts <- list(
  Ref = list(isoform = c("ENSMUST00000106992.9","ENSMUST00000100347.10")),
  Mapt0N3R = MaptES %>% filter(Gencode_3 == "Yes" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "Yes" & Gencode_14 == "No" & Gencode_15 == "No"),
  Mapt0N4R = MaptES %>% filter(Gencode_3 == "Yes" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "No" & Gencode_14 == "No" & Gencode_15 == "No"),
  Mapt1N3R = MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "Yes" & Gencode_14 == "No" & Gencode_15 == "No"),
  Mapt1N4R = MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "Yes" & Gencode_11 == "No" & Gencode_12 == "No" & Gencode_14 == "No" & Gencode_15 == "No"),
  Mapt2N3R = MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "No" & Gencode_11 == "No" & Gencode_12 == "Yes" & Gencode_14 == "No" & Gencode_15 == "No"),
  Mapt2N4R = MaptES %>% filter(Gencode_3 == "No" & Gencode_4 == "No" & Gencode_11 == "No" & Gencode_12 == "No" & Gencode_14 == "No" & Gencode_15 == "No")
)
MaptRTranscripts <- lapply(MaptRTranscripts, function(x) as.character(x[["isoform"]]))
names(MaptRTranscripts) <- str_remove(names(MaptRTranscripts),"Mapt")
# create a dataframe for downstream subsetting <MaptType><isoform>
MaptRTranscriptsdf <- do.call(rbind, MaptRTranscripts) %>% reshape2::melt() %>% select(Var1, value) %>% `colnames<-`(c("MaptType", "isoform")) %>%
  mutate(R = ifelse(grepl("3R", MaptType), "3R","4R")) 

# ggtranscript of Mapt transcripts with highlights of N and R regions
pMaptRTranscripts <- generalggTranPlots(MaptRTranscripts, gtf$targ_merged, class.files$targ_filtered, "Mapt") + 
  annotate("rect", xmin = c(104286000, 104309000), xmax = c(104291000, 104325000), ymin = -Inf, ymax = Inf, alpha = .1, fill = c("green"))

# ggtranscript of Mapt transcripts with respective ORF
MaptRProtein  <- lapply(MaptRTranscripts, function(x) unique(gtf$ptarg_merged %>% filter(transcript %in% x,) %>% .[["gene_id"]]))
MaptRProtein$Ref <-  c("ENSMUST00000106992.9","ENSMUST00000100347.10")
pMaptRProtein <- generalggTranPlots(MaptRProtein, gtf$targ_merged, class.files$targ_filtered, "Mapt", cpat = protein$cpat, species = "mouse") +
  annotate("rect", xmin = c(104286000, 104309000), xmax = c(104291000, 104325000), ymin = -Inf, ymax = Inf, alpha = .1, fill = c("green"))

# expression of Mapt transcripts (ONT normalised counts)
pMaptRTrascriptExp <- TargetedDESeq$ontResTranAnno$wald$norm_counts_all %>% 
  merge(., MaptRTranscriptsdf, by = "isoform", all.y = TRUE) %>% 
  filter(MaptType != "Reference") %>%
  mutate(group = factor(ifelse(group == "CONTROL","WT","TG"), levels = c("WT","TG"))) %>%
  mutate(age = as.factor(time)) %>%  
  mutate(MaptType = str_remove(MaptType,"Mapt")) %>%
  group_by(isoform, group, age, MaptType, R) %>% 
  summarise(meanCounts = mean(normalised_counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = group, y = log10(meanCounts), colour = age)) + geom_boxplot() + 
  geom_point(aes(fill = age), size = 1, shape = 21, position = position_jitterdodge())  + 
  facet_nested(~R + MaptType, nest_line = element_line(linetype = 2)) +
  theme(strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "grey")) + mytheme +
  labs(x = "Genotype") + theme(legend.position = "top") + 
  scale_colour_manual(values = c("black","#CFCFCF","#777777","red"), name = "Age (months)") +
  scale_fill_manual(values = c("black","#CFCFCF","#777777","red"), name = "Age (months)") 

# ratio of Mapt 4R vs 3R transcript by mean expression
meanMaptExp <- function(normcounts, isoList){
  dat <- normcounts%>% filter(isoform %in% isoList) %>% 
    # + 1 as 0 for many of the isoforms
    mutate(normalised_counts = normalised_counts + 1) %>%
    # take the average of all the normalised counts of a subset of isoforms for each sample
    group_by(sample) %>% 
    summarise(meanCounts = mean(normalised_counts))
  return(dat)
}

# ratio of mean 4R isoform expression/ mean 3R isoform expresssion for each sample
dat <- merge(meanMaptExp(TargetedDESeq$ontResTranAnno$wald$norm_counts_all,MaptRTranscriptsdf[MaptRTranscriptsdf$R == "4R","isoform"]),
             meanMaptExp(TargetedDESeq$ontResTranAnno$wald$norm_counts_all,MaptRTranscriptsdf[MaptRTranscriptsdf$R == "3R","isoform"]), by = "sample") %>% 
  mutate(Ratio = meanCounts.x/meanCounts.y) %>% mutate(sampleID = word(sample,c(2),sep=fixed("_"))) %>% 
  merge(., phenotype$targ_ont, by.x = "sampleID", by.y = "sample") %>%
  mutate(group = factor(ifelse(group == "CONTROL","WT","TG"), levels = c("WT","TG")))
# plot of ratio across age and genotype
pMaptRTrascriptRatio1 <- ggplot(dat, aes(x = as.factor(time), y = Ratio, colour = group)) + geom_point() + mytheme +
  labs(x = "Age (months)", y = "Ratio") + theme(legend.position = "right") +
  scale_colour_manual(values = c(label_colour("WT"),"red"), name = "Genotype") + 
  stat_summary(data=dat, aes(x=as.factor(time), y=Ratio, group=group), fun ="mean", geom="line", linetype = "dotted") 
# plot of ratio across genotype 
pMaptRTrascriptRatio2 <- ggplot(dat, aes(x = group, y = Ratio)) + geom_boxplot() + mytheme +
  labs(x = "Genotype", y = "Ratio") +
  geom_point()  


## ---------- Trem2 -----------------

# visaulisation of Trem2 transcripts with same ORF LR.Trem2.54
Trem2ProteinSameORF <- list(
  Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
  `RNA Transcript` = class.files$ptarg_filtered %>% filter(corrected_acc == "PB.20818.54") %>% .[["isoform"]],
  `RNA Isoform` = paste0(class.files$ptarg_filtered %>% filter(corrected_acc == "PB.20818.54") %>% .[["isoform"]],"_ORF_1")
)
pTrem2ProteinSameORF <- generalggTranPlots(Trem2ProteinSameORF, gtf$targ_merged, class.files$targ_filtered, "Trem2") 

# different protein sequence annotated to Trem2
Trem2pID <- unique(class.files$ptarg_filtered[class.files$ptarg_filtered$associated_gene == "Trem2","corrected_acc"])
pTremFinal <- list(
  Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
  `RNA Transcript` = Trem2pID
) %>% generalggTranPlots(., gtf$targ_merged, class.files$targ_filtered, "Trem2")
pTremFinalProtein <- list(
  Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
  `RNA Isoform` = unique(gtf$ptarg_merged %>% filter(transcript %in% Trem2pID) %>% .[["gene_id"]])
) %>% generalggTranPlots(., gtf$targ_merged, class.files$targ_filtered, "Trem2", protein$cpat, "mouse")

# exon skipping 
Trem2ESTranscripts <- list(
  Ref = list(X = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"])),
  `E3-E4+` = Trem2ES %>% filter(Gencode_4 == "Yes" & Gencode_5 == "No"),
  `E3+E4-` = Trem2ES %>% filter(Gencode_4 == "No" & Gencode_5 == "Yes"),
  `E3-E4-` = Trem2ES %>% filter(Gencode_4 == "Yes" & Gencode_5 == "Yes")
)
Trem2ESTranscripts <- lapply(Trem2ESTranscripts, function(x) as.character(x[["X"]]))
pTrem2ESTranscripts <- generalggTranPlots(Trem2ESTranscripts, gtf$targ_merged, class.files$targ_filtered, "Trem2")  +
    annotate("rect", xmin = c(48351000), xmax = c(48352000), ymin = -Inf, ymax = Inf, alpha = .1, fill = c("green"))
Trem2ESProtein  <- lapply(Trem2ESTranscripts , function(x) unique(gtf$ptarg_merged %>% filter(transcript %in% x,) %>% .[["gene_id"]]))
Trem2ESProtein$Ref <-  unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"])
pTrem2ESProtein <- generalggTranPlots(Trem2ESProtein, gtf$targ_merged, class.files$targ_filtered, "Trem2", cpat = protein$cpat, species = "mouse") +
  annotate("rect", xmin = c(48351000), xmax = c(48352000), ymin = -Inf, ymax = Inf, alpha = .1, fill = c("green"))

# cryptic exon
# in frame or nmd and low-quality ORFs
Trem2CE <- list(
  Reference = c("ENSMUST00000024791.14","ENSMUST00000113237.3"),
  `NE - not NMD` = as.character(paste0("PB.20818.", c("80","192","573","1074"))),
  `NE - NMD` = as.character(paste0("PB.20818.", c("493","362","1096")))
) %>% generalggTranPlots(., gtf$targ_merged, class.files$targ_filtered, "Trem2") + 
  annotate("rect", xmin = c(48347000), xmax = c(48347400), ymin = -Inf, ymax = Inf, alpha = .1, fill = c("green"))

Trem2CEProtein <- list(
  Reference = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(unique(paste0("PB.20818.", c("54","62")))),) %>% .[["gene_id"]]),
  `NE - not NMD` = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(unique(paste0("PB.20818.", c("80","192","573","1074")))),) %>% .[["gene_id"]]),
  `NE - NMD` = unique(gtf$ptarg_merged %>% filter(transcript %in% unique(unique(paste0("PB.20818.", c("493","362","1096")))),) %>% .[["gene_id"]])
)%>% generalggTranPlots(., gtf$targ_merged, class.files$targ_filtered,"Trem2",cpat = protein$cpat, species = "mouse") + 
  annotate("rect", xmin = c(48347000), xmax = c(48347400), ymin = -Inf, ymax = Inf, alpha = .1, fill = c("green"))


# nonsense-mediated decay 
nmdList <- list(
  Trem2True = nmd[nmd$is_nmd == "True" & nmd$pr_gene == "Trem2","corrected_acc"],
  Trem2False =  nmd[nmd$is_nmd == "False" & nmd$pr_gene == "Trem2","corrected_acc"]
)
Trem2nmd <- list(
  Ref = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
  `Transcript NMD` = as.character(nmdList$Trem2True),
  `Isoform NMD` = unique(gtf$ptarg_merged %>% filter(transcript %in% nmdList$Trem2True) %>% .[["gene_id"]])
)
Trem2notnmd <- list(
  Ref = unique(gtf$ref_target[gtf$ref_target$gene_name == "Trem2" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
  `Transcript Not NMD` = as.character(nmdList$Trem2False),
  `Protein Not NMD` = unique(gtf$ptarg_merged %>% filter(transcript %in% nmdList$Trem2False) %>% .[["gene_id"]])
)
generalggTranPlots(Trem2nmd, gtf$targ_merged, class.files$targ_filtered, "Trem2")  
generalggTranPlots(Trem2notnmd, gtf$targ_merged, class.files$targ_filtered, "Trem2")  


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


## ---------- cryptic exons -----------------

InternalNovelExons <- lapply(FICLENE, function(x) x[x$classification == "Internal_NovelExon",])
InternalNovelExons <- lapply(InternalNovelExons, function(x) merge(x, 
                                                                   protein$cpat[,c("pb_acc","coding_score","orf_calling_confidence")], 
                                                                   by.x = "transcriptID", by.y = "pb_acc"))
InternalNovelExons <- lapply(InternalNovelExons, function(x) x %>% mutate(coding = ifelse(coding_score < 0.44, "noncoding","coding")))
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
    ref <-data.frame("PB.40586.2026",NA,NA,NA,NA,NA,"coding","False","True")
    colnames(ref) <- colnames(df)
    df <- rbind(df, ref)
  }else if(gene == "Bin1"){
    df <- df[df$transcriptID %in% paste0("PB.22007.",c("925","1554","1033","4967","1470","14222","1014")),]
  }else if(gene == "Clu"){
    df <-  df[df$transcriptID %in% paste0("PB.14646.",c("6992")),]
    # PB.40586.2026 = reference match ("PB.14646.139", "PB.14646.483")
    ref <- data.frame("PB.14646.139",NA,NA,NA,NA,NA,"coding","False","True")
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

pdf(paste0(dirnames$targ_output,"/CrypticExonsApoe.pdf"), width = 10, height = 8)
pInternalNovelExonsProteinBoth$Apoe
dev.off()

pdf(paste0(dirnames$targ_output,"/CrypticExonsClu.pdf"), width = 10, height = 8)
pInternalNovelExonsProteinBoth$Clu
dev.off()

pdf(paste0(dirnames$targ_output,"/CrypticExonsTrem2.pdf"), width = 10, height = 8)
pInternalNovelExonsProteinBoth$Trem2
dev.off()

pdf(paste0(dirnames$targ_output,"/CrypticExonsBin1.pdf"), width = 10, height = 8)
pInternalNovelExonsProteinBoth$Bin1
dev.off()

pdf(paste0(dirnames$targ_output,"/CrypticExonsSnca.pdf"), width = 10, height = 4)
pInternalNovelExonsProteinBoth$Snca
dev.off()


## ---------- output (pdf) -----------------

pdf(paste0(dirnames$targ_output,"/MaptIsoforms.pdf"), width = 14, height = 17)
#pdf("MaptIsoforms.pdf", width = 14, height = 17)
plot_grid(plot_grid(pMaptRTranscripts,pMaptRProtein, labels = c("i","ii")),
          plot_grid(pMaptRTrascriptExp),
          plot_grid(pMaptRTrascriptRatio2,pMaptRTrascriptRatio1,rel_widths = c(0.4,0.6),labels = c("i","ii")),
          ncol=1, rel_heights = c(0.5,0.3,0.2), labels = c("A","B","C"), scale = 0.95)
dev.off()


pdf("Trem2IsoformsIdentialOrf.pdf", width = 10, height = 10)
plot_grid(pTrem2ProteinSameORF)
dev.off()
pdf("Trem2IsoformsUniqueOrf.pdf", width = 10, height = 20)
plot_grid(pTremFinal,pTremFinalProtein, nrow=1)
dev.off()
pdf("Trem2IsoformsESOrf.pdf", width = 10, height = 10)
plot_grid(pTrem2ESTranscripts, pTrem2ESProtein, nrow = 1) 
dev.off()
pdf("Trem2IsoformsCEOrf.pdf", width = 10, height = 10)
plot_grid(Trem2CE,Trem2CEProtein,nrow=1)
dev.off()