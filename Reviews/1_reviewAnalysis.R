#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## simple analysis rewiew 
##    find reference genes
##    correlation
##
## --------------------------------

library("dplyr")
library("stringr")
library("cowplot")
library("data.table")
library("ggplot2")
library(ggrepel)
suppressMessages(library("gridExtra"))
suppressMessages(library("grid"))
source("/lustre/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Reviews/rTg4510_config.R")

options(scipen = 999)

## ---------------------------  find reference genes
housekeepingGenes = c("Las1l","Rrp1","Gusb","Polr2b","Cyc1","Tbp","Xpnpep1","Gapdh","Actb","Rpl13a","Sdha")
class.files$glob_iso[class.files$glob_iso$associated_gene %in% housekeepingGenes,] %>% 
  group_by(associated_gene) %>% tally()

class.files$targ_all[class.files$targ_all$associated_gene %in% housekeepingGenes,] %>% 
  group_by(associated_gene) %>% tally()


# gene expression plots 
hkdf <- class.files$glob_iso[class.files$glob_iso$associated_gene %in% housekeepingGenes,] %>% arrange(associated_gene)
HKGeneExpPlots <- list()
PbGeneId = unique(word(hkdf$isoform,c(2),sep=fixed(".")))
for(i in PbGeneId){
  HKGeneExpPlots[[i]] <- plot_trans_exp_individual_overtime(i,GlobalDESeq$resGeneAnno$wald$norm_counts,type="gene") + labs(x = NULL, y = NULL)
}
names(HKGeneExpPlots) <- unique(hkdf$associated_gene)

# most abundant top-ranked expression plot
HKTransExpPlots <- list()
for(i in unique(hkdf$associated_gene)){
  HKTransExpPlots [[i]] <- plot_transexp_overtime(i,GlobalDESeq$resTranAnno$wald$norm_counts_all,show="toprank",rank=3,setorder=c("CONTROL","CASE")) + 
    labs(x = NULL, y = NULL) + theme(legend.position = "None")
}

HKGeneExpPlotsList <- plot_grid(plotlist = HKGeneExpPlots)
HKTransExpPlotsList <- plot_grid(plotlist = HKTransExpPlots)
y.grob <- textGrob("Normalized counts", gp=gpar(fontsize=18), rot=90)
x.grob <- textGrob("Age (months)", gp=gpar(fontsize=18))
pHKGeneExpPlotsList  <- grid.arrange(arrangeGrob(HKGeneExpPlotsList, left = y.grob, bottom = x.grob))
pHKTransExpPlotsList  <- grid.arrange(arrangeGrob(HKTransExpPlotsList, left = y.grob, bottom = x.grob))



## ---------------------------  correlation 

## correlate the number of novel isoforms to number of all reads
# determine the total number of reads per sample (target + offtargets)
# rawExp = demux_fl_count.csv
class.files$targ_offtargets <- merge(class.files$targ_offtargets, rawExp$targ_ont_all, by = "isoform", all.x = T)
allReadCountSample <- colSums(class.files$targ_offtargets %>% select(contains(c("ONT", "Iso.Seq")))) %>% reshape2::melt(., value.name = "totalReads") %>% tibble::rownames_to_column(., var = "sample")

# fiilter only novel isoforms in the target gene panel (before expression filter)
novelIsoforms <- class.files$targ_all %>% filter(associated_transcript == "novel") %>% select(contains(c("ONT", "Iso.Seq")))

# sum the number of novel isoforms (frequency) for each sample; i.e. count the number of occurence where FL read count != 0 
novelIsoformsSample <- apply(novelIsoforms, 2, function(c)sum(c!=0)) %>% reshape2::melt(., value.name = "novelNum") %>% tibble::rownames_to_column(., var = "sample")

# count the total number of target reads per sample (before expression filter)
allTargetReadCountSample <- colSums(class.files$targ_all %>% select(contains(c("ONT", "Iso.Seq")))) %>% reshape2::melt(., value.name = "totalReads") %>% tibble::rownames_to_column(., var = "sample")

# create df for correlationl merging the number of novel isoforms per sample and the number of total reads per sample 
# remove ONT_sum_FL and Iso-Seq_sum_FL columns
# separate platforms 
corrNovelNumReads <- merge(novelIsoformsSample,allReadCountSample) %>% filter(!sample %in% c("ONT_sum_FL","Iso.Seq_sum_FL")) %>% 
  mutate(platform = word(sample,c(1),sep=fixed("_")), sampleID = word(sample,c(2),sep=fixed("_")))
PBcorrNovelNumReads <- corrNovelNumReads[corrNovelNumReads$platform == "Iso.Seq",]
ONTcorrNovelNumReads <- corrNovelNumReads[corrNovelNumReads$platform == "ONT",]

# correlation test (Pearson's)
cor.test(PBcorrNovelNumReads$novelNum,PBcorrNovelNumReads$totalReads)
cor.test(ONTcorrNovelNumReads$novelNum,ONTcorrNovelNumReads$totalReads)

pCorr1 <- ggplot(corrNovelNumReads, aes(x = novelNum, y = totalReads, colour = platform)) + geom_point() + 
  mytheme + labs(x = "Number of novel isoforms", y = "Number of reads") + 
  scale_colour_discrete(name = "Platform") +
  theme(legend.position = "top") +
  geom_smooth(method = "lm", se = FALSE)

# take into consideration that the samples were batched; determining which sample was in which batch
# if comparing the number of novel isoforms per library (i.e. per batch) using t-test
batch1 <- c("K19", "K23", "K21", "K18", "K20", "K17")
batch2 <- c("S19", "K24", "L22", "M21", "O18", "O23", "O22", "P19", "T20")
batch3 <- c("Q20", "Q21", "S18", "S23", "Q18", "Q17", "L18", "Q23", "T18")
batchedSamples <- data.frame(sampleID = c(batch1,batch2,batch3), batch = c(rep("1", length(batch1)),rep("2", length(batch2)),rep("3", length(batch3))))
corrNovelNumReads <- merge(corrNovelNumReads,batchedSamples,by = "sampleID")

# only ONT reads
CorrNovelNumReads <- corrNovelNumReads %>% mutate(ratio = novelNum/totalReads) 
pCorr2 <- CorrNovelNumReads %>% filter(batch != c(1)) %>% 
  ggplot(., aes(x = batch, y = ratio)) + geom_boxplot() + geom_point() +
  mytheme + labs(x = "Sequencing library", y = "Ratio of novel isoforms: total reads") + facet_grid(~platform) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "white", fill = NA))
t.test(data = CorrNovelNumReads %>% filter(platform == "ONT"), ratio ~ batch)
t.test(data = CorrNovelNumReads %>% filter(platform == "Iso.Seq" & batch != "1"), ratio ~ batch)

# whole transcriptome 
allReadCountSample <- colSums(class.files$glob_iso %>% select(contains("FL."))) %>% reshape2::melt(., value.name = "totalReads") %>% tibble::rownames_to_column(., var = "sample")
novelIsoforms <- class.files$glob_iso %>% filter(associated_transcript == "novel") %>% select(contains(c("FL.")))
novelIsoformsSample <- apply(novelIsoforms, 2, function(c)sum(c!=0)) %>% reshape2::melt(., value.name = "novelNum") %>% tibble::rownames_to_column(., var = "sample")
corrNovelNumReads <- merge(novelIsoformsSample,allReadCountSample)
cor.test(corrNovelNumReads$novelNum, corrNovelNumReads$totalReads)
pCorr3 <- ggplot(corrNovelNumReads, aes(x = novelNum, y = totalReads)) + geom_point(colour = "#F8766D") + 
  mytheme + labs(x = "Number of novel isoforms", y = "Number of reads") + 
  theme(legend.position = "top") 

# gene expression 
corrNovelIsoformGeneExpression <- function(classfiles, normCounts, platform){
  
  if(!platform %in% c("ONT","Iso.Seq")){
    message("Platform argument: ONT, Iso.Seq")
  }
  
  # aggregate the counts by associated gene 
  meanGeneExpression <- aggregate(normalised_counts ~ associated_gene, data = normCounts, FUN = mean)
  colnames(meanGeneExpression) <- c("associated_gene", "mean_gene_counts")
  meanGeneExpression <<- meanGeneExpression
  
  # number of novel isoforms per sample in dataset
  novelIsoform <- classfiles %>%
    # novel isoform
    filter(associated_transcript == "novel") %>%
    # select columns that are FL read counts
    select(contains(c("ONT", "Iso.Seq", "associated_gene"))) %>%
    group_by(associated_gene) %>%
    # for each sample (i.e ONT_S19), count the number of isoforms with reads > 0
    summarise_all(~ sum(. != 0)) %>% 
    select(-contains("sum_FL")) %>% 
    reshape2::melt(., variable.name = "sample", value.name = "num_novel_iso") 
  novelIsoform <<- novelIsoform
  
  # the average number of isoforms for each gene across all samples
  meanNumNovel <- novelIsoform %>% filter(grepl(platform,sample)) %>% group_by(associated_gene) %>% summarise(mean = mean(num_novel_iso))
  meanGeneNumNovel <- merge(meanGeneExpression,meanNumNovel, by = "associated_gene")
  cortest <- cor.test(meanGeneNumNovel$mean_gene_counts, meanGeneNumNovel$mean)
  print(cortest)
  corr.value <- cor(meanGeneNumNovel$mean_gene_counts, meanGeneNumNovel$mean)
  
  corr <- grobTree(textGrob(paste("r = ", round(corr.value, 2)), 
                            x = 0.05, y = 0.80, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 14, fontface = "italic",family="CM Roman")))
  
  
  p <- ggplot(meanGeneNumNovel, aes(y = mean_gene_counts, x = mean, label = as.factor(associated_gene))) + geom_point() + 
    annotation_custom(corr) +
    mytheme + labs(y = "Mean gene expression", x = "Mean number of novel isoforms") + 
    geom_text_repel()#+
    #geom_smooth(method=lm, colour = "gray", se = FALSE)
  
  return(p)
}

pCorrGeneONT <- corrNovelIsoformGeneExpression(class.files$targ_all, TargetedDESeq$ontResGeneAnno$waldgenotype$norm_counts,"ONT") 
pCorrGeneIso <- corrNovelIsoformGeneExpression(class.files$targ_all, TargetedDESeq$isoResGeneAnno$wald$norm_counts,"Iso.Seq")



## ---------------------------  output (pdf)

#pdf(paste0(dirnames$targ_output,"/housekeepingGenes.pdf"), width = 14, height = 17)
pdf("housekeepingGenes.pdf", width = 12, height = 14)
plot_grid(pHKGeneExpPlotsList,pHKTransExpPlotsList,nrow=2, labels = c("A","B"), scale = 0.95)
dev.off()

pdf("Correlation.pdf", width = 10, height = 10)
plot_grid(plot_grid(pCorr3, NULL,rel_widths = c(0.6,0.4), labels = c("A","")), 
          plot_grid(pCorr1,pCorr2, rel_widths = c(0.6,0.4), labels = c("B","C")), nrow = 2, scale = 0.95)
dev.off()
plot_grid(pCorrGeneONT, pCorrGeneIso, labels = c("A","B"))
