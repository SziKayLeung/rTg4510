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
suppressMessages(library("gridExtra"))
suppressMessages(library("grid"))

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

class.files$targ_offtargets <- merge(class.files$targ_offtargets, rawExp$targ_ont_all, by = "isoform", all.x = T)
novelIsoforms <- class.files$targ_all %>% filter(associated_transcript == "novel") %>% select(contains(c("ONT", "Iso.Seq")))
novelIsoformsSample <- apply(novelIsoforms, 2, function(c)sum(c!=0)) %>% reshape2::melt(., value.name = "novelNum") %>% tibble::rownames_to_column(., var = "sample")
allReadCountSample <- colSums(Filter(is.numeric, class.files$targ_all %>% select(contains(c("ONT", "Iso.Seq"))))) %>% reshape2::melt(., value.name = "totalReads") %>% tibble::rownames_to_column(., var = "sample")

corrNovelNumReads <- merge(novelIsoformsSample,allReadCountSample) %>% filter(!sample %in% c("ONT_sum_FL","Iso.Seq_sum_FL")) %>% 
  mutate(platform = word(sample,c(1),sep=fixed("_")), sampleID = word(sample,c(2),sep=fixed("_")))
PBcorrNovelNumReads <- corrNovelNumReads[corrNovelNumReads$platform == "Iso.Seq",]
ONTcorrNovelNumReads <- corrNovelNumReads[corrNovelNumReads$platform == "ONT",]
cor.test(PBcorrNovelNumReads$novelNum,PBcorrNovelNumReads$totalReads)
cor.test(ONTcorrNovelNumReads$novelNum,ONTcorrNovelNumReads$totalReads)

library("ggplot2")
ggplot(corrNovelNumReads, aes(x = novelNum, y = totalReads, colour = platform)) + geom_point() + mytheme


## library
batch1 <- c("K19", "K23", "K21", "K18", "K20", "K17")
batch2 <- c("S19", "K24", "L22", "M21", "O18", "O23", "O22", "P19", "T20")
batch3 <- c("Q20", "Q21", "S18", "S23", "Q18", "Q17", "L18", "Q23", "T18")
batchedSamples <- data.frame(sampleID = c(batch1,batch2,batch3), batch = c(rep("1", length(batch1)),rep("2", length(batch2)),rep("3", length(batch3))))
corrNovelNumReads <- merge(corrNovelNumReads,batchedSamples,by = "sampleID")

ontCorrNovelNumReads <- corrNovelNumReads[corrNovelNumReads$platform == "ONT",] %>% mutate(ratio = novelNum/totalReads) 
t.test(data = ontCorrNovelNumReads, ratio ~ batch)

batchCorrNovelNumReads <- corrNovelNumReads %>% 
  group_by(platform, batch) %>% 
  summarise(across(c(novelNum, totalReads), list(sum = sum)))

ggplot(ontCorrNovelNumReads, aes(x = batch, y = ratio)) + geom_boxplot()


# gene expression 
meanGeneExpression <- aggregate(normalised_counts ~ associated_gene, data = TargetedDESeq$ontResGeneAnno$waldgenotype$norm_counts, FUN = mean)
colnames(meanGeneExpression) <- c("associated_gene", "mean_gene_counts")

dat <- class.files$targ_all %>%
  filter(associated_transcript == "novel") %>%
  select(contains(c("ONT", "Iso.Seq", "associated_gene")))

result <- dat %>%
  group_by(associated_gene) %>%
  summarise_all(~ sum(. != 0))

meanNumNovel <- reshape2::melt(result , variable.name = "sample", value.name = "num_novel_iso") %>% 
  group_by(associated_gene) %>% summarise(mean = mean(num_novel_iso))
  
merge(meanGeneExpression,meanNumNovel, by = "associated_gene") %>% ggplot(., aes(x = mean_gene_counts, y = mean)) + geom_point()


## ---------------------------  output (pdf)

#pdf(paste0(dirnames$targ_output,"/housekeepingGenes.pdf"), width = 14, height = 17)
pdf("housekeepingGenes.pdf", width = 12, height = 14)
plot_grid(pHKGeneExpPlotsList,pHKTransExpPlotsList,nrow=2, labels = c("A","B"), scale = 0.95)
dev.off()