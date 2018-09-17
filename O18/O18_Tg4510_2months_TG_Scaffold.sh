---
title: "Comparison with RNASeqData"
output: html_document
Date: 17 September 2018
---
AD mouse model, age: Tg4510, 2months (Baseline)
Phenotype: TG 
Brain region: Entorinhal Cortex
Sample: O18

Objective: To tabulate the number of full-length reads obtained per gene from Isoseq and order genes from high to low
for comparison with RNAseq data for exact sample 

Rationale: To evaluate whether Isoseq output comparable to RNAseq output 

Analysis: 
1) Downloaded raw subread.bam file from Sequel output 
2) CCS and Isoseq3 command line (Lima, Cluster, Polish)
3) Mapped to mouse genome using GMAP 
4) Tofu Cupcake 
5) Sqanti for isoform characterisation  

```{r setup, include=FALSE}
# knitr::opts_knit$set(root.dir ="//isad.isadroot.ex.ac.uk/UOE/User/My Documents/PacBio/AD_Mouse_Tg4510/2 months/O18")
knitr::opts_knit$set(root.dir ="~")
```
## Input RNA-seq Count
```{r}
# Input: Isabel's Tg4510 short-read ordered by raw counts 
#RNAseq <- read.csv("//isad.isadroot.ex.ac.uk/UOE/User/My Documents/PacBio/AD_Mouse_Tg4510/countdata_ordered_Tg4510.csv", row.names = 1, header = TRUE)
RNAseq <- read.csv("/home/sLeung/IsoSeq3_Tg4510/countdata_ordered_Tg4510.csv", row.names = 1, header = TRUE)

head(RNAseq)

# only include RNAseq raw counts from Sample O18(Tg4510) for direct comparisons
# row.names of gene name to copy across from RNASeq as lost during splicing 
RNAseq_O18 <- data.frame(RNAseq[,"O18"],row.names = (rownames(RNAseq)))
colnames(RNAseq_O18)[1] = "RNASeq Raw Counts"

head(RNAseq_O18)

```
## Input PacBio Count
Input PacBio FL data fromm cupcake with isoform identification from Sqanti
```{r}
# Input: Isoseq O18 processed from Count.O18.Rmd  
# O18_PacBio_complete <- read.csv("//isad.isadroot.ex.ac.uk/UOE/User/My Documents/PacBio/AD_Mouse_Tg4510/2 months/O18/Isoseq3/O18_ordered_unique.csv", header = TRUE)[,-1]
O18_PacBio_complete <- read.csv("/mnt/data1/Szi/O18/Comparisons_RNAseq/O18_ordered_unique.csv", header = TRUE)[,-1]

# For Gene column with novel genes, replace NA with associated gene 
# first to convert both columns into character 
O18_PacBio_complete$mgi_symbol <- as.character(O18_PacBio_complete$mgi_symbol)
O18_PacBio_complete$associated_gene <- as.character(O18_PacBio_complete$associated_gene)
# if gene column is NA, replace from text in associated gene column in respective row
O18_PacBio_complete$mgi_symbol <- ifelse(is.na(O18_PacBio_complete$mgi_symbol) == TRUE, O18_PacBio_complete$associated_gene, O18_PacBio_complete$mgi_symbol)


# Only take 3 columns forwarded for comparison with RNA-seq data: 1) gene names 2) PacBio ID for identification of different isoforms per gene 3) Number of FL counts 
PacBio_O18 <- O18_PacBio_complete[,c("mgi_symbol","isoform","count_fl")]
# As example: different isoforms for same gene 
PacBio_O18[which(PacBio_O18 == "Arf3"),]
```
## Merge RNAseq and PacBio Count
```{r}
# Merge RNAseq and PacBio count by RNAseq's rownames (column = 0) and PacBio's column name "mgi"
# keep all rows even if NA values in respective dataframes i.e For gene A, where RNASeq rawcount = X but PacBio FL count = NA 
MCount <- merge(RNAseq_O18,PacBio_O18,by.x = 0,by.y = "mgi_symbol",all = TRUE)
colnames(MCount) <- c("Genes","RNASeq_Raw_Counts","PacBio_ID","PacBio_FL_Counts")

# checking for correct merging
head(MCount[which(MCount$Genes == "Arf3"),])

# remove rows(genes) that are not detected in both PacBio and RNAseq data 
# only interested in genes that are both or either detected in technology for comparison of expression/threshold for detection
# splice out the rows that are 0 counts for both technology and name as none
none <-  MCount[which(MCount$RNASeq_Raw_Counts == 0 &  is.na(MCount$PacBio_FL_Count)),]
# number of genes before cleaning 
dim(MCount)[1]
# remove rows in MCount dataframe if the gene column matches the gene column in none column i.e. filter out the genes that are not detected
MCount <- MCount[!MCount$Genes %in% none$Genes, , drop = FALSE]
dim(MCount)[1]

# to check the genes where there are RNAseq counts but no PacBio counts 
head(MCount[which(is.na(MCount$PacBio_FL_Counts)),])
```
## Sum of FL Count of all isoforms per gene 
As in example above, PacBio FL count presented per isoform rather than per gene. However, RNA-seq count from Isabel presented for counts per gene (and not isoform). Therefore for direct comparison of count per gene, summed PacBio FL counts of isofors per gene. 
Did also consider choosing only the isoform with the highest number of FL counts, yet biased results especially given if many isoforms with similar or slgihtly smaller number of FL-counts. Assuming RNA-seq captures expression of all RNA transcripts irrespective of isoforms, decided to also compare with total number of FL counts
```{r message=FALSE}
# PacBio Output with one gene but multiple isoforms, thus complicates correlation with RNA-seq gene expression data 
# sum of FLcounts of all isoforms across same gene 
#  %>% group_by(Genes) %>% = group/splice MCount dataframe by genes
# In the groups, keep RNAseq Counts and PacBio ID columns as they are 
# sum(PacBio_FL_Counts)
library(dplyr)
MCount_all <- MCount %>% group_by(Genes) %>%
   summarize(RNASeq_Raw_Counts = RNASeq_Raw_Counts[1], PacBio_ID = PacBio_ID[1], 
             PacBio_FL_Counts = sum(PacBio_FL_Counts))
# validation
MCount_all[which(MCount_all$Genes == "Arf3"),]
# write.csv(MCount_all,file = '//isad.isadroot.ex.ac.uk/UOE/User/My Documents/PacBio/AD_Mouse_Tg4510/2 months/O18/MergeCount_allIsoforms.csv')
write.csv(MCount_all,file = '/mnt/data1/Szi/O18/Comparisons_RNAseq/MergeCount_allIsoforms.csv')

```
## Normalise Read Counts for correlation 
```{r}
# Check normality of counts - not normalised 
hist(MCount_all$PacBio_FL_Counts)
hist(MCount_all$RNASeq_Raw_Counts)
# Try normalising counts with different methods 
hist(log2(MCount_all$RNASeq_Raw_Counts))
hist(log(MCount_all$PacBio_FL_Counts))
hist(sqrt(MCount_all$PacBio_FL_Counts))

# Normalise counts using log2 and cbind 
# Convert level of Counts from character to numeric for log2
# Replace all 0 counts with NA as otherwise generate infinity from log2(0)
MCount_all$RNASeq_Raw_Counts <- as.numeric(as.character(MCount_all$RNASeq_Raw_Counts))
MCount_all[MCount_all == 0] <- NA
MCount_all_norm <- cbind(MCount_all, log2(MCount_all$RNASeq_Raw_Counts),log2(MCount_all$PacBio_FL_Counts))

colnames(MCount_all_norm) <- c("Genes","RNASeq_Raw_Counts","PacBio_ID","PacBio_FL_Counts","log2RNaSeq_Raw_Counts",
                                    "log2PacBio_FL_Counts")
head(MCount_all_norm)
                                    
```
## Correlation 
To check correlation of RNAseq raw counts with PacBio FL counts, hypothesise high gene expression (Rawcounts) from RNAseq to also give similarly high expression of FL-counts, even if PacBio semi-quantiative. 
For information on missing values: http://bwlewis.github.io/covar/missing.html
```{r}
# complete.obs to not include rows that have missing observations/counts in either technology 
cor(MCount_all_norm$log2RNaSeq_Raw_Counts,MCount_all_norm$log2PacBio_FL_Counts,method = "pearson", use = "complete.obs")
```
## Correlation Plot of normalised data
```{r message = FALSE}
# Before plotting correlation, want to detect the genes that are detected by RNAseq and not by PacBio and vice-versa 
# ensure working with dataframe 
MCount_all_norm <- data.frame(MCount_all_norm)
# recoverted NA values back to 0, otherwise not detect the genes that are not detected in either technology 
MCount_all_norm[is.na(MCount_all_norm)] <- 0
# Create new column as detection and as default fill with "Both" i.e. detected by both technology 
MCount_all_norm$Detection ="Both"
# For rows with 0 RNASeq raw counts, populate detection column with "noRSeq"
MCount_all_norm$Detection[MCount_all_norm$log2RNaSeq_Raw_Counts == "0"] = "noRSeq"
# For rows with 0 PacBio raw counts, populate detection column with "noPB"
MCount_all_norm$Detection[MCount_all_norm$log2PacBio_FL_Counts == "0"] <- "noPB"
# For rows with 0 PacBio raw counts and 0 , populate detection column with "noPB", should be 0 as removed earlier before
MCount_all_norm$Detection[MCount_all_norm$PacBio_FL_Counts == "0" & MCount_all_norm$RNASeq_Raw_Counts == "0"] <- "none"

dim(MCount_all_norm[which(MCount_all_norm$Detection == "noPB"),])
dim(MCount_all_norm[which(MCount_all_norm$Detection == "Both"),])
dim(MCount_all_norm[which(MCount_all_norm$Detection == "noRSeq"),])
dim(MCount_all_norm)[1]

# Correlation plot
# if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library("ggpubr")
Corplot_norm <- ggscatter(MCount_all_norm, y = "log2RNaSeq_Raw_Counts", x = "log2PacBio_FL_Counts", 
          color = "Detection",palette = c("black", "red", "blue","yellow"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "log2RNaSeq_Raw_Counts", xlab = "log2PacBio_FL_Counts")

Corplot_norm
```
# Correlation plot of raw data
```{r message = FALSE}
# open directory with saved file of Mcount_all, rather than running through code again
# MCount_all <- read.csv("//isad.isadroot.ex.ac.uk/UOE/User/My Documents/PacBio/AD_Mouse_Tg4510/2 months/O18/MergeCount_allIsoforms.csv", header = TRUE)
# Before plotting correlation, want to detect the genes that are detected by RNAseq and not by PacBio and vice-versa 
# ensure working with dataframe 
library("data.table")
MCount_all <- data.frame(MCount_all)
MCount_all$Genes <- as.character(MCount_all$Genes)
MCount_all$PacBio_ID <- as.character(MCount_all$PacBio_ID)
# recoverted NA values back to 0, otherwise not detect the genes that are not detected in either technology 
MCount_all[is.na(MCount_all)] <- 0
# Create new column as detection and as default fill with "Both" i.e. detected by both technology 
MCount_all$Detection ="Both"
# For rows with 0 RNASeq raw counts, populate detection column with "noRSeq"
MCount_all$Detection[MCount_all$RNASeq_Raw_Counts == "0" & !MCount_all$Genes %like% "novel" ] = "noRSeq"
# For rows with 0 PacBio raw counts, populate detection column with "noPB"
MCount_all$Detection[MCount_all$PacBio_FL_Counts == "0"] <- "noPB"
# For rows with 0 PacBio raw counts and 0 , populate detection column with "noPB", should be 0 as removed earlier before
MCount_all$Detection[MCount_all$PacBio_FL_Counts == "0" & MCount_all$RNASeq_Raw_Counts == "0"] <- "none"
# For rows with novel genes, populate detection column with "novel"
# use datatable function %like% as gene names are novel.x, therefore only partial selection
MCount_all$Detection[MCount_all$Genes %like% "novel"] <- "novel PB"

dim(MCount_all[which(MCount_all$Detection == "noPB"),])[1]
dim(MCount_all[which(MCount_all$Detection == "Both"),])[1]
dim(MCount_all[which(MCount_all$Detection == "noRSeq"),])[1]
dim(MCount_all)[1]

# Correlation plot
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library("ggpubr")
Corplot <- ggscatter(MCount_all, y = "RNASeq_Raw_Counts", x = "PacBio_FL_Counts", 
          color = "Detection",
          palette = c("black", "red"," blue","green","yellow"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.use = "complete.obs",
          ylab = "RNASeq_Raw_Counts", xlab = "PacBio_FL_Counts")
Corplot
```
## Histograms of RNAseq genes detected by PacBio vs undetected
Of raw counts i.e not logged
```{r message=FALSE}
# subset data that is not detected by PacBio by noPB and data detected by both
noPB <- MCount_all[MCount_all$Detection == "noPB",]
Both <- MCount_all[MCount_all$Detection == "Both",]
noRseq <- MCount_all[MCount_all$Detection == "noRSeq",]
# very long tail, therfore removed RNAseq counts of 50
hist(noPB$RNASeq_Raw_Counts,breaks = 50)
noPB_50 <- noPB[which(noPB$RNASeq_Raw_Counts > 50),]
hist(noPB_50$RNASeq_Raw_Counts,breaks = 100)
hist(Both$RNASeq_Raw_Counts)

noPB_50_5000 <- noPB[which(noPB$RNASeq_Raw_Counts > 50 & noPB$RNASeq_Raw_Counts < 5000),]
hist(noPB_50_5000$RNASeq_Raw_Counts,breaks = 100)

hist(Both$RNASeq_Raw_Counts)
Both_50_5000 <- Both[which(Both$RNASeq_Raw_Counts > 50 & Both$RNASeq_Raw_Counts < 5000),]
hist(Both_50_5000$RNASeq_Raw_Counts,breaks = 100)

hist(noRseq$PacBio_FL_Counts)
noRseq_0_300 <- Both[which(noRseq$PacBio_FL_Counts > 50 & noRseq$PacBio_FL_Counts < 300),]
hist(noRseq_0_300$PacBio_FL_Counts,breaks = 50)

#komogorpf-smirnof test
```

```{r}
# HeatMap of normalised data
# To gather data in a format that can be mapped i.e from wide to long 
library("tidyr")
# second column is called Technology which essentially has whether the counts come from PacBio or RNASeq 
# third column are the counts 
# -Genes = keep the column 
MCount_all_norm_long<- gather(data = MCount_all_norm, key = Technology, value = Counts,
                              -Genes)

# need to convert from character to factor for heatmap 
class(MCount_all_norm_long$Genes)
class(MCount_all_norm_long$Technology)
MCount_all_norm_long$Genes <- as.factor(MCount_all_norm_long$Genes)
MCount_all_norm_long$Technology <- as.factor(MCount_all_norm_long$Technology)
MCount_all_norm_long$Counts <- as.numeric(as.character(MCount_all_norm_long$Counts))

library("ggplot2")
library("digest")
MCount_all.heatmap <- ggplot(data = MCount_all_norm_long, 
                                  mapping = aes(x = Technology, y = Genes,
                                                fill = Counts)) + 
  geom_tile() + xlab(label = "Sequencing Platform")

MCount_all.heatmap

```