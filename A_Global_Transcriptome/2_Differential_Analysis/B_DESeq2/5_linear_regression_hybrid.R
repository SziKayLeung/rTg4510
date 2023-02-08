## ---------- Script -----------------
##
## Purpose: perform differential analysis on mouse rTg4510 ONT targeted datasets using linear regression
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))
library(pheatmap)
library(RColorBrewer)

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/transcriptome_stats/read_sq_classification.R")

## ---------- input -----------------

dirnames <- list(
  rTg4510 = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome/",
  whole = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/"
)

input_files <- list(
  phenotype = paste0(dirnames$whole, "3_differential/1_Input/B_hybrid_59/WholeAllMouse_PhenotypeTAPPAS.txt"), 
  expression = paste0(dirnames$whole, "3_differential/1_Input/B_hybrid_59/WholeMouseRNASeq_sqantisubset.expression.txt"),
  classfiles = paste0(dirnames$whole, "2_post_isoseq3/9_sqanti3/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt")
)

input <- list()
input$phenotype <- read.table(input_files$phenotype, sep = "\t", header = T)
input$expression <- read.table(input_files$expression, sep = "\t", header = T)  %>% tibble::column_to_rownames("X")
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")


# input phenoptype characteristics as factor
# input phenotpye sample ID has to match the sample ID in expression matrix 
# samples include "ont_" in expression matrix
input$phenotype$time <- as.factor(input$phenotype$time)
input$phenotype$group <- as.factor(input$phenotype$group)
input$phenotype$group <- relevel(input$phenotype$group,"CONTROL")
input$phenotype$sample <- str_remove(as.character(input$phenotype$sample),"FL.")

# expression matrix to include only ONT samples
# ensure expression column and phenotype rows are in the same order
input$expression <- input$expression %>% select(input$phenotype$sample)
print(ncol(input$expression))


## ---------- Creating DESeq2 object and analysis -----------------

# create DESeq2 data object
# normalises reads within DESeq2
# design is linear model – last expression is tested for differential expression
dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(input$expression)), 
                              colData = input$phenotype, 
                              design = ~group + time + group:time)


# estimate size factors to account for differences in sequencing depth
# if all samples have exact same sequencing depth, size factor should be ~ 1
dds <- estimateSizeFactors(dds)

# pre-filtering the data set
# remove the rows of the DESeqDataSet that have no counts, or only a single count across all samples
# arbitrary threshold = 4 
# minimum 2 samples per group, and expect minimum 2 FL reads per isoform
cat("No of isoforms before filtering:", nrow(dds),"\n")
dds <- dds[rowSums(counts(dds)) > 4, ] 
cat("No of isoforms after filtering:", nrow(dds),"\n")

p <- cbind(reshape2::melt(sizeFactors(dds)), reshape2::melt(colSums(counts(dds)))) %>% 
  `colnames<-`(c("sizefactors", "nreads")) %>% 
  tibble::rownames_to_column("sample") %>% 
  mutate(sample = str_remove(sample,"ont_")) %>%
  ggplot(., aes(x = sizefactors, y = nreads)) + geom_point() +
  geom_label_repel(aes(label = sample), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  theme_bw() + labs(y = "Number of reads", x = "Size Factors")
p  

# normalization to stabilize variance (regularized logarithm)
rld <- rlog(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(rld, intgroup = c("time", "group"), returnData = TRUE)
pcaData$sample <- sapply(pcaData$name, function(x) str_remove(x, "ont_"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = time, shape = group.1)) +
  geom_point(size =3) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_label_repel(aes(label = sample), segment.color = 'grey50') 
p1  

# run differential analysis
colData(dds)$group <- relevel(colData(dds)$group, "CONTROL")
design(dds) <- ~ group + time + group:time
dds <- DESeq(dds, test="Wald")
res <- results(dds)
summary(res)
res_df <- as.data.frame(res)

# annotate results 
res_df <- merge(res_df, input$classfiles[,c("isoform","associated_gene")], by.x = 0, by.y = "isoform") %>% arrange(padj)

# below 0.1 threshold (p-value)
sig0.1 <- res_df %>% filter(padj < 0.1) 
nrow(sig0.1)
