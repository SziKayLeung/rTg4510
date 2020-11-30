# Szi Kay Leung
# Aim: To run differential gene expression using DESeq2 on Iso-Seq data alone 
# 02/11/2020: Using output from SQANTI3 (chained data as input) with FL read count for differential gene expression analysis

library("dplyr")
library("reshape2")
library("tidyr")
library("DESeq2")
library("tibble")
library("vsn")

# read in SQANTI2 classification file of all merged data
class <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN_OLD/SQANTI3/all_samples.chained.rep_classification.filtered_lite_classification.txt",header=T)

# datawrangle for input to DESeq2 
gene_level <- class %>% 
  # keep only columns with isoform FL counts
  select(associated_gene,FL.K17,FL.K18,FL.K23,FL.K24,FL.L22,FL.M21,FL.O18,FL.O23,FL.Q20,FL.Q21,FL.S18,FL.S23) %>% 
  #select(associated_gene,FL_TPM.K17,FL_TPM.K18,FL_TPM.K23,FL_TPM.K24,FL_TPM.L22,FL_TPM.M21,FL_TPM.O18,FL_TPM.O23,FL_TPM.Q20,FL_TPM.Q21,FL_TPM.S18,FL_TPM.S23) %>% 
  # summarise the number of isoform FL counts per gene per sample: wide to long
  melt() %>% group_by(associated_gene, variable) %>% tally(value) %>%
  # reformat for count matrix input to DESeq2: long to wide
  spread(., variable, n) %>% 
  # ensure order of columns = order of coldata 
  .[,c("associated_gene","FL.O18","FL.K18","FL.S18","FL.L22","FL.Q20","FL.K24","FL.Q21","FL.K17","FL.M21","FL.O23","FL.S23","FL.K23")] %>% 
  #.[,c("associated_gene","FL_TPM.O18","FL_TPM.K18","FL_TPM.S18","FL_TPM.L22","FL_TPM.Q20","FL_TPM.K24","FL_TPM.Q21","FL_TPM.K17","FL_TPM.M21","FL_TPM.O23","FL_TPM.S23","FL_TPM.K23")] %>%
  column_to_rownames(., var = "associated_gene")

# convert FL counts to TPM and round to for integer
gene_level <- as.data.frame(apply(gene_level,2, function(x) round(x/sum(x)*1000000,0)))
coldata <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/Diff_Analysis/condition.txt", row.names = NULL) %>% select(-age)
colnames(coldata)[1] <- ""

### DESeq2 ###
# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(gene_level), colData = coldata, design = ~ genotype)

# set WT as reference level for comparison
dds$genotype <- factor(dds$genotype, levels = c("WT","TG"))

# first attempt with no prefiltering of reads 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Differential gene expression
dds <- DESeq(dds)
res <- results(dds)
res
reference <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.vM22_gene_annotation_table.txt", as.is = T, sep = "\t", header = T)

final_res <- merge(as.data.frame(res),reference[,c("gene_id","GeneSymbol")], by.x = 0, by.y = "gene_id")
final_resSig <- subset(final_res, padj < 0.05)
write.csv(as.data.frame(final_resSig), file="condition_treated_results.csv")

# Log fold change shrinkage for visualization and ranking (using normal)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="genotype_TG_vs_WT",type="apeglm")
resLFC

# plot 
# Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res, ylim=c(-2,2))
# More useful to visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds
plotMA(resLFC, ylim=c(-2,2))


# heatmap
library("pheatmap")
coldata <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/Diff_Analysis/condition.txt", row.names = NULL) 
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- coldata[,c("genotype","age")]
rownames(df) <- colnames(assay(ntd)[select,])
pheatmap(assay(ntd)[select,],show_rownames=FALSE,cluster_rows=FALSE, annotation_col=df)


# data transformation
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))
meanSdPlot(assay(ntd))
