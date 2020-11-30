# Szi Kay Leung
# Aim: To run differential gene expression using DESeq2 on Iso-Seq data alone 
# 02/11/2020: Using output from SQANTI3 (chained data as input) with FL read count for differential gene expression analysis

library("dplyr")
library("reshape2")
library("tidyr")
library("DESeq2")
library("tibble")
library("vsn")
library("pheatmap")
library("gplots")


# read in SQANTI2 classification file of all merged data
class <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN_OLD/SQANTI3/all_samples.chained.rep_classification.filtered_lite_classification.txt",header=T)

# datawrangle for input to DESeq2 
class[["transcript_name_id"]] <- paste0(class$associated_transcript,"_", class$isoform)
countdata <- class[,c("transcript_name_id","FL.O18","FL.K18","FL.S18","FL.L22","FL.Q20","FL.K24","FL.Q21","FL.K17","FL.M21","FL.O23","FL.S23","FL.K23")] %>% column_to_rownames(., var = "transcript_name_id")
countdata <- as.matrix(countdata)
(condition <- factor(c(rep("TG", 6), rep("WT", 6))))

# Make DESeq dataset
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Optional filtering step to remove very low counts (chosen minimum of 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2 pipeline
dds <- DESeq(dds)
res <- DESeq2::results(dds)
res

# shrink 
resLFC <- lfcShrink(dds, coef="condition_WT_vs_TG", type="apeglm")
resLFC

# Save normalised counts separately, used to average technical replicates
#norm_counts <- as.data.frame(counts(dds, normalized=TRUE))
#write.csv(norm_counts, file="norm_counts.txt")

#DESeq2 results
table(res$padj<0.05)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Transcript"
head(resdata)

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot", genecol = "black", fitcol = "cyan", finalcol = "blue", legend = TRUE)

# RLD for viewing
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

# Plot residual p-values
hist(res$pvalue, breaks=50, col="grey")

# We recommend instead using the varianceStabilizingTransformation or shifted log (see vignette).
# plotSparsity(dds)

mycols <- brewer.pal(8, "Accent")[1:length(unique(condition))]

# Heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")

# PCA
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomright", main="Principal Component Analysis", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1: ",as.character(pc1var),"% variance")
  pc2lab <- paste0("PC2: ",as.character(pc2var),"% variance")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-15, 0), ylim=c(-2, 3))


# MA plot
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="blue", pch=20, cex=1))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot-0.05.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

plotMA(res)
plotMA(resLFC)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

library("pheatmap")
ntd <- normTransform(dds)
library("vsn")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df)[1] <- "condition"
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE)

# clustering
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# Volcano Plot
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, xlab="log2(Fold Change)", legendpos="topright", labelsig=FALSE, textcx=1.5, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlab=xlab, cex.axis=1.8, cex.lab=1.5, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p-adj<",sigthresh,sep=""), paste("|log2(FC)|>",lfcthresh,sep=""), "both"), cex=1.5, pch=20, col=c("blue","orange","green"))
}
volcanoplot(resdata, lfcthresh=2, sigthresh=0.05, xlim=c(-30, 30), ylim=c(0,30), legendpos="topright")

mm10_reference_file <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.vM22_gene_annotation_table.txt"
mm10_reference <- read.table(mm10_reference_file,as.is = T, header=T, sep = "\t")
DTU <- merge(DTU, mm10_reference[c("gene_id","GeneSymbol")], by.x = ("associated_gene"), by.y= "gene_id", all.x = T) %>% 
  .[,c("isoform","GeneSymbol","associated_gene","associated_transcript","structural_category","subcategory",
       "u","WT_mean_exp","TG_mean_exp","WT_median_exp","TG_median_exp","mean_exp_diff","Direction",FL_count_colnames)]


write.csv(as.data.frame(res), file="DESeq2_IsoseqFL_TranscriptExpression.csv")
