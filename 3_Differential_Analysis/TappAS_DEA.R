# Szi Kay Leung: sl693@exeter.ac.uk
# Differential Gene Expression Analysis (Whole Transcriptome, Targeted Transcriptome)

#### Analysis ################
# 1. Whole Transcriptome: Iso-Seq Scaffold + RNASeq Expression abundance 
# 2. Whole Transcriptome: Iso-Seq Scaffold + IsoSeq FL Read count 
# 3. Targeted Transcriptome (20 Genes): Iso-Seq Scaffold + RNASeq Expression abundance
# 4. Targeted Transcriptome (20 Genes): Iso-Seq Scaffold + IsoSeq FL Read count

# Prerequisite: run tappAS and transfer gene_matrix.tsv, transcript_matrix.tsv, time_factors.txt from project directory (knight) to ISCA 
####

# library 
suppressMessages(library("maSigPro"))
suppressMessages(library("mclust"))
suppressMessages(library(xlsx))


# input and output directory paths
wholeindir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/TAPPAS/TAPPAS_output"
targetedindir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAPPAS/TAPPAS_output"
indir = list(paste0(wholeindir,"/IsoSeq"),paste0(wholeindir,"/RNASeq"),
             paste0(targetedindir,"/IsoSeq"),paste0(targetedindir,"/RNASeq"))
names(indir) = c("WholeIso","WholeRNA","TargetedIso","TargetedRNA")
outdir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures"

# run time course DEA analysis
time_analysis <- function(dtname, degree, siglevel, r2cutoff, knum, usemclust, groups, indir) {
  # read expression factors
  cat("\nReading factors file data...")
  factors=read.table(file.path(indir, "time_factors.txt"), row.names=1, sep="\t", header=TRUE)
  groups = ncol(factors) - 2
  times = length(unique(factors[,1]))
  
  # read expression matrix
  cat("\nReading normalized ", dtname, " matrix file data...")
  expMatrix = read.table(file.path(indir, paste0(dtname, "_matrix.tsv")), row.names=1, sep="\t", header=TRUE)
  cat("\nRead ", nrow(expMatrix), " normalized ", dtname, " expression data rows")
  
  # create a regression matrix for regression model 
  design <- make.design.matrix(factors, degree)
  minobs = (degree + 1) * groups
  #cat("\ngroups: ", groups, ", degree: ", degree, ", minobs: ", minobs)
  #cat("\nsiglevel: ", siglevel, ", r2cutoff: ", r2cutoff, ", knum: ", knum, ", useMclust: ", usemclust)
  #print(design$edesign)
  #cat("\nRunning p.vector()...\n")
  
  ### Finding significant genes
  # regression fit for each gene for each condition 
  # p value for each regression fit associated with F statistic after correction with BH FDR
  # counts = TRUE uses a GLM and applies NB distribution
  # A gene with different profile (regression fit) between 2 conditions will have statistically significant coefficient
  fit <- p.vector(expMatrix, design, Q=siglevel, MT.adjust="BH", min.obs=minobs, counts=TRUE)
  cat("\nFound ", fit$i, " out of ", fit$g, " items to be significant\n")
  #cat("\ndesign:\n")
  print(head(fit$SELEC))
  #fit$p.vector
  
  ### Finding significant differences
  # stepwise regression - iterative regression approach that selects from a pool of potential variables the ‘best’ ones (according to a specified criterion) to fit the available data.
  # Apply a threshold of Rsqaured coefficient cutoff = goodness of fit for meaningful differences
  
  cat("\nRunning T.fit() using step.method='backward' and alfa = ", siglevel, "...\n")
  tstep <- T.fit(fit, step.method="backward", alfa=siglevel)
  tstep$sol
  
  cat("\nRunning get.siggenes() using rsq=", r2cutoff, " and vars='groups'...\n")
  sigs <- get.siggenes(tstep, rsq=r2cutoff, vars="groups")

  for(i in 1:length(names(sigs$sig.genes))) {
    gname <- names(sigs$sig.genes)[i]
    cat("\n", gname, " significant ", dtname, ": ")
    print(sigs$sig.genes[[gname]]$g)
    df <- data.frame(sigs$sig.genes[[gname]]$sig.pvalues[, 1:2])
    df <- df[order(df$p.value), ]
    #df <<- df
    print(nrow(df))
    print(head(df))
    cat("\nSaving ", dtname, " DEA results to file...")
  }
  return(tstep$sol)
}

degree = 3
knum = 9
mvalue = 1
siglevel = 0.05
r2cutoff = 0.5
usemclust = TRUE
outdir = indir

gene_sigs = list()
count = 1
for(dir in c(indir$TargetedIso,indir$TargetedRNA,indir$WholeIso,indir$WholeRNA)){
  gene_sigs[[count]] = time_analysis(dtname="gene",degree=degree, siglevel=siglevel, r2cutoff=r2cutoff, knum=knum, usemclust=usemclust, groups=groups, indir=dir)
  count = count + 1
}
names(gene_sigs) = c("TargetedIso","TargetedRNA","WholeIso","WholeRNA")


transcript_sigs = time_analysis(dtname="transcript",factors=myfactors, degree=degree, siglevel=siglevel, r2cutoff=r2cutoff, knum=knum, usemclust=usemclust, groups=groups, outdir=outdir)

write.table(gene_sigs,paste0(outdir,"/GeneDEA_Results.tsv"))
write.table(transcript_sigs,paste0(outdir,"/TranscriptsDEA_Results.tsv"))

# save results as workbook
wb <- createWorkbook()
datas <- c(as.data.frame(gene_sigs$TargetedIso),as.data.frame(gene_sigs$TargetedRNA),as.data.frame(gene_sigs$WholeIso),as.data.frame(gene_sigs$WholeRNA))
sheetnames <- c("TargetedIso_Genexp","TargetedRNA_Genexp","WholeIso_Genexp","WholeRNA_Genexp")
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"))

write.xlsx(gene_sigs$TargetedIso, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedIso_Genexp")
write.xlsx(gene_sigs$TargetedRNA, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedRNA_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeIso, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIso_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeRNA, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeRNA_Genexp", append=TRUE)

# find the isoforms that are differentially expressed
targeted.class.files = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt", sep = "\t",as.is = T, header = T)

IsoformDEA = row.names(transcript_sigs[transcript_sigs$`R-squared` > r2cutoff,])
targeted.class.files[targeted.class.files$isoform %in% IsoformDEA,]


