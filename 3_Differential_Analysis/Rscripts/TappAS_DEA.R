# Szi Kay Leung: sl693@exeter.ac.uk
# Differential Gene Expression Analysis (Whole Transcriptome, Targeted Transcriptome)

#### Analysis ################
# 1. Whole Transcriptome: Iso-Seq Scaffold + RNASeq Expression abundance 
# 2. Whole Transcriptome: Iso-Seq Scaffold + IsoSeq FL Read count 
# 3. Targeted Transcriptome (20 Genes): Iso-Seq Scaffold + RNASeq Expression abundance
# 4. Targeted Transcriptome (20 Genes): Iso-Seq Scaffold + IsoSeq FL Read count
# 5. Targeted Transcriptome (20 Genes): Iso-Seq Scaffold + IsoSeq FL Read count Collapsed from partial transcripts

# Prerequisite: run tappAS and transfer gene_matrix.tsv, transcript_matrix.tsv, time_factors.txt from project directory (knight) to ISCA 
####

# library 
suppressMessages(library("maSigPro"))
suppressMessages(library("mclust"))
suppressMessages(library("xlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))


# input and output directory paths
wholeindir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT"
targetedindir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_OUTPUT"
indir = list(paste0(wholeindir,"/IsoSeq_Expression"),paste0(wholeindir,"/RNASeq_Expression"),
             paste0(targetedindir,"/IsoSeq_Expression"),paste0(targetedindir,"/RNASeq_Expression"),
             paste0(targetedindir,"/IsoSeq_Expression_Collapsed"),paste0(targetedindir,"/RNASeq_Expression_Collapsed"))
names(indir) = c("WholeIso","WholeRNA","TargetedIso","TargetedRNA","TargetedIsoCol","TargetedRNACol")
outdir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures"

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
targeted.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt"
targeted.class.files <- SQANTI_class_preparation(targeted.class.names.files,"standard")

targeted.collapsed.class.names.files = "//gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.classification.txt"
targeted.collapsed.class.files <- SQANTI_class_preparation(targeted.collapsed.class.names.files,"standard")


whole.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantitamafiltered.classification.txt"
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

# target genes 
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")

# run time course DEA analysis
time_analysis <- function(dtname, degree, siglevel, r2cutoff, knum, usemclust, groups, indir) {
  # read expression factors
  cat("\nReading factors file data...")
  factors=read.table(file.path(indir, "/Data/time_factors.txt"), row.names=1, sep="\t", header=TRUE)
  groups = ncol(factors) - 2
  times = length(unique(factors[,1]))
  
  # read expression matrix
  cat("\nReading normalized ", dtname, " matrix file data...")
  expMatrix = read.table(file.path(indir, paste0("/Data/",dtname, "_matrix.tsv")), row.names=1, sep="\t", header=TRUE)
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

### Apply function
# gene expression 
gene_sigs = list()
count = 1
for(dir in c(indir$TargetedIso,indir$TargetedRNA,indir$WholeIso,indir$WholeRNA,indir$TargetedIsoCol, indir$TargetedRNACol)){
  
  # different degree of freedom depending on expression input (WholeIso = 1 as only 2 age groups, RNA & TargetedIso = 3 as 4 age groups)
  if(count %in% c("3")){
    cat("Processing IsoSeq:", dir, "\n")
    gene_sigs[[count]] = time_analysis(dtname="gene",degree=1, siglevel=0.05, r2cutoff=0.5, knum=9, usemclust=TRUE, groups=groups, indir=dir)
  }else{
    cat("Processing RNASeq:", dir, "\n")
    gene_sigs[[count]] = time_analysis(dtname="gene",degree=3, siglevel=0.05, r2cutoff=0.5, knum=9, usemclust=TRUE, groups=groups, indir=dir)
  }
  count = count + 1
}
names(gene_sigs) = c("TargetedIso","TargetedRNA","WholeIso","WholeRNA","TargetedIsoCol","TargetedRNACol")
gene_sigs$WholeIso_Target = subset(gene_sigs$WholeIso, rownames(gene_sigs$WholeIso) %in% TargetGene)
gene_sigs$WholeRNA_Target = subset(gene_sigs$WholeRNA, rownames(gene_sigs$WholeRNA) %in% TargetGene)

### transcript expression
trans_sigs = list()
count = 1
for(dir in c(indir$TargetedIso,indir$TargetedRNA,indir$WholeIso,indir$WholeRNA,indir$TargetedIsoCol, indir$TargetedRNACol)){
  # WholeIso = Whole Transcriptome with Iso-Seq expression across 2 time points ==> 1 degree of freedom
  # WholeRNA = Whole Transcriptome with RNA-Seq expression across 4 time points ==> 3 degrees of freedom
  # TargetedRNA, TargetedIso, TargetedIsoCol, TargetedRNACol = Targeted Transcriptome with RNA-Seq expression across 4 time points 
  if(count %in% c("3")){degreespec=1}else{degreespec=3}
    
  trans_sigs[[count]] = time_analysis(dtname="transcript",degree=degreespec, siglevel=0.05, r2cutoff=0.5, knum=9, usemclust=TRUE, groups=groups, indir=dir)
  count = count + 1
}
names(trans_sigs) = c("TargetedIso","TargetedRNA","WholeIso","WholeRNA","TargetedIsoCol","TargetedRNACol")

# annotate transcripts from sqanti classification file 
trans_sigsanno = list()
count = 1
for(dataset in trans_sigs){
  
  if(count %in% c("1","2")){
    # Targeted TAMA filtering 
    cat("Annotation with TAMA Filtering SQANTI file: Targeted","\n")
    dat = targeted.class.files
  }else if(count %in% c("3","4")){
    cat("Annotation with TAMA Filtering SQANTI file: Whole","\n")
    dat = whole.class.files
  }else if(count %in% c("5","6")){
    cat("Annotation with Subset Collapsed file after SQANTI, TAMA filtering: Targeted","\n")
    dat = targeted.collapsed.class.files
  }else{print("NULL")}
  
  trans_sigsanno[[count]] = merge(dat[,c("isoform","associated_gene","associated_transcript","structural_category")],dataset, by.x = "isoform", by.y = 0) %>% filter(`p-value` != "NA", `R-squared` > 0.5)
  count = count + 1
}
names(trans_sigsanno) = c("TargetedIso","TargetedRNA","WholeIso","WholeRNA","TargetedIsoCol", "TargetedRNACol")
trans_sigsanno$WholeIso_Target = subset(trans_sigsanno$WholeIso, associated_gene %in% TargetGene)
trans_sigsanno$WholeRNA_Target = subset(trans_sigsanno$WholeRNA, associated_gene %in% TargetGene)

              
### Save results as workbook
write.xlsx(gene_sigs$TargetedIso, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedIso_Genexp")
write.xlsx(gene_sigs$TargetedRNA, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedRNA_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeIso, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIso_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeRNA, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeRNA_Genexp", append=TRUE)
write.xlsx(gene_sigs$TargetedIsoCol, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedIsoCol_Geneexp", append=TRUE)
write.xlsx(gene_sigs$TargetedRNACol, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedRNACol_Geneexp", append=TRUE)
write.xlsx(gene_sigs$WholeIso_Target, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIsoTarget_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeRNA_Target, paste0(outdir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeRNATarget_Genexp", append=TRUE)



# transcript expression 
write.xlsx(trans_sigsanno$TargetedIso, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="TargetedIso_Transexp")
write.xlsx(trans_sigsanno$TargetedRNA, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="TargetedRNA_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$WholeIso, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeIso_Transexp",append=TRUE)
write.xlsx(trans_sigsanno$WholeRNA, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeRNA_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$TargetedIsoCol, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="TargetedIsoCol_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$WholeIso_Target, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeIsoTargeted_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$WholeRNA_Target, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeRNATargeted_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$WholeRNA_Target, paste0(outdir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeRNATargeted_Transexp", append=TRUE)

