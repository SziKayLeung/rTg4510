# Szi Kay Leung: sl693@exeter.ac.uk
# Differential Gene, Transcript Expression and Isoform Usage Analysis (Whole Transcriptome, Targeted Transcriptome)


#### Analysis ################
## Differential Gene and Transcript Expression Analysis
# 1. Whole Transcriptome: Iso-Seq Scaffold + RNASeq Expression abundance 
# 2. Whole Transcriptome: Iso-Seq Scaffold + IsoSeq FL Read count 
# 3. Targeted Transcriptome (20 Genes): Iso-Seq Scaffold + RNASeq Expression abundance
# 4. Targeted Transcriptome (20 Genes): Iso-Seq Scaffold + IsoSeq FL Read count

## Differential Isoform Usage Analysis 
# 5. Whole Transcriptome: Iso-Seq Scaffold + RNASeq Expression abundance 
# 6. Whole Transcriptome: Iso-Seq Scaffold + IsoSeq FL Read count 

# Prerequisite: run tappAS and ensure gene_matrix.tsv, transcript_matrix.tsv, time_factors.txt from project directory (knight) transferred to ISCA 
####

# library 
suppressMessages(library("maSigPro"))
suppressMessages(library("mclust"))
suppressMessages(library("xlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("MASS"))


# input and output directory paths
wholeindir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_noRNASEQ"
targetedindir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ"
indir = list(paste0(wholeindir,"/TAPPAS_OUTPUT/IsoSeq_Expression"),paste0(wholeindir,"/TAPPAS_OUTPUT/RNASeq_Expression"),
             paste0(targetedindir,"/TAPPAS_OUTPUT/IsoSeq_Expression"),paste0(targetedindir,"/TAPPAS_OUTPUT/RNASeq_Expression"),
             paste0(targetedindir,"/TAPPAS_OUTPUT/IsoSeq_Expression_All"))
names(indir) = c("WholeIso","WholeRNA","TargetedIso","TargetedRNA","TargetedIsoAll")
outdir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures"
funcdir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Rscripts"

# sqanti files
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
targeted.class.names.files = paste0(targetedindir, "/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.classification.txt")
targeted.class.files <- SQANTI_class_preparation(targeted.class.names.files,"standard")

whole.class.names.files = paste0(wholeindir, "/COLLAPSE_FILTER/WholeIsoSeq_sqantisubset.classification.txt")
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

# target genes 
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")

# source function scripts
source(paste0(funcdir,"/TappAS_DEA.R"))
source(paste0(funcdir,"/TappAS_DIU.R"))

#### Differential Gene Expression Analysis ##############################
gene_sigs = list()
count = 1
for(dir in c(indir$TargetedIso,indir$TargetedRNA,indir$WholeIso,indir$WholeRNA,indir$TargetedIsoAll)){
  
  # different degree of freedom depending on expression input (WholeIso = 1 as only 2 age groups, RNA & TargetedIso = 3 as 4 age groups)
  if(count %in% c("3")){
    cat("Processing IsoSeq:", dir, "\n")
    gene_sigs[[count]] = DEA_time_analysis(dtname="gene",degree=1, siglevel=0.05, r2cutoff=0.5, knum=9, usemclust=TRUE, groups=groups, indir=dir)
  }else{
    cat("Processing RNASeq:", dir, "\n")
    gene_sigs[[count]] = DEA_time_analysis(dtname="gene",degree=3, siglevel=0.05, r2cutoff=0.5, knum=9, usemclust=TRUE, groups=groups, indir=dir)
  }
  count = count + 1
}
names(gene_sigs) = c("TargetedIso","TargetedRNA","WholeIso","WholeRNA","TargetedIsoAll")
gene_sigs$WholeIso_Target = subset(gene_sigs$WholeIso, rownames(gene_sigs$WholeIso) %in% TargetGene)
gene_sigs$WholeRNA_Target = subset(gene_sigs$WholeRNA, rownames(gene_sigs$WholeRNA) %in% TargetGene)

### Gene expression Output
write.xlsx(gene_sigs$TargetedIso, paste0(outdir,"/Final_DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedIso_Genexp")
write.xlsx(gene_sigs$TargetedRNA, paste0(outdir,"/Final_DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedRNA_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeIso, paste0(outdir,"/Final_DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIso_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeRNA, paste0(outdir,"/Final_DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeRNA_Genexp", append=TRUE)
write.xlsx(gene_sigs$TargetedIsoAll, paste0(outdir,"/Final_DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedIsoAll_Geneexp", append=TRUE)
write.xlsx(gene_sigs$WholeIso_Target, paste0(outdir,"/Final_DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIsoTarget_Genexp", append=TRUE)
write.xlsx(gene_sigs$WholeRNA_Target, paste0(outdir,"/Final_DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeRNATarget_Genexp", append=TRUE)


#### Differential Transcript Expression Analysis ##############################
trans_sigs = list()
count = 1
for(dir in c(indir$TargetedIso,indir$TargetedRNA,indir$WholeIso,indir$WholeRNA,indir$TargetedIsoAll)){
  # WholeIso = Whole Transcriptome with Iso-Seq expression across 2 time points ==> 1 degree of freedom
  # WholeRNA = Whole Transcriptome with RNA-Seq expression across 4 time points ==> 3 degrees of freedom
  # TargetedRNA, TargetedIso, TargetedIsoCol, TargetedRNACol = Targeted Transcriptome with RNA-Seq expression across 4 time points 
  if(count %in% c("3")){degreespec=1}else{degreespec=3}
  
  trans_sigs[[count]] = DEA_time_analysis(dtname="transcript",degree=degreespec, siglevel=0.05, r2cutoff=0.5, knum=9, usemclust=TRUE, groups=groups, indir=dir)
  count = count + 1
}
names(trans_sigs) = c("TargetedIso","TargetedRNA","WholeIso","WholeRNA","TargetedIsoAll")

# annotate transcripts from sqanti classification file 
trans_sigsanno = list()
count = 1
for(dataset in trans_sigs){
  
  if(count %in% c("1","2","6")){
    # Targeted TAMA filtering 
    cat("Annotation with SQANTI file: Targeted","\n")
    dat = targeted.class.files
  }else if(count %in% c("3","4")){
    cat("Annotation with SQANTI file: Whole","\n")
    dat = whole.class.files
  }else{print("NULL")}
  
  trans_sigsanno[[count]] = merge(dat[,c("isoform","associated_gene","associated_transcript","structural_category")],dataset, by.x = "isoform", by.y = 0) %>% 
    filter(`p-value` != "NA", `R-squared` > 0.5)
  count = count + 1
}
names(trans_sigsanno) = c("TargetedIso","TargetedRNA","WholeIso","WholeRNA","TargetedIsoAll")
trans_sigsanno$WholeIso_Target = subset(trans_sigsanno$WholeIso, associated_gene %in% TargetGene)
trans_sigsanno$WholeRNA_Target = subset(trans_sigsanno$WholeRNA, associated_gene %in% TargetGene)

# Transcript expression Output
write.xlsx(trans_sigsanno$TargetedIso, paste0(outdir,"/Final_DifferentialTransExpression_Analysis.xlsx"), sheetName="TargetedIso_Transexp")
write.xlsx(trans_sigsanno$TargetedRNA, paste0(outdir,"/Final_DifferentialTransExpression_Analysis.xlsx"), sheetName="TargetedRNA_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$WholeIso, paste0(outdir,"/Final_DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeIso_Transexp",append=TRUE)
write.xlsx(trans_sigsanno$WholeRNA, paste0(outdir,"/Final_DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeRNA_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$TargetedIsoAll, paste0(outdir,"/Final_DifferentialTransExpression_Analysis.xlsx"), sheetName="TargetedIsoAll_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$WholeIso_Target, paste0(outdir,"/Final_DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeIsoTargeted_Transexp", append=TRUE)
write.xlsx(trans_sigsanno$WholeRNA_Target, paste0(outdir,"/Final_DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeRNATargeted_Transexp", append=TRUE)


#### Differential Isoform Usage ##############################
DIU_Isoseq = DIU_time_analysis(indir$WholeIso,degreevalue=1)
DIU_RNAseq = DIU_time_analysis(indir$WholeRNA,degreevalue=3)

gene_matrix=read.table(file.path(indir$WholeRNA,"/Data/gene_matrix.tsv")) 
myfactors=read.table(file.path(indir$WholeRNA, "Data/time_factors.txt"), row.names=1, sep="\t", header=TRUE)
groups = ncol(myfactors) - 2
times = length(unique(myfactors[,1]))
if(groups==1){timepoints = times}
# cat("\nGroups: ", groups, ", times: ", times)


IsoRNA_transMatrix_fc = DIU_prefilter(indir$WholeRNA, mff=2.5, filteringType = "FOLD")
DIU_rnaseq_fc = run_DIU_analysis(degree = "3", triple = FALSE, IsoRNA_transMatrix_fc, myfactors)
DIU_rnaseq_fc = merge_exp(DIU_rnaseq_fc,gene_matrix)

### DIU output
write.xlsx(DIU_Isoseq$DIU_prop, paste0(outdir,"/Final_DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_isoseq_prop")
write.xlsx(DIU_Isoseq$DIU_fc, paste0(outdir,"/Final_DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_isoseq_fc", append=TRUE)
write.xlsx(DIU_RNAseq$DIU_prop, paste0(outdir,"/Final_DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_rnaseq_prop", append=TRUE)
write.xlsx(DIU_rnaseq_fc, paste0(outdir,"/Final_DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_rnaseq_fc", append=TRUE)

