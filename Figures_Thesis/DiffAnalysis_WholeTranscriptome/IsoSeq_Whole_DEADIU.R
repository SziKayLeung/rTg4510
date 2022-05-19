# Szi Kay Leung: sl693@exeter.ac.uk
# Differential Gene, Transcript Expression and Isoform Usage Analysis (Whole Transcriptome)


#### Analysis ################
## Differential Gene and Transcript Expression Analysis
# 1. Whole Transcriptome: Iso-Seq Scaffold + RNASeq Expression abundance 
# 2. Whole Transcriptome: Iso-Seq Scaffold + IsoSeq FL Read count 

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
suppressMessages(library("tibble"))

args = commandArgs(trailingOnly=TRUE) 
output_tappas_dir <- args[1]
input.classfile <- args[2]
output_dir <- args[3]

output_tappas_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/TAPPAS_OUTPUT"
input.classfile = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/SQANTI3/WholeIsoSeq_ISMrem.classification.txt"
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/TAPPAS_OUTPUT"


# input and output directory paths
indir = list(paste0(output_tappas_dir,"/IsoSeq_Expression"),
             paste0(output_tappas_dir,"/RNASeq_Expression"))
names(indir) = c("WholeIso","WholeRNA")
funcdir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Rscripts"

# sqanti files
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
whole.class.files <- SQANTI_class_preparation(input.classfile,"standard")

# source function scripts
source(paste0(funcdir,"/TappAS_DEA.R"))
source(paste0(funcdir,"/TappAS_DIU.R"))

##### Differential Analysis
# WholeIso = Whole Transcriptome with Iso-Seq expression across 2 time points ==> 1 degree of freedom
# WholeRNA = Whole Transcriptome with RNA-Seq expression across 4 time points ==> 3 degrees of freedom

#### Differential Gene Expression Analysis ##############################c
cat("Processing Differential Gene Expression Analysis \n")
g_sigs = list()
g_sigs[["WholeIso"]] = DEA_time_analysis(dtname="gene",degree=1,siglevel=0.05,r2cutoff=0.5,knum=9,usemclust=TRUE,groups=groups,indir=indir$WholeIso)
g_sigs[["WholeRNA"]] = DEA_time_analysis(dtname="gene",degree=3,siglevel=0.05,r2cutoff=0.5,knum=9,usemclust=TRUE,groups=groups,indir=indir$WholeRNA)
g_sigs[["WholeIso_All"]] = DEA_time_analysis(dtname="gene",degree=1,siglevel=1,r2cutoff=0.02,knum=9,usemclust=TRUE,groups=groups,indir=indir$WholeIso)

# target genes 
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu",
                "Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")
g_sigs$WholeIso_Target = subset(g_sigs$WholeIso, rownames(g_sigs$WholeIso) %in% TargetGene)
g_sigs$WholeRNA_Target = subset(g_sigs$WholeRNA, rownames(g_sigs$WholeRNA) %in% TargetGene)

### Gene expression Output
write.xlsx(g_sigs$WholeIso, paste0(output_dir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIso_Genexp")
write.xlsx(g_sigs$WholeRNA, paste0(output_dir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeRNA_Genexp", append=TRUE)
write.xlsx(g_sigs$WholeIso_Target, paste0(output_dir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIsoTarget_Genexp", append=TRUE)
write.xlsx(g_sigs$WholeRNA_Target, paste0(output_dir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeRNATarget_Genexp", append=TRUE)
write.xlsx(g_sigs$WholeIso_All, paste0(output_dir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="WholeIsoAll_Genexp", append=TRUE)


#### Differential Transcript Expression Analysis ##############################
cat("Processing Differential Transcript Expression Analysis \n")
t_sigs = list()
t_sigs[["WholeIso"]]=DEA_time_analysis(dtname="transcript",degree=1,siglevel=0.05,r2cutoff=0.5,knum=9, usemclust=TRUE,groups=groups,indir=indir$WholeIso)
t_sigs[["WholeRNA"]]=DEA_time_analysis(dtname="transcript",degree=3,siglevel=0.05,r2cutoff=0.5,knum=9, usemclust=TRUE,groups=groups,indir=indir$WholeRNA)

# annotate transcripts from sqanti classification file 
whole.class.files = whole.class.files[,c("isoform","associated_gene","associated_transcript","structural_category")]
t_sigsanno = lapply(t_sigs, function(x) merge(whole.class.files,x, by.x = "isoform", by.y = 0) %>% filter(`p-value` != "NA", `R-squared` > 0.5))

t_sigsanno$WholeIso_Target = subset(t_sigsanno$WholeIso, associated_gene %in% TargetGene)
t_sigsanno$WholeRNA_Target = subset(t_sigsanno$WholeRNA, associated_gene %in% TargetGene)

# Transcript expression Output
write.xlsx(t_sigsanno$WholeIso, paste0(output_dir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeIso_Transexp")
write.xlsx(t_sigsanno$WholeRNA, paste0(output_dir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeRNA_Transexp", append=TRUE)
write.xlsx(t_sigsanno$WholeIso_Target, paste0(output_dir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeIsoTargeted_Transexp", append=TRUE)
write.xlsx(t_sigsanno$WholeRNA_Target, paste0(output_dir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="WholeRNATargeted_Transexp", append=TRUE)


#### Differential Transcript Usage ##############################
cat("Processing Differential Transcript Usage Analysis \n")
DIU_Isoseq = DIU_time_analysis(indir$WholeIso,degreevalue=1)
DIU_RNAseq = DIU_time_analysis(indir$WholeRNA,degreevalue=3)

gene_matrix=read.table(file.path(indir$WholeRNA,"/Data/gene_matrix.tsv")) 
myfactors=read.table(file.path(indir$WholeRNA, "Data/time_factors.txt"), row.names=1, sep="\t", header=TRUE)
groups = ncol(myfactors) - 2
times = length(unique(myfactors[,1]))
if(groups==1){timepoints = times}

IsoRNA_transMatrix_fc = DIU_prefilter(indir$WholeRNA, mff=2.5, filteringType = "FOLD")
DIU_rnaseq_fc = run_DIU_analysis(degree = "3", triple = FALSE, IsoRNA_transMatrix_fc, myfactors)
DIU_rnaseq_fc = merge_exp(DIU_rnaseq_fc,gene_matrix)

### DIU output
write.xlsx(DIU_Isoseq$DIU_prop, paste0(output_dir,"/DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_isoseq_prop")
write.xlsx(DIU_Isoseq$DIU_fc, paste0(output_dir,"/DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_isoseq_fc", append=TRUE)
write.xlsx(DIU_RNAseq$DIU_prop, paste0(output_dir,"/DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_rnaseq_prop", append=TRUE)
write.xlsx(DIU_rnaseq_fc, paste0(output_dir,"/DifferentialTranscriptUsage_Analysis.xlsx"), sheetName="DIU_rnaseq_fc", append=TRUE)

cat("All Done!\n")