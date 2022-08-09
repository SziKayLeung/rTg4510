# Szi Kay Leung: sl693@exeter.ac.uk
# Differential Gene, Transcript Expression and Isoform Usage Analysis (Whole Transcriptome)


#### Analysis ################
## Differential Gene and Transcript Expression Analysis
# 1. Targeted Transcriptome: Iso-Seq Scaffold + IsoSeq FL Read count 

## Differential Isoform Usage Analysis 
# 2. Targeted Transcriptome: Iso-Seq Scaffold + IsoSeq FL Read count 

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

output_tappas_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ/TAPPAS_OUTPUT"
input.classfile = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ/SQANTI3/AllMouseTargeted_ISMrem.classification.txt"

# input and output directory paths
indir = list(paste0(output_tappas_dir,"/IsoSeq_Expression"))
names(indir) = c("TargetedIso")
funcdir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Rscripts"

# sqanti files
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
class.files <- SQANTI_class_preparation(input.classfile,"standard")

# source function scripts
source(paste0(funcdir,"/TappAS_DEA.R"))
source(paste0(funcdir,"/TappAS_DIU.R"))

##### Differential Analysis
# Targeted Transcriptome = Targeted Transcriptome with Iso-Seq expression across 4 time points ==> 3 degrees of freedom

#### Differential Gene Expression Analysis ##############################c
cat("Processing Differential Gene Expression Analysis \n")
g_sigs = list()
g_sigs[["TargetedIso"]] = DEA_time_analysis(dtname="gene",degree=3,siglevel=0.05,r2cutoff=0.5,knum=9,usemclust=TRUE,groups=groups,indir=indir$TargetedIso)

### Gene expression Output
write.xlsx(g_sigs$TargetedIso, paste0(output_tappas_dir,"/DifferentialGeneExpression_Analysis.xlsx"), sheetName="TargetedAllIso_Genexp")

#### Differential Transcript Expression Analysis ##############################
cat("Processing Differential Transcript Expression Analysis \n")
t_sigs = list()
t_sigs[["TargetedIso"]]=DEA_time_analysis(dtname="transcript",degree=3,siglevel=0.05,r2cutoff=0.5,knum=9, usemclust=TRUE,groups=groups,indir=indir$TargetedIso)

# annotate transcripts from sqanti classification file 
class.files = class.files[,c("isoform","associated_gene","associated_transcript","structural_category")]
t_sigsanno = lapply(t_sigs, function(x) merge(class.files,x, by.x = "isoform", by.y = 0) %>% filter(`p-value` != "NA", `R-squared` > 0.5))

# Transcript expression Output
write.xlsx(t_sigsanno$TargetedIso, paste0(output_tappas_dir,"/DifferentialTransExpression_Analysis.xlsx"), sheetName="TargetedIso_Transexp")

cat("All Done!\n")
