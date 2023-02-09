## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
## 
##   
##
##

## ---------- Directory -----------------

WKD_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/D_Methylation/"
METADIR_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/D_methylation/"
METH_ROOT <- "/lustre/projects/Research_Project-191406/EmmaW/RRBSAnnotatedResults/rTg4510/DMPs/"

## ---------- Target Genes -----------------

TGENE_DIR <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/B_isoseq_targeted"
TargetGene <- read.table(paste0(TGENE_DIR, "/TargetGenes.tsv"))[["V1"]]

## ---------- Methyation Input files -----------------

# Isabel's RRBS data (Output of load: RRBS_completebetas)
#load(file = paste0(WKD_ROOT, "1_RRBS/rTg4510_RRBSbetasComplete.RData"))
#RRBS_completebetas = RRBS_completebetas %>% tibble::rownames_to_column(var = "position") %>% 
#  mutate(base = word(position,c(2),sep = fixed(":")), chrom = word(position,c(1),sep = fixed(":"))) 
#save(RRBS_completebetas, file = paste0(WKD_ROOT, "1_RRBS/rTg4510_RRBSbetasComplete_positionmod.RData"))
load(file = paste0(WKD_ROOT, "1_RRBS/rTg4510_RRBSbetasComplete_positionmod.RData"))

# Phenotype file
RRBS_Phenotype <- read.csv(paste0(METADIR_ROOT, "Tg4510_phenotype_RRBS.csv"))
#RRBS_Phenotype <- read.csv(paste0(METADIR_ROOT, "Tg4510_coldata_RRBS.csv"), 
#                           row.names=1, stringsAsFactors=FALSE)
RRBS_Phenotype$Age_months <- as.factor(RRBS_Phenotype$Age_months)
RRBS_Phenotype$Genotype <- as.factor(RRBS_Phenotype$Genotype)
RRBS_Phenotype$Genotype <- relevel(RRBS_Phenotype$Genotype, "WT")

# Differentially Methylated Positions for rTg4510 
# Read in files downloaded from RRBS_Paper dropbox folder (Isabel's analysis)
DMPs = list.files(path = paste0(METH_ROOT),pattern = "DMPsInteractionModel", full.names = T)
Whole_DMP = lapply(DMPs, function(x) read.table(x, sep = ",", as.is = T, header = T)) 
names(Whole_DMP) = lapply(list.files(path = paste0(METH_ROOT),pattern = "DMPsInteractionModel"), 
                          function(x) word(x,c(4), sep = fixed("_")))

# Pathology 
Whole_DMP$pathology <- read.csv(paste0(METH_ROOT,"/DMPsPathology_rTg4510_sig_pathology_1500bptssAnno.csv"))

# Differentially Methylated Regions
DMRs <- read.csv(paste0(WKD_ROOT, "3_DMRs/DMRsGenotype_annotated.csv"))
