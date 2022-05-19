# Szi Kay Leung: sl693@exeter.ac.uk
# Global files for analysis 

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")

# root dir
diffroot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/"
rawroot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome"
fig_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/"
tabroot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/TAPPAS_OUTPUT"

# Final classification file
class.names.files = paste0(diffroot_dir, "/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt")
class.files <- SQANTI_class_preparation(class.names.files,"standard")

# Isabel's supplementary table of differentially expressed genes in rTg4510 
Isabel_gene_Tg4510AgeGenotypeDEG = read.csv(paste0(fig_dir, "/Isabel_Supp4_Tg4510AgeGenotypeDEG.csv"), header = T, as.is = T)
Isabel_gene_Tg4510GenotypeDEG = read.csv(paste0(fig_dir, "/Isabel_Supp2_Tg4510GenotypeDEG.csv"), header = T, as.is = T)


#### TAPPAS (Differential Analysis) ###############################################
tappasiso_input_dir <- paste0(diffroot_dir, "TAPPAS_OUTPUT/IsoSeq_Expression")
tappasrna_input_dir <- paste0(diffroot_dir, "TAPPAS_OUTPUT/RNASeq_Expression")
tappasiso_phenotype <- read.table(paste0(rawroot_dir, "/WholeIsoSeq_PhenotypeTAPPAS.txt"), header = T) %>% mutate(variable = paste0(sample))
tappasrna_phenotype <- read.table(paste0(rawroot_dir,"/WholeAllMouse_PhenotypeTAPPAS.txt"), header = T) %>% mutate(col_names = paste0(group,".",sample),age = paste0(time, "_mos"), pheno = ifelse(group == "CONTROL", "WT", "TG"), variable = paste0("FL.", sample,"_",pheno,"_",age))

#### DGE/DTE - Differential Gene and Transcript Expression #########
# Output from Tappas_DEA.R
tappassiggene <- excel_sheets(paste0(tabroot_dir,"/DifferentialGeneExpression_Analysis.xlsx")) %>% set_names() %>% purrr::map(read_excel, path = paste0(tabroot_dir,"/DifferentialGeneExpression_Analysis.xlsx"))

tappassigtrans <- excel_sheets(paste0(tabroot_dir, "/DifferentialTransExpression_Analysis.xlsx")) %>%  set_names() %>% purrr::map(read_excel, path = paste0(tabroot_dir, "/DifferentialTransExpression_Analysis.xlsx"))

#### DIU - Differential Isoform Usage #########
tappasDIU <- excel_sheets(paste0(tabroot_dir, "/DifferentialTranscriptUsage_Analysis.xlsx")) %>%  set_names() %>% purrr::map(read_excel, path = paste0(tabroot_dir, "/DifferentialTranscriptUsage_Analysis.xlsx")) 
# only keep the genes that show significant differential transcript usage
tappasDIU <- lapply(tappasDIU, function(x) x %>% filter(adjPValue < 0.05))

tappasDIU$rnaseq <- read.table(paste0(tabroot_dir,"/tappAS_RNASeq_DIUGene_Transcripts.tsv"), as.is = T, sep = "\t")
colnames(tappasDIU$rnaseq) <- c("gene","description","DIU","adjPValue","podiumChange","Ctrl_2","Ctr_4","Ctrl_6","Ctrl_8","TG_2","TG_4","TG_6","TG_8")

#### FDA - Functional Diversity Analysis #########
# By Presence
#FDA_FP = read.table(paste0(tappasrna_input_dir,"/FDA_FeaturePresence_Genes_presence.tsv"), as.is = T, sep = "\t")
#colnames(FDA_FP) = c("Gene","NMD","repeat","3UTR Motif","5UTR Motif","PAS","uORF","miRNA Binding",
#                     "Domain","Signal","ACT Site","Binding","Coiled","Compbias","Intramem","Motif","PTM","Transmem")
# By Genomic Position 
#FDA_GP = read.table(paste0(tappasrna_input_dir,"/tappAS_FDA_Category_GenomicPosition_Genes_genomic.tsv"), as.is = T, sep = "\t")
#colnames(FDA_GP) = c("Gene","Repeat","3UTR Length","5UTR Length","CDS","PolyA Site","3UTR Motif","5UTR Motif","PAS","uORF","miRNA Binding"#,"Domain","Signal","ACT Site","Binding","Coiled","Compbias","Intramem","Motif","PTM","Transmem")

# FDA by position of the Genes varying by ID
#FDA_lengths = lapply(list.files(path = tappasrna_input_dir, pattern = "FDA.tsv", full.names = T), function(x) read.table(x) %>% `colnames<-`(c("Gene", "Status")))
#names(FDA_lengths) = word(list.files(path = tappasrna_input_dir, pattern = "FDA.tsv"),c(1),sep = fixed("_"))

#### DFI - Differential Feature Inclusion #########
#DFI = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression/tappAS_DFI_Results_Presence.tsv", header = T, as.is = T, sep = "\t")
#colnames(DFI) = c("Gene","Feature","FeatureElements","DFI_Status","Pvalue","MajorIsoformSwitching","FavouredCondition")
#DFI = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression/dfi_matrix.01952066090.tsv", header = T, as.is = T, sep = "\t")
#DFI_counts = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression/dfi_total_features.01952066090.tsv", header = T, as.is = T, sep = "\t")

#DFI_co = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression/tappAS_coDFI_FeatureAssociations_Presence.tsv")
#colnames(DFI_co) = c("FeatureID1","FeatureID2","GeneswithBothDFI","Coinclusion","MutualExclusive")

####  Alternative splicing events by groups ###############################################
# sqanti files subsetted on the groups in whole transcriptome --> "WT_2mos","WT_8mos","TG_2mos","TG_8mos"
# suppa files gneerated according to the groups (using files generated from group_sqanti_dir)
AS_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/AS/"

## Group classification files 
group.class.files = lapply(list.files(path = AS_dir, pattern = "sqantisubset.classification.txt", full.names = T), 
                           function(x)  SQANTI_class_preparation(x,"standard"))
names(group.class.files) = lapply(list.files(path = AS_dir, pattern = "sqantisubset.classification.txt"), function(x) substr(x,1,5))
all.group.class.files = bind_rows(group.class.files, .id = "Dataset")

#### DMP/DMR #########
# Isabel's RRBS data (Output of load: RRBS_completebetas)
#load("/gpfs/mrc0/projects/Research_Project-MRC148213/EmmaW/DementiaMouseRRBSmethylation/rTg4510_RRBSbetasComplete.RData") 
load(file = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP/Updated/rTg4510_RRBSbetasComplete.RData")
RRBS_completebetas = RRBS_completebetas %>% tibble::rownames_to_column(var = "position") %>% mutate(base = word(position,c(2),sep = fixed(":")), chrom = word(position,c(1),sep = fixed(":"))) 

#RRBS_Phenotype= read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP/Tg4510_phenotype_RRBS.csv")
RRBS_Phenotype <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP/Updated/Tg4510_coldata_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
RRBS_Phenotype$Age_months <- as.factor(RRBS_Phenotype$Age_months)
RRBS_Phenotype$Genotype <- as.factor(RRBS_Phenotype$Genotype)
RRBS_Phenotype$Genotype <- relevel(RRBS_Phenotype$Genotype, "WT")

# Differentially Methylated Positions or Regions for rTg4510 
# Read in files downloaded from RRBS_Paper dropbox folder (Isabel's analysis)
DMPs = list.files(path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP",pattern = "DMPsInteractionModel", full.names = T)
Whole_DMP = lapply(DMPs, function(x) read.table(x, sep = ",", as.is = T, header = T)) 
names(Whole_DMP) = lapply(list.files(path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP",pattern = "DMPsInteractionModel"), function(x) word(x,c(4), sep = fixed("_")))

DMRs = read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP/DMRsGenotype_annotated.csv")

#### WGCNA #########
WGCNA <- excel_sheets(paste0(fig_dir,"/WGCNA.xlsx")) %>% set_names() %>% purrr::map(read_excel, path = paste0(fig_dir, "/WGCNA.xlsx")) %>% lapply(., function(x) x %>% `colnames<-`(c("Gene", "Connectivity", "Module_Membership")))

                      