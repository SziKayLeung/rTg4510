# Szi Kay Leung: sl693@exeter.ac.uk
# Global files for analysis 

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
root_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/"

# SQANTI filtered file 
sqanti.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"
sqanti.class.files <- SQANTI_class_preparation(sqanti.names.files,"standard")

# SQANTI, TAMA filtered file
class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt"
class.files <- SQANTI_class_preparation(class.names.files,"standard")

# SQANTI, TAMA further collapsed files
collapsed.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/COLLAPSE_FILTER/Length/AllMouseTargeted_sqantisubset.classification.txt"
coll.class.files <- SQANTI_class_preparation(collapsed.class.names.files,"standard")

sqanti.filtered.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"
sqantifil.class.files <- read.table(sqanti.filtered.names.files, header = T)

# SQANTI file of merged whole and targeted transcriptome 
merged.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/Whole_Targeted/merged_whole_targeted_classification.txt"
merged.class.files <- read.table(merged.names.files, header = T)

# Tama-merged corresponding isoforms
merged_targetedid <- read.table(paste0(root_dir,"Whole_Targeted/merged_targetedPBID.txt"), header = T)
merged_wholeid <- read.table(paste0(root_dir,"Whole_Targeted/merged_wholePBID.txt"), header = T)

#### TAPPAS (Differential Analysis) ###################
tappasiso_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_OUTPUT/IsoSeq_Expression"
tappasrna_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression"
tappasisocol_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_OUTPUT/IsoSeq_Expression_Collapsed"
tappasrnacol_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression_Collapsed"
tappasrnacolexp_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression_Collapsed_byFSMEXp"
tappasrnamerged_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression_WholeTargeted"
tappasiso_phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_PhenotypeTAPPAS.txt", header = T) %>% mutate(variable = paste0(sample))
tappasrna_phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_RNASeqPhenotypeTAPPAS.txt", header = T) %>% mutate(col_names = paste0(group,".",sample),age = paste0(time, "_mos"), pheno = ifelse(group == "CONTROL", "WT", "TG"), variable = sample)

# Output from Tappas_DEA.R
tappassiggene <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx") %>%  set_names() %>% purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx")

tappassigtrans <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTransExpression_Analysis.xlsx") %>%  set_names() %>% purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTransExpression_Analysis.xlsx")

#tappasDIU <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTranscriptUsage_Analysis.xlsx") %>%  set_names() %>% purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTranscriptUsage_Analysis.xlsx")
# only keep the genes that show significant differential transcript usage
#tappasDIU <- lapply(tappasDIU, function(x) x %>% filter(adjPValue < 0.05))

## Group classification files 
#group.class.files.dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/GROUPS"
#group.class.files = lapply(list.files(path = group.class.files.dir, pattern = "sqantitamafiltered.classification.txt", full.names = T), 
#                           function(x)  SQANTI_class_preparation(x,"standard"))
#names(group.class.files) = lapply(list.files(path = group.class.files.dir, pattern = "sqantitamafiltered.classification.txt"), function(x) substr(x,1,5))
#all.group.class.files = bind_rows(group.class.files, .id = "Dataset")
