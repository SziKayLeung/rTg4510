# Szi Kay Leung: sl693@exeter.ac.uk
# Global files for analysis 

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")

# root dir
diffroot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ"
rawroot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome"
root_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ"
tabroot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/TAPPAS_OUTPUT"

# Final classification file
class.names.files = paste0(diffroot_dir, "/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt")
class.files <- SQANTI_class_preparation(class.names.files,"standard")

#### TAPPAS (Differential Analysis) ###############################################
tappasiso_input_dir <- paste0(diffroot_dir, "TAPPAS_OUTPUT/IsoSeq_Expression")
tappasiso_phenotype <- read.table(paste0(rawroot_dir, "/TargetedMouse_PhenotypeTAPPAS.txt"), header = T) %>% mutate(variable = paste0(sample))

#### DGE/DTE - Differential Gene and Transcript Expression #########
# Output from Tappas_DEA.R
tappassiggene <- excel_sheets(paste0(tabroot_dir,"/DifferentialGeneExpression_Analysis.xlsx")) %>% set_names() %>% purrr::map(read_excel, path = paste0(tabroot_dir,"/DifferentialGeneExpression_Analysis.xlsx"))

tappassigtrans <- excel_sheets(paste0(tabroot_dir, "/DifferentialTransExpression_Analysis.xlsx")) %>%  set_names() %>% purrr::map(read_excel, path = paste0(tabroot_dir, "/DifferentialTransExpression_Analysis.xlsx"))