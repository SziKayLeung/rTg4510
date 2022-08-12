#### TAPPAS (Differential Analysis) ###################
tappasiso_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAPPAS/Results/IsoSeq"
tappasrna_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAPPAS/Results/RNASeq"
tappas_phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_PhenotypeTAPPAS.txt", header = T) %>% mutate(col_names = paste0(group,".",sample))

# Output from Tappas_DEA.R
tappassig <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx") %>%  set_names() %>%
  purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx")
