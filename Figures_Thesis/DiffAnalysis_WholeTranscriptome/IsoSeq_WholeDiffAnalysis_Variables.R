# Szi Kay Leung: sl693@exeter.ac.uk
# Global files for analysis 

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")


# SQANTI, TAMA filtered file
class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantitamafiltered.classification.txt"
class.files <- SQANTI_class_preparation(class.names.files,"standard")


#### TAPPAS (Differential Analysis) ###################
tappasiso_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/IsoSeq_Expression"
tappasrna_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/RNASeq_Expression"
tappasiso_phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome/WholeIsoSeq_PhenotypeTAPPAS.txt", header = T) %>% mutate(variable = paste0(sample))
tappasrna_phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome/WholeAllMouse_PhenotypeTAPPAS.txt", header = T) %>% mutate(col_names = paste0(group,".",sample),age = paste0(time, "_mos"), pheno = ifelse(group == "CONTROL", "WT", "TG"), variable = paste0("FL.", sample,"_",pheno,"_",age))

# Isabel's supplementary table of differentially expressed genes in rTg4510 
Isabel_gene_Tg4510AgeGenotypeDEG = read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/Isabel_Supp4_Tg4510AgeGenotypeDEG.csv", header = T, as.is = T)
Isabel_gene_Tg4510GenotypeDEG = read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/Isabel_Supp2_Tg4510GenotypeDEG.csv", header = T, as.is = T)

# Output from Tappas_DEA.R
tappassiggene <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx") %>%  set_names() %>% purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx")

tappassigtrans <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTransExpression_Analysis.xlsx") %>%  set_names() %>% purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTransExpression_Analysis.xlsx")

tappasDIU <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTranscriptUsage_Analysis.xlsx") %>%  set_names() %>% purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialTranscriptUsage_Analysis.xlsx")
# only keep the genes that show significant differential transcript usage
tappasDIU <- lapply(tappasDIU, function(x) x %>% filter(adjPValue < 0.05))

## Alternative splicing events by groups
# sqanti files subsetted on the groups in whole transcriptome --> "WT_2mos","WT_8mos","TG_2mos","TG_8mos"
group_sqanti_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/GROUPS/"
# suppa files gneerated according to the groups (using files generated from group_sqanti_dir)
suppa_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/AS"

## Group classification files 
group.class.files.dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/GROUPS"
group.class.files = lapply(list.files(path = group.class.files.dir, pattern = "sqantitamafiltered.classification.txt", full.names = T), 
                           function(x)  SQANTI_class_preparation(x,"standard"))
names(group.class.files) = lapply(list.files(path = group.class.files.dir, pattern = "sqantitamafiltered.classification.txt"), function(x) substr(x,1,5))
all.group.class.files = bind_rows(group.class.files, .id = "Dataset")
