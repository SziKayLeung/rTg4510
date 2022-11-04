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


output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"


## ---------- TappAS input files -----------------

# input directory
ISOSEQ_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/thesis_dump/DiffAnalysis_noRNASEQ/"
ONT_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/thesis_dump/TALON/"
MERGED_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/"
TAPPAS_PHENOTYPE_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/"

# output files from running tappAS
TAPPAS_INPUT_DIR = list(
  iso = paste0(ISOSEQ_DIR,"TAPPAS_OUTPUT/IsoSeq_Expression"),
  ont = paste0(ONT_DIR,"TAPPAS_OUTPUT")
)

# phenotype 
TAPPAS_PHENOTYPE = list(
  iso = paste0(TAPPAS_PHENOTYPE_DIR,"B_isoseq_targeted/TargetedMouse_PhenotypeTAPPAS.txt"),
  ont = paste0(TAPPAS_PHENOTYPE_DIR,"F_ont_targeted/ONT_phenotype.txt")
)
phenotype <- lapply(TAPPAS_PHENOTYPE, function(x) read.table(x, header = T))


## ---------- SQANTI classification files -----------------

# Classification file
class.names.files <- list(
  iso = paste0(ISOSEQ_DIR, "SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"),
  ont = paste0(ONT_DIR, "All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon_classification.txt"),
  merged = paste0(MERGED_DIR, "3_sqanti3/IsoSeqONT_final_genename_classification_noISM.txt")
)

class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

## ---------- TappAS output -----------------

read_dea_files <- function(dir){
  f = excel_sheets(dir) %>% purrr::set_names() %>% purrr::map(read_excel, path = dir)
  return(f)
}

tappassiggene <- lapply(list("iso" = paste0(ISOSEQ_DIR,"/TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx"),
                             "ont" = paste0(ONT_DIR,"TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx")), 
                        function(x) read_dea_files(x))

tappassigtrans <- lapply(list("iso" = paste0(ISOSEQ_DIR,"/TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx"),
                              "ont" = paste0(ONT_DIR,"TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx")), 
                         function(x) read_dea_files(x))
