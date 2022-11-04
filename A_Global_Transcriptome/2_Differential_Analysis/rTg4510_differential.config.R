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


output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Whole_Transcriptome"


## ---------- TappAS input files -----------------

# input directory
WKD_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/DiffAnalysis_Final/"
TAPPAS_PHENOTYPE_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/A_isoseq_whole/Tappas/"
TAPPAS_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/DiffAnalysis_Final/TAPPAS_OUTPUT"

# output files from running tappAS
TAPPAS_INPUT_DIR = list(
  iso = paste0(WKD_ROOT,"TAPPAS_OUTPUT/IsoSeq_Expression"),
  rna = paste0(WKD_ROOT,"TAPPAS_OUTPUT/RNASeq_Expression")
)

# phenotype 
TAPPAS_PHENOTYPE = list(
  iso = paste0(TAPPAS_PHENOTYPE_DIR,"WholeIsoSeq_PhenotypeTAPPAS.txt"),
  rna = paste0(TAPPAS_PHENOTYPE_DIR,"WholeAllMouse_PhenotypeTAPPAS.txt")
)
phenotype <- lapply(TAPPAS_PHENOTYPE, function(x) read.table(x, header = T))


## ---------- SQANTI classification files -----------------

# Classification file
class.names.files = paste0(WKD_ROOT, "SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt")
class.files <- SQANTI_class_preparation(class.names.files,"standard")


## ---------- TappAS output -----------------

tappassigtrans <- excel_sheets(paste0(TAPPAS_ROOT, "/DifferentialTransExpression_Analysis.xlsx")) %>% 
  set_names() %>% purrr::map(read_excel, path = paste0(TAPPAS_ROOT, "/DifferentialTransExpression_Analysis.xlsx"))

# DIU - Differential Isoform Usage 
tappasDIU <- excel_sheets(paste0(TAPPAS_ROOT, "/DifferentialTranscriptUsage_Analysis.xlsx")) %>%  
  set_names() %>% purrr::map(read_excel, path = paste0(TAPPAS_ROOT, "/DifferentialTranscriptUsage_Analysis.xlsx")) 
# only keep the genes that show significant differential transcript usage
tappasDIU <- lapply(tappasDIU, function(x) x %>% filter(adjPValue < 0.05))

tappasDIU$rnaseq <- read.table(paste0(TAPPAS_ROOT,"/tappAS_RNASeq_DIUGene_Transcripts.tsv"), as.is = T, sep = "\t")
colnames(tappasDIU$rnaseq) <- c("gene","description","DIU","adjPValue","podiumChange","Ctrl_2","Ctr_4","Ctrl_6","Ctrl_8","TG_2","TG_4","TG_6","TG_8")


