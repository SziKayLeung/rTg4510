# Szi Kay Leung: sl693@exeter.ac.uk
# Global files for analysis 

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")

# Input Sequencing
targetedpheno <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/Targeted_Sample_Demographics.csv") 
tg4510_samples <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Sample_Description/Tg4510_fullsample.csv")[,c("Genotype","Sample.ID","RIN","ng.ul")]

# Input IsoSeq Paths
CCS_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/IsoSeq/CCS"
LiMA_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/IsoSeq/LIMA"
REFINE_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/IsoSeq/REFINE"
CLUSTER_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/IsoSeq/CLUSTER"
CLUSTER_Merge <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/All_Targeted_Merged/CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv")


# Target Genes and Groups for plots
TargetGene <- str_to_title(c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF"))
ADReg_Genes = c("App","Mapt","Fyn","Trpa1","Vgf")
GWAS_Genes = c("Apoe","Abca1","Abca7","Bin1","Cd33","Picalm","Ptk2b","Sorl1","Trem2")
FTD_Genes = c("Snca","Fus","Tardbp")
EWAS_Genes = c("Ank1","Rhbdf2")


Mapping_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/Post_IsoSeq/MAPPING"
Probes_input <- list.files(path = Mapping_dir, pattern = "fasta.sam.probe_hit.txt", full.names = T)
Probes_files <- lapply(Probes_input, function(x) read.table(x, header=T, as.is=T, sep="\t"))
names(Probes_files ) <- list.files(path = Mapping_dir, pattern = "fasta.sam.probe_hit.txt")

Merged_Mapping_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/All_Targeted_Merged/MAPPING/"
Merged_probes_files <- read.table(paste0(Merged_Mapping_dir,"All_Targeted_Merged.fasta.sam.probe_hit.txt"),header=T, as.is=T, sep="\t")

barcode <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/Barcode_Configs/"
All_Targeted_Merged <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/All_Targeted_Merged"

# SQANTI file 
subsettargeted.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SUBSET/SQANTI3_noRNASEQ/SubsetAllMouseTargeted.collapsed_classification.filtered_lite_classification.txt" 
subsettargeted.class.files <- SQANTI_class_preparation(subsettargeted.class.names.files,"standard") %>% 
  mutate(ss_category = paste0(structural_category,"_", subcategory)) %>%
  filter(ss_category != "ISM_3prime_fragment") %>%
  mutate(ADGene = ifelse(associated_gene %in% TargetGene, "Target Genes","Not Target Genes"))

IsoSeq_rootdir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ"
IsoSeq_filtered_class =  SQANTI_class_preparation(paste0(IsoSeq_rootdir,"/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"),"standard") %>% filter(associated_gene %in% TargetGene) %>% mutate(Dataset = "Iso-Seq") %>%  
  mutate(within_50_cage = ifelse(abs(dist_to_cage_peak) <= 50, "Within 50bp","Not within 50bp"),
         RNASeq_supported = ifelse(min_cov <= 1, "Supported","Not Supported"))

# SQANTI, TAMA filtered file
targeted.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.classification.txt"
targeted.gtf.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/All_Targeted_Merged/Alternative_Pipeline/SQANTI_TAMA_FILTER/All_Targeted_Merged_sqantitamafiltered.classification.gtf"
targeted.class.files <- SQANTI_class_preparation(targeted.class.names.files,"standard")

# SQANTI, TAMA further collapsed files
#collapsed.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.classification.txt"
#coll.class.files <- SQANTI_class_preparation(collapsed.class.names.files,"standard")

# SQANTI filtered files
targeted_sqfil_reason <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI2/AllMouseTargeted.collapsed_classification.filtered_lite_reasons.txt", sep = ",", header = T) %>%  mutate(PBId = paste0("PB.",word(.$filtered_isoform, c(2), sep = fixed("."))))

diff_root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse"

targeted.class.names.files_nornaseq <- paste0(diff_root, "/DiffAnalysis_noRNASEQ/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt")
targeted.preclass.files_nornaseq <- SQANTI_class_preparation(targeted.class.names.files_nornaseq,"standard")

targeted.class.names.files_rnaseq <- paste0(diff_root, "/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt")
targeted.preclass.files_rnaseq <- SQANTI_class_preparation(targeted.class.names.files,"standard")

# Input mapping paths
mapping_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/MAP/PAF"

# PAF mapping stats 
map <- read.table(paste0(mapping_input_dir,"/AllMouseTargeted_reads_with_alignment_statistics.txt"))
colnames(map) <- c("name_of_read","chrom","start","read_length","alignment_length","alignment_length_perc","length_of_target_seq_alignment","alignment_identity","strand","num_mismatches","num_insertions","num_deletions")


#whole.gtf.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.gtf"
#whole.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.txt"
#whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

# Comparison from gffcompare
cuff_tmap = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SUBSET/SQANTI3_noRNASEQ/Whole_Targeted.SubsetAllMouseTargeted.collapsed_classification.filtered_lite.gtf.tmap", header = T)

whole.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_noRNASEQ/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt"
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard") %>% 
  mutate(ss_category = paste0(structural_category,"_", subcategory)) %>%
  filter(ss_category != "ISM_3prime_fragment") %>% 
  mutate(ADGene = ifelse(associated_gene %in% TargetGene, "Target Genes","Not Target Genes")) %>%
  mutate(within_50_cage = ifelse(abs(dist_to_cage_peak) <= 50, "Within 50bp","Not within 50bp"))

# Input Human MAPT Paths
humanmapt_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/HUMANMAPT"


# Human MAPT inputs
# headers of the clustered read names that contain human MAPT sequences
# Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
hMAPT_input <- read.table(paste0(humanmapt_input_dir,"/hMAPT_counts.txt")) %>% `colnames<-`(c("Sample", "Counts"))
mMAPT_input <- read.table(paste0(humanmapt_input_dir,"/hMAPT_counts.txt")) %>% `colnames<-`(c("Sample", "Counts"))
cluster_reads <- list.files(CLUSTER_input_dir, pattern = "cluster_report.csv", full.names = T)
hMAPT_input_all <- read.table(paste0(humanmapt_input_dir,"/AllMouseTargeted_Ids.txt"))


Probe_file <- function(){
  colnames(Merged_probes_files)[6] <- "Probes"
  Merged_probes_files$Gene <- word(Merged_probes_files$Probes,c(3),  sep = fixed ('_'))
  Merged_probes_files$Gene <- word(Merged_probes_files$Gene,c(1),  sep = fixed ('('))
  Merged_probes_files$Gene[Merged_probes_files$Gene == "27642461"] <- "MAPT"
  Merged_probes_files$Gene[Merged_probes_files$Gene == "27642460"] <- "ANK1"
  
  Merged_probes_files <<- Merged_probes_files
}

Filtered_data_staged <- function(){
  Merged_PostIso <- list(read.csv(paste0(All_Targeted_Merged,"/CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv")), # Cluster
                         read.table(paste0(All_Targeted_Merged,"/TOFU/All_Targeted_Merged.collapsed.abundance.txt"), header = T), # Cupcake collapse
                         read.table(paste0(All_Targeted_Merged,"/TOFU/All_Targeted_Merged.collapsed.read_stat.txt"), header = T), 
                         read.table(paste0(All_Targeted_Merged,"/TOFU/All_Targeted_Merged.collapsed.filtered.abundance.txt"), as.is = T, sep = "\t", header = T), 
                         read.table(paste0(All_Targeted_Merged,"/SQANTI/All_Targeted_Merged.collapsed.filtered_classification.filtered_lite_classification.txt"),
                                    as.is = T, sep = "\t", header = T)
  )
  
  names(Merged_PostIso) <- c("cluster","cupcake","cupcake_readstat","cupcake_filter","sqanti_filter")
  Merged_PostIso <<- Merged_PostIso
}

#### TAPPAS (Differential Analysis) ###################
tappasiso_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAPPAS/Results/IsoSeq"
tappasrna_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAPPAS/Results/RNASeq"
tappas_phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_PhenotypeTAPPAS.txt", header = T) %>% mutate(col_names = paste0(group,".",sample))

# Output from Tappas_DEA.R
tappassig <- excel_sheets("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx") %>%  set_names() %>%
  purrr::map(read_excel, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Tables4Figures/DifferentialGeneExpression_Analysis.xlsx")




