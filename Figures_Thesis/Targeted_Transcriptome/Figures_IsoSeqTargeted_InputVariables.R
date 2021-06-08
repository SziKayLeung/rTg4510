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
TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")
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

# SQANTI, TAMA filtered file
targeted.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt"
targeted.gtf.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/OLD2021/All_Targeted_Merged/Alternative_Pipeline/SQANTI_TAMA_FILTER/All_Targeted_Merged_sqantitamafiltered.classification.gtf"
targeted.class.files <- SQANTI_class_preparation(targeted.class.names.files,"standard")


# SQANTI filtered files
targeted_sqfil_reason <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI2/AllMouseTargeted.collapsed_classification.filtered_lite_reasons.txt", sep = ",", header = T) %>%  mutate(PBId = paste0("PB.",word(.$filtered_isoform, c(2), sep = fixed("."))))
targeted.preclass.names.files <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI2/AllMouseTargeted.collapsed_classification.txt"
targeted.preclass.files <- SQANTI_class_preparation(targeted.preclass.names.files,"standard")

targeted.preclass.names.files.final <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI2/TAMAFILTERSQANTI/AllMouseTargeted_sqantitamafiltered.classification.txt"
targeted.preclass.files.final <- SQANTI_class_preparation(targeted.preclass.names.files.final,"standard")

# Input mapping paths
mapping_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/MAP/PAF"

# PAF mapping stats 
map <- read.table(paste0(mapping_input_dir,"/AllMouseTargeted_reads_with_alignment_statistics.txt"))
colnames(map) <- c("name_of_read","chrom","start","read_length","alignment_length","alignment_length_perc","length_of_target_seq_alignment","alignment_identity","strand","num_mismatches","num_insertions","num_deletions")


whole.gtf.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.gtf"
whole.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.txt"
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

# output file from TAMA merge 
TAMA_transfile <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Whole_vs_Targeted/merged_whole_targeted_trans_report.txt", header = T) %>% mutate(gene_id = word(transcript_id, c(1), sep = fixed("."))) %>% mutate(gene_name = word(word(all_source_trans, c(3), sep = fixed("_")),c(1), sep = ","))

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


