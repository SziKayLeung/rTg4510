# Szi Kay Leung: sl693@exeter.ac.uk
# Global files for analysis 

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")

# Phenotypes
twomos <- c("K18","O18","S18", "K17","M21","Q21")
WT <- c("K17","M21","Q21","K23","O23","S23")

# Input Sequencing
sequencing_output <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Sample_Description/Sequenced_mouse_output.csv")
tg4510_samples <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Sample_Description/Tg4510_fullsample.csv")[,c("Genotype","Sample.ID","RIN","ng.ul")]
sequenced <- merge(sequencing_output, tg4510_samples, by.x = "Sample", by.y = "Sample.ID", all.x = TRUE)
sequenced$Genotype <- factor(sequenced$Genotype, levels=c("WT", "TG"))

# Input IsoSeq Paths
CCS_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq/CCS"
LiMA_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq/LIMA"
REFINE_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq/REFINE"
CLUSTER_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq/CLUSTER"

# Input Lengths Paths
lengths_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq/CCS/Lengths"

# Input mapping paths
mapping_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/MAP"
all_mapping_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/MAPPING"

# Input Human MAPT Paths
humanmapt_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/HUMANMAPT"

# sample run 
sample_run <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/Mouse/All_Tg4510_Demultiplex.csv")

# SQANTI, TAMA filtered file
class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.txt"
class.files <- SQANTI_class_preparation(class.names.files,"standard")

# Mouse ERCC
root <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/"
class_suffix <- "_sqantitamafiltered.classification.txt"
ERCC_sqanti_dir <- paste0(root,"Whole_Transcriptome/All_Tg4510/ERCC/WholeIsoSeq")
mouse_lnc_sqanti_dir <- paste0(root,"Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/LNCRNA/WholeIsoSeq")

# Human MAPT inputs
# headers of the clustered read names that contain human MAPT sequences
# Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
hMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "hMAPT2.header", full.names = T)
mMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "mMAPT1.header", full.names = T)
cluster_reads <- list.files(CLUSTER_input_dir, pattern = "cluster_report.csv", full.names = T)
hMAPT_input_all <- read.table(paste0(humanmapt_input_dir,"/WholeIsoSeq_Ids.txt"))

# PAF mapping stats 
map <- read.table(paste0(mapping_input_dir,"/WholeIsoSeq_reads_with_alignment_statistics.txt"))
colnames(map) <- c("name_of_read","chrom","start","read_length","alignment_length","alignment_length_perc","length_of_target_seq_alignment","alignment_identity","strand","num_mismatches","num_insertions","num_deletions")


ERCC_sqanti_files <- function(){
  Ercc.class.file <- read.table(paste0(ERCC_sqanti_dir,class_suffix), header=T, as.is=T, sep="\t")
  Ercc.class.file <<- Ercc.class.file
}

# RNASeq Stringtie Transcriptome
rnaseq_transcriptome <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/rnaseq_stringtie_merged_final_classification.filtered_lite_classification.txt"
cuffrefmap_input <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.refmap"
cufftmap_input <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.tmap"

rnaseq_sqanti_files <- function(){
  rnaseq.class.files <- SQANTI_class_preparation(rnaseq_transcriptome,"rnaseq") %>% mutate(associated_gene = toupper(.$associated_gene))
  rnaseq.class.files$Sample <- "RNA-Seq"
  rnaseq.class.files <<- rnaseq.class.files
}

# INPUT: lncRNA Classification files 
lncrna_class_files <- function(){
  mouse.lnc.file <- paste0(mouse_lnc_sqanti_dir, class_suffix)
  lnc.files <- SQANTI_class_preparation(mouse.lnc.file ,"standard")
  lnc.files <<- lnc.files
}


# RAREFACTION at gene level and isoform level 
rarefaction_files <- function(){
  # list input directory of rarefaction files (prepared)
  Mouse_rarefaction_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/RAREFACTION/WholeIsoSeq"
  
  all_rarefaction_genes <- read.table(paste0(Mouse_rarefaction_dir,'.rarefaction.by_refgene.min_fl_2.txt'),sep=' ',header=T,skip=1) %>% mutate(type = "Genes")
  all_rarefaction_isoforms <- read.table(paste0(Mouse_rarefaction_dir,'.rarefaction.by_refisoform.min_fl_2.txt'),sep=' ',header=T,skip=1) %>% mutate(type = "Isoforms")

  all_rarefaction_isoforms_category <- read.table(paste0(Mouse_rarefaction_dir,'.rarefaction.by_refisoform.min_fl_2.by_category.txt'),sep=' ',header=T,skip=1)
  all_rarefaction_isoforms_category$category <- factor(all_rarefaction_isoforms_category$category, 
                                                            levels = c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", 
                                                                       "novel_not_in_catalog", "fusion", "antisense", "genic", "intergenic"),
                                                            labels = c("FSM", "ISM", "NIC", "NNC", "Fusion", "Antisense", "Genic", "Intergenic"))
  
  for(j in 1:nrow(all_rarefaction_isoforms_category)){
    all_rarefaction_isoforms_category$type[j] <- if(all_rarefaction_isoforms_category$category[j] %in% c("FSM", "ISM")){
      "Annotated" 
    } else if (all_rarefaction_isoforms_category$category[j] %in% c("NIC", "NNC")){
      "Novel"
    } else {
      "Others"
    }
  }

  
  all_rarefaction_genes <<- all_rarefaction_genes
  all_rarefaction_isoforms <<- all_rarefaction_isoforms
  all_rarefaction_isoforms_category <<- all_rarefaction_isoforms_category 
}

##### FeatureCounts 
# FeatureCounts output from RNA2IsoSeq (aligned RNASeq to Iso-Seq defined transcriptome)
#IsoSeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsovsRnaseq/IsoSeq_Def/RNA2IsoSeq.transcript_id.tsv", header = T)
# FeatureCounts output from RNA2RNASeq (aligned RNASeq to RNA-Seq defined transcriptome)
#RNASeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsovsRnaseq/RNASeq_Def/RNA2RNASeq.transcript_id.tsv", header = T)
##### Kallisto
RNASeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/KALLISTO/WholeIsoSeq.abundance.tsv", header = T) 
IsoSeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/KALLISTO/WholeIsoSeq.abundance.tsv", header = T) %>% mutate(target_id = word(target_id, c(1), sep = fixed("<")))

# RNASeq Stringtie Transcriptome
rnaseq_transcriptome <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/rnaseq_stringtie_merged_final_classification.filtered_lite_classification.txt"
cuffrefmap_input <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.refmap"
cufftmap_input <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.tmap"