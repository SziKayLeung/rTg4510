## ---------- Script -----------------
##
## Purpose: Input variables for rTg4510 Iso-Seq global dataset 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk_
##


## ---------- Directory and input files -----------------
dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/",
  meta = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/A_isoseq_whole/",
  scripts = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome/1_IsoSeq_Pipeline/",
  P2021 = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/2_post_isoseq3/9_paper2021/",
  ref = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/references/",
  as_events = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/Whole_Transcriptome_Paper/Output/Tables/AS_IR/"
)


# Phenotypes
twomos <- c("K18","O18","S18", "K17","M21","Q21")
control_samples <- c("K17","M21","Q21","K23","O23","S23")
case_samples <- c("O18","K18","S18","L22","Q20","K24")


# Input IsoSeq Paths
spec_dirnames <- list(
  ccs = paste0(dirnames$root, "1_isoseq3/1_ccs"),
  lima = paste0(dirnames$root, "1_isoseq3/2_lima"),
  refine = paste0(dirnames$root, "1_isoseq3/3_refine"),
  cluster = paste0(dirnames$root, "1_isoseq3/4_cluster"),
  lengths = paste0(dirnames$root, "1_isoseq3/1_ccs/Lengths")
)

# Miscellaneous input files 
misc_input <- list(
  
  # sample run 
  sample_run = read.table(paste0(dirnames$scripts, "/rTg4510_samples.tsv")) %>% `colnames<-`(c("sample", "dir", "run_id")),
  
  # Input Sequencing
  sequencing_output = read.csv(paste0(dirnames$meta, "Sequenced_mouse_output.csv")),
  tg4510_samples = read.csv(paste0(dirnames$meta, "Tg4510_fullsample.csv"))[,c("Genotype","Sample.ID","RIN","ng.ul")],
  
  # cluster csv file
  ccs_output = read.csv(paste0(spec_dirnames$ccs, "/WholeIsoSeqAll_CCS_output.csv"), header = T),
  lima_output = read.csv(paste0(spec_dirnames$lima, "/WholeIsoSeqAll_LIMA_summary.csv"), header = T),
  CLUSTER_merge = read.csv(paste0(dirnames$root, "1_isoseq3/5_merged_cluster/WholeIsoSeq.clustered.cluster_report.csv")),
  
  # Mapping statistics
  map = read.table(paste0(dirnames$root,"2_post_isoseq3/6_minimap/WholeIsoSeq_reads_with_alignment_statistics.txt")) %>%
    `colnames<-`(c("name_of_read","chrom","start","read_length","alignment_length",
                   "alignment_length_perc","length_of_target_seq_alignment","alignment_identity",
                   "strand","num_mismatches","num_insertions","num_deletions")),
  
  # rarefaction 
  rarefaction = read_rarefaction_files(paste0(dirnames$P2021,"RAREFACTION/"), "WholeIsoSeq"),
  
  # kallisto alignment 
  RNASeq_Def = read.table(paste0(dirnames$P2021, "RNASeq/KALLISTO/WholeIsoSeq.abundance.tsv"), header = T), 
  IsoSeq_Def = read.table(paste0(dirnames$P2021, "/SQANTI_TAMA_FILTER/GENOME/KALLISTO/WholeIsoSeq.abundance.tsv"), header = T) %>% 
    mutate(target_id = word(target_id, c(1), sep = fixed("<"))),
  
  # Gffcompare output
  cuffrefmap_input = read.table(paste0(dirnames$P2021, "RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.refmap"),
                                header =T),
  cufftmap_input = read.table(paste0(dirnames$P2021, "RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.tmap"),
                              header = T),
  
  # ERCC calcuations
  ercc_cal = read.table(paste0(dirnames$ref, "ERCC/ERCC_whole_transcriptome_calculations.txt"), as.is = T, sep = "\t", header = T),
  
  # Alternative splicing events
  AS_IR_genes = paste0(dirnames$as_events,"ALL_SUPPA2_Genes_Output_Updated2.csv"),
  AS_IR_knowngenes = paste0(dirnames$as_events,"Mouse_AS_KnownGenes.csv"),
  AS_IR_isoforms = paste0(dirnames$as_events,"Mouse_AS_NovelAnnoIsoforms.csv")
)


# important to have a "Genotype column"
# Total.Bases..GB, Genotype, RIN 
misc_input$sequenced = merge(misc_input$sequencing_output, misc_input$tg4510_samples, by.x = "Sample", by.y = "Sample.ID", all.x = TRUE) %>% 
  mutate(Genotype = factor(Genotype, levels=c("WT", "TG"))) %>% filter(Genotype != "NA")

# sample run 
#sample_run <- read.csv(paste0(dirnames$scripts, "IsoSeq_Processing/Mouse/All_Tg4510_Demultiplex.csv"))


## ---------- SQANTI classification file -----------------

input.class.names.files <- list(
  isoseq = paste0(dirnames$root, "2_post_isoseq3/9_paper2021/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.txt"),
  ercc = paste0(dirnames$P2021, "ERCC/WholeIsoSeq_sqantitamafiltered.classification.txt"),
  lncRNA = paste0(dirnames$P2021, "SQANTI_TAMA_FILTER/LNCRNA/WholeIsoSeq_sqantitamafiltered.classification.txt"),
  rnaseq = paste0(dirnames$P2021, "RNASeq/SQANTI2/rnaseq_stringtie_merged_final_classification.filtered_lite_classification.txt")
)

input.class.files <- list(
  isoseq = SQANTI_class_preparation(input.class.names.files$isoseq,"standard"),
  ercc = read.table(input.class.names.files$ercc, header=T, as.is=T, sep="\t"),
  lncRNA = SQANTI_class_preparation(input.class.names.files$lncRNA,"standard"),
  rnaseq = SQANTI_class_preparation(input.class.names.files$rnaseq,"rnaseq") %>% 
    mutate(associated_gene = toupper(.$associated_gene))
)


## ---------- Transgene  -----------------

# Human MAPT inputs
# headers of the clustered read names that contain human MAPT sequences
# Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
humanmapt_input_dir <- paste0(dirnames$P2021, "HUMANMAPT") 
hMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "hMAPT2.header", full.names = T)
mMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "mMAPT1.header", full.names = T)
cluster_reads <- list.files(spec_dirnames$cluster, pattern = "cluster_report.csv", full.names = T)
hMAPT_input_all <- read.table(paste0(humanmapt_input_dir,"/WholeIsoSeq_Ids.txt"))

