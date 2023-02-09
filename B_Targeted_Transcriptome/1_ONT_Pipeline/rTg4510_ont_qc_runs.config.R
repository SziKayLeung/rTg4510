## ---------- Script -----------------
##
## Purpose: input variables for characterising QC of ONT targeted mouse transcriptome datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- Packages -----------------

suppressMessages(library("dplyr"))


## ---------- LOGEN modules -----------------

LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "longread_QC/number_ont_reads.R"))


## ---------- Directory and input files -----------------

dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/",
  raw = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/1_raw/F_ont_targeted/",
  meta = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/F_ont_targeted/"
)


# Miscellaneous input files 
misc_input <- list(
  
  # list of samples
  tg4510_samples = read.csv(paste0(dirnames$meta, "Tg4510_fullsample.csv"))[,c("Genotype","Sample.ID","RIN","ng.ul")],
  
  # list of target genes
  TargetGene = read.table(paste0(dirnames$meta, "TargetGenes.tsv"))[["V1"]],
  
  # list of barcodes for each sample
  BarcodedPhenotype = read.csv(paste0(dirnames$meta, "ONTBarcoded_Phenotype.csv")),
  
  # Batch 2 sequencing summary
  B2_seq_summary  = paste0(dirnames$raw, "sequencing_summary_FAO06462_2abef59d.txt"), 
  
  # Batch 3 sequencing summary 
  B3_seq_summary  = paste0(dirnames$raw, "/sequencing_summary_FAO06635_ecb1ad97.txt")#,
  
  # read counts 
  #ReadCounts = read_batch_readcounts(dirnames$raw, dirnames$root, 2, 3)

)

samples_removed = c("BC10")
num_barcodes = 10


## ---------- Sequencing raw files -----------------

sequencing_data <- list(
  B2 = prepare_summary_file(misc_input$B2_seq_summary), 
  B3 = prepare_summary_file(misc_input$B3_seq_summary)
)

# Summarise the number of sequencing counts in batch 2 and 3 (sequenced and basecalled)
Sequencing_Counts <- batch_summarise_counts(sequencing_data)


## ---------- Target Rate -----------------

Probes_input <- list(
  Batch2 = read_target_probes(paste0(dirnames$root,"/4b_target_rate/Batch2"), "fasta.sam.probe_hit.txt"),
  Batch3 = read_target_probes(paste0(dirnames$root,"/4b_target_rate/Batch3"), "fasta.sam.probe_hit.txt")
)
Probes_input <- unlist(Probes_input, recursive=FALSE)