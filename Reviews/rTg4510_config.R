# Szi Kay Leung: sl693@exeter.ac.uk

suppressMessages(library("data.table"))
LOGEN <- "/lustre/projects/Research_Project-MRC148213/sl693/scripts/LOGen"
source(paste0(LOGEN,"/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"/target_gene_annotation/summarise_gene_stats.R"))
source(paste0(LOGEN,"/compare_datasets/dataset_identifer.R"))
source(paste0(LOGEN, "/differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN, "/aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "/differential_analysis/run_DESeq2.R"))


TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                        "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")

## --------------------------- 
# directory names
root_dir <- "/lustre/projects/Research_Project-MRC148213/sl693/"
dirnames <- list(
  # global transcriptome (Iso-Seq, Iso-Seq + RNA-Seq)
  glob_root = paste0(root_dir, "rTg4510/A_IsoSeq_Whole"),
  targ_root = paste0(root_dir, "rTg4510/G_Merged_Targeted"),
  glob_output = paste0(root_dir, "rTg4510/01_figures_tables/Whole_Transcriptome"),
  targ_output = paste0(root_dir, "rTg4510/01_figures_tables/Targeted_Transcriptome"),
  
  # proteogeonomics
  protein = paste0(root_dir, "rTg4510/G_Merged_Targeted/4_proteogenomics/6_refined_database")
)

## ------------- Phenotype files -------------------

phenotype  <- list(
  targ_ont = read.table(paste0(root_dir, "rTg4510/0_metadata/TargetedOntPhenotype.txt"), header = T), 
  targ_iso = read.table(paste0(root_dir, "rTg4510/0_metadata/TargetedIsoSeqPhenotype.txt"), header = T)
)


## --------------------------- 
# Final classification file
class.names.files <- list(
  glob_iso = paste0(dirnames$glob_root, "/2_sqanti3/WholeIsoSeq.collapsed_RulesFilter_result_classification.txt"),
  targ_offtargets = paste0(dirnames$targ_root, "/2_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.txt"),
  targ_all = paste0(dirnames$targ_root, "/2_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts.txt"),
  targ_filtered = paste0(dirnames$targ_root, "/2_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt")
) 
class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))


## ---------------------------

# Expression 
rawExp <- list(
  targ_ont_all = read.csv(paste0(dirnames$targ_root, "/1_cupcake_collapse/demux_fl_count.csv"))
)

## ---------- DESeq2 results -----------------
cat("Input DESEQ2 results\n")
GlobalDESeq <- list(
  resTranAnno = readRDS(file = paste0(dirnames$glob_output, "/IsoSeq_DESeq2TranscriptLevel.RDS")),
  resGeneAnno = readRDS(file = paste0(dirnames$glob_output, "/IsoSeq_DESeq2GeneLevel.RDS"))
) 

TargetedDESeq <- list(
  ontResTranAnno = readRDS(file = paste0(dirnames$targ_output, "/Ont_DESeq2TranscriptLevel.RDS")),
  isoResTranAnno = readRDS(file = paste0(dirnames$targ_output, "/IsoSeq_DESeq2TranscriptLevel.RDS")),
  ontResGeneAnno = readRDS(file = paste0(dirnames$targ_output, "/Ont_DESeq2GeneLevel.RDS")),
  isoResGeneAnno = readRDS(file = paste0(dirnames$targ_output, "/IsoSeq_DESeq2GeneLevel.RDS"))
) 
TargetedDESeq$ontResTranAnno$wald$anno_res <- TargetedDESeq$ontResTranAnno$wald$anno_res %>% filter(padj < 0.05)
TargetedDESeq$ontResTranAnno$wald8mos$anno_res <- TargetedDESeq$ontResTranAnno$wald8mos$anno_res %>% filter(padj < 0.05)
TargetedDESeq$ontResTranAnno$wald$norm_counts <- merge(TargetedDESeq$ontResTranAnno$wald$norm_counts, class.files$targ_filtered[,c("isoform","structural_category")], by = "isoform")


## -------- Proteogenomics -------------------

protein = list(
  t2p.collapse = read.table(paste0(dirnames$protein,"/all_iso_ont_orf_refined.tsv"), sep = "\t", header = T)
)
