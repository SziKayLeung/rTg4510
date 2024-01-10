# Szi Kay Leung: sl693@exeter.ac.uk

suppressMessages(library("data.table"))
LOGEN <- "/lustre/projects/Research_Project-MRC148213/sl693/scripts/LOGen"
source(paste0(LOGEN,"/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"/target_gene_annotation/summarise_gene_stats.R"))
source(paste0(LOGEN,"/compare_datasets/dataset_identifer.R"))
source(paste0(LOGEN, "/differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN, "/aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "/differential_analysis/run_DESeq2.R"))
source(paste0(LOGEN, "/merge_characterise_dataset/run_ggtranscript.R"))
source(paste0(LOGEN, "/aesthetics_basics_plots/pthemes.R"))


wholesamples <- c("K17","K18","K23","K24","L22","M21","O18","O23","Q20","Q21","S18","S23")
wholeTG <- c("K18","K24", "L22","O18","Q20","S18")
wholeWT <- setdiff(wholesamples, wholeTG)
whole2mos <- c("K17","M21","Q21","K18","O18","S18")
whole8mos <- setdiff(wholesamples, whole2mos)

targetedWT <- c("K19","K23","K21","K17","S19","M21","O23","P19","Q21","S23","Q17","Q23")
targetedTG <- c("K18","K20","K24","L22","O18","O22","T20","Q20","S18","Q18","L18","T18")

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
  protein = paste0(root_dir, "rTg4510/G_Merged_Targeted/4_proteogenomics/"),
  
  # reference
  references = paste0(root_dir,"reference/annotation")
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
  cpat = read.table(paste0(dirnames$protein,"5_calledOrfs/all_iso_ont_best_orf.tsv"), sep ="\t", header = T),
  t2p.collapse = read.table(paste0(dirnames$protein,"6_refined_database/all_iso_ont_orf_refined.tsv"), sep = "\t", header = T)
)


## --------------------------
# gtf
gtf = list(
  targ_merged = rtracklayer::import(paste0(dirnames$targ_root,"/2_sqanti3/all_iso_ont_collapsed.filtered_counts_filtered.gtf")),
  ptarg_merged = rtracklayer::import(paste0(dirnames$targ_root,"/4_proteogenomics/5_calledOrfs/all_iso_ont.gtf")),
  ref_target = rtracklayer::import(paste0(dirnames$references,"/gencode.M22.annotation.20Targets.gtf"))
)
gtf <- lapply(gtf, function(x) as.data.frame(x))
gtf$ptarg_merged <- gtf$ptarg_merged %>% mutate(transcript = word(transcript_id,c(1),sep=fixed("_")))

gtf$targ_merged <- rbind(gtf$targ_merged[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                         gtf$ptarg_merged[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                         gtf$ref_target[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])


## -------- FICLE output -------------------
#Maptprotein <- unique(class.files$ptarg_filtered[class.files$ptarg_filtered$associated_gene == "Mapt","corrected_acc"])
MaptES <- read.csv("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510_FICLE/Mapt/Stats/Mapt_general_exon_level.csv")
Trem2ES <- read.csv("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510_FICLE/FICLE/TargetGenes/Trem2/Stats/Trem2_Exonskipping_generaltab.csv") 
Trem2NE <- read.csv("/lustre/projects/Research_Project-MRC148213/sl693/rTg4510_FICLE/FICLE/TargetGenes/Trem2/Stats/Trem2_NE.csv") 
