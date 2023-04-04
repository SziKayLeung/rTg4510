## ---------- Script -----------------
##
## Purpose: perform differential analysis on mouse rTg4510 ONT and Iso-Seq targeted datasets using linear regression
## Transcript level separate analysis
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))


## ---------- source functions -----------------

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))


## ---------- input -----------------

# directory names
dirnames <- list(
  rTg4510 = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Targeted_Transcriptome"
)

# read input files
input_files <- list(
  ontPhenotype = paste0(dirnames$rTg4510, "0_metadata/F_ont_targeted/ONT_phenotype.txt"), 
  isoPhenotype = paste0(dirnames$rTg4510, "0_metadata/B_isoseq_targeted/TargetedMouse_PhenotypeTAPPAS.txt"), 
  
  # merged data
  expression = paste0(dirnames$rTg4510, "G_Merged_Targeted/B_cupcake_pipeline/2_collapse/demux_fl_count.csv"),
  classfiles = paste0(dirnames$rTg4510, "G_Merged_Targeted/B_cupcake_pipeline/3_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts.txt")
)


input <- list()
input$ontPhenotype <- read.table(input_files$ontPhenotype, sep = "\t", header = T) %>% mutate(sample = paste0("ONT_",sample))
input$isoPhenotype <- read.table(input_files$isoPhenotype, sep = "\t", header = T) %>% mutate(sample = str_replace(sample, "FL.", "Iso.Seq_"))
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")
input$expression <- read.csv(input_files$expression) %>% .[.$isoform != "0",] %>% filter(isoform %in% input$classfiles$isoform)


# datawrangle for input to run_DESeq2()
# keep only ONT counts or Iso-Seq counts
# keep only targeted genes counts
input$ontExpression <- input$expression %>% select(isoform, contains("ONT"), -ONT_sum_FL)
input$isoExpression <- input$expression %>% select(isoform, contains("Iso.Seq"), -Iso.Seq_sum_FL)
input$ontmos8Phenotype <- input$ontPhenotype %>% filter(time == "8")
input$isomos8Phenotype <- input$isoPhenotype %>% filter(time == "8")

input$gene_expression <- input$expression %>% mutate(associated_gene = word(isoform,c(2), sep = fixed(".")))
input$gene_expression <- aggregate(. ~ associated_gene, input$gene_expression %>% select(-"isoform"), sum)
rownames(input$gene_expression) <- input$gene_expression$associated_gene


## ---------- ONT: Creating DESeq2 object and analysis -----------------

# run DESeq2
ontResTran <- list(
  wald = run_DESeq2(test="Wald",input$ontExpression,input$ontPhenotype,threshold=10,exprowname="isoform",controlname="CONTROL",design="time_series",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$ontExpression,input$ontPhenotype,threshold=10,exprowname="isoform",controlname="CONTROL",design="time_series",interaction="On"),
  wald8mos = run_DESeq2(input$ontExpression,input$ontmos8Phenotype ,exprowname="isoform",controlname="CONTROL",design="case_control",interaction="On",test="Wald")
)

ontResTranAnno <- lapply(ontResTran, function(x) anno_DESeq2(x,input$classfiles,input$ontPhenotype,controlname="CONTROL",level="transcript",sig=0.1))

# split by genotype and age effects
ontResTranEffects <- do.call(rbind,dissect_DESeq2(wald=ontResTranAnno$wald$anno_res,lrt=ontResTranAnno$lrt$anno_res))
ontResTranEffects

## ---------- Iso-Seq: Creating DESeq2 object and analysis -----------------

# run DESeq2
isoResTran <- list(
  wald = run_DESeq2(test="Wald",input$isoExpression,input$isoPhenotype,threshold=10,exprowname="isoform",controlname="CONTROL",design="time_series",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$isoExpression,input$isoPhenotype,threshold=10,exprowname="isoform",controlname="CONTROL",design="time_series",interaction="On"),
  wald8mos = run_DESeq2(input$isoExpression,input$isomos8Phenotype,exprowname="isoform",controlname="CONTROL",design="case_control",interaction="On",test="Wald")
)

isoResTranAnno <- lapply(isoResTran, function(x) anno_DESeq2(x,input$classfiles,input$isoPhenotype,controlname="CONTROL",level="transcript",sig=0.1))

# split by genotype and age effects
isoResTranEffects <- do.call(rbind,dissect_DESeq2(wald=isoResTranAnno$wald$anno_res,lrt=isoResTranAnno$lrt$anno_res))
isoResTranEffects


## ---------- Output -----------------

saveRDS(ontResTranAnno, file = paste0(dirnames$output, "/Ont_DESeq2TranscriptLevel.RDS"))
saveRDS(isoResTranAnno, file = paste0(dirnames$output, "/IsoSeq_DESeq2TranscriptLevel.RDS"))


## ---------- Iso-Seq Differential gene expression -----------------

# run DESeq2
resGene <- list(
  wald = run_DESeq2(test="Wald",input$gene_expression,input$isoPhenotype,exprowname="associated_gene",threshold=10,controlname="CONTROL",design="time_series",
                    interaction="On"),
  lrt = run_DESeq2(test="LRT",input$gene_expression,input$ontPhenotype,exprowname="associated_gene",threshold=10,controlname="CONTROL",design="time_series",
                   interaction="On")
)

# annotate results
resGeneAnno <- lapply(resGene, function(x) anno_DESeq2(x,input$classfiles,input$phenotype,controlname="CONTROL",level="gene",sig=0.1))