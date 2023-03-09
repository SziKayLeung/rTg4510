## ---------- Script -----------------
##
## Purpose: perform differential analysis on mouse rTg4510 ONT targeted datasets using linear regression
## Transcript level - Iso-Seq, RNA-Seq (hybrid, already aligned to Iso-Seq)
## Gene level - Iso-Seq, RNA-Seq
## 1/ run DESeq2 (Wald and LRT)
## 2/ annotate results using class.files
## 3/ split by genotype and age effects
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
SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome/1_IsoSeq_Pipeline/"
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_transcript_level.R"))
source(paste0(SC_ROOT, "02_source_characterise_functions.R"))


## ---------- input -----------------

dirnames <- list(
  rTg4510 = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/Whole_Transcriptome",
  whole = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/",
  rnaseq = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/1_RNASeq_Isabel/"
)

input_files <- list(
  phenotype = paste0(dirnames$rTg4510, "0_metadata/A_isoseq_whole/WholeIsoSeq_PhenotypeTAPPAS.txt"),
  expression = paste0(dirnames$whole, "2_post_isoseq3/7_tofu/WholeIsoSeq_fl_count.csv"),
  classfiles = paste0(dirnames$whole, "2_post_isoseq3/9_sqanti3/WholeIsoSeq.collapsed_RulesFilter_result_classification.txt"),
  rnaseq_expression = paste0(dirnames$whole, "2_post_isoseq3/10_rnaseq2isoseq/WholeIsoSeq_rnaseq.expression.txt"),
  rnaseq_phenotype = paste0(dirnames$whole, "3_differential/1_Input/B_hybrid_59/WholeAllMouse_PhenotypeTAPPAS.txt"),
  rnaseq_genotype = paste0(dirnames$rnaseq, "Isabel_Supp2_Tg4510GenotypeDEG.csv"),
  rnaseq_progressive = paste0(dirnames$rnaseq, "Isabel_Supp4_Tg4510AgeGenotypeDEG.csv")
)

input <- list()
# phenotype
input$phenotype <- read.table(input_files$phenotype, sep = "\t", header = T)
input$phenotype$sample <- word(str_remove(as.character(input$phenotype$sample),"FL."),c(1),sep=fixed("_"))
# expression
input$expression <- read.csv(input_files$expression)
input$gene_expression <- input$expression %>% mutate(associated_gene = word(id,c(2), sep = fixed(".")))
input$gene_expression <- aggregate(. ~ associated_gene, input$gene_expression %>% select(-"id"), sum)
rownames(input$gene_expression) <- input$gene_expression$associated_gene
# SQANTI class files
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")
# RNA-Seq differential results
input$rnaseq_genotype <- read.csv(input_files$rnaseq_genotype)
input$rnaseq_progressive <- read.csv(input_files$rnaseq_progressive)
input$rnaseq_phenotype <- read.table(input_files$rnaseq_phenotype, sep = "\t", header = T)
input$rnaseq_expression <- read.table(input_files$rnaseq_expression, sep = "\t", header = T)
input$rnaseq_gene_expression <- input$rnaseq_expression %>% mutate(associated_gene = word(X,c(2), sep = fixed(".")))
input$rnaseq_gene_expression <- aggregate(. ~ associated_gene, input$rnaseq_gene_expression %>% select(-"X"), sum)
rownames(input$rnaseq_gene_expression) <- input$rnaseq_gene_expression$associated_gene


## ---------- Iso-Seq Differential gene expression -----------------

# run DESeq2
resGene <- list(
  wald = run_DESeq2(test="Wald",input$gene_expression,input$phenotype,exprowname="associated_gene",threshold=10,controlname="CONTROL",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$gene_expression,input$phenotype,exprowname="associated_gene",threshold=10,controlname="CONTROL",interaction="On")
)

# annotate results
resGeneAnno <- lapply(resGene, function(x) anno_DESeq2(x,input$classfiles,input$phenotype,controlname="CONTROL",level="gene",sig=0.1))


## ---------- Iso-Seq Differential transcript expression -----------------

# run DESeq2
resTran <- list(
  wald = run_DESeq2(test="Wald",input$expression,input$phenotype,threshold=10,exprowname="id",controlname="CONTROL",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$expression,input$phenotype,threshold=10,exprowname="id",controlname="CONTROL",interaction="On")
)

# annotate results and filter at 0.05
resTranAnno <- lapply(resTran, function(x) anno_DESeq2(x,input$classfiles,input$phenotype,controlname="CONTROL",level="transcript",sig=0.1))

# split by genotype and age effects
resTranEffects <- dissect_DESeq2(wald=resTranAnno$wald$anno_res,lrt=resTranAnno$lrt$anno_res)
resTranEffects <- do.call(rbind,resTranEffects)
resTranEffects


## ---------- RNA-Seq Differential transcript expression -----------------

# results from running DESeq2
RresTran <- list(
  wald = run_DESeq2(test="Wald",input$rnaseq_expression,input$rnaseq_phenotype,threshold=10,exprowname="X",controlname="CONTROL",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$rnaseq_expression,input$rnaseq_phenotype,threshold=10,exprowname="X",controlname="CONTROL",interaction="On")
)

# annotate results
RresTranAnno <- lapply(RresTran, function(x) anno_DESeq2(x,input$classfiles,input$rnaseq_phenotype,controlname="CONTROL",level="transcript"))

# split by genotype and age effects
RresTranEffects <- dissect_DESeq2(RresTranAnno$wald$anno_res,RresTranAnno$lrt$anno_res)


## ---------- RNA-Seq Differential gene expression -----------------

# results from running DESeq2
RresGene <- list(
  wald = run_DESeq2(test="Wald",input$rnaseq_gene_expression,input$rnaseq_phenotype,threshold=10,exprowname="associated_gene",controlname="CONTROL",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$rnaseq_gene_expression,input$rnaseq_phenotype,threshold=10,exprowname="associated_gene",controlname="CONTROL",interaction="On")
)

# annotate results
RresGeneAnno <- lapply(RresGene, function(x) anno_DESeq2(x, input$classfiles, input$rnaseq_phenotype, controlname="CONTROL", level="gene", sig=0.1))


## ---------- Gene level comparison -----------------

# comparison of effect size (statistics) at the gene level between Iso-Seq and RNA-Seq (Isabel's results)
resGeneComparison <- list(
  # genotype level, like-to-like comparison with wald test 
  genotype = merge(resGeneAnno$wald$stats_Wald[,c("associated_gene","WaldStatistic_group_CASE_vs_CONTROL")],
                   input$rnaseq_genotype[,c("Gene", "WaldStatistic_Genotype_TG_vs_WT")],
                   by.x = "associated_gene", by.y = "Gene"),
  
  # interaction level, like-to-like copmarison with lrt test
  interaction = merge(resGeneAnno$lrt$stats_LRT[,c("associated_gene","LRTStatistic")],
                      input$rnaseq_progressive[,c("Gene", "LRTStatistic")],
                      by.x = "associated_gene", by.y = "Gene")
)

## ---------- Output -----------------

saveRDS(resTranAnno, file = paste0(dirnames$output, "/IsoSeq_DESeq2TranscriptLevel.RDS"))
saveRDS(resGeneAnno, file = paste0(dirnames$output, "/IsoSeq_DESeq2GeneLevel.RDS"))
saveRDS(RresTranAnno, file = paste0(dirnames$output, "/RNASeqHybrid_DESeq2TranscriptLevel.RDS"))
saveRDS(RresGeneAnno, file = paste0(dirnames$output, "/RNASeqHybrid_DESeq2GeneLevel.RDS"))
saveRDS(resGeneComparison, file = paste0(dirnames$output, "/Comparison_DESeq2GeneLevel.RDS"))
