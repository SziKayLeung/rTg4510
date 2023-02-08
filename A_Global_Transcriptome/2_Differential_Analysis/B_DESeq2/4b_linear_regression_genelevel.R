## ---------- Script -----------------
##
## Purpose: perform differential analysis on mouse rTg4510 ONT targeted datasets using linear regression
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
library(pheatmap)
library(RColorBrewer)

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/transcriptome_stats/read_sq_classification.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/aesthetics_basics_plots/draw_density.R")

## ---------- input -----------------

dirnames <- list(
  rTg4510 = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome/",
  whole = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/2_post_isoseq3/",
  rnaseq = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/1_RNASeq_Isabel/"
)

input_files <- list(
  phenotype = paste0(dirnames$rTg4510, "0_metadata/A_isoseq_whole/WholeIsoSeq_PhenotypeTAPPAS.txt"), 
  expression = paste0(dirnames$whole, "7_tofu/WholeIsoSeq.Demultiplexed_Abundance.txt"),
  classfiles = paste0(dirnames$whole, "9_sqanti3/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt"),
  rnaseq_genotype = paste0(dirnames$rnaseq, "Isabel_Supp4_Tg4510AgeGenotypeDEG.csv")
)

input <- list()
input$phenotype <- read.table(input_files$phenotype, sep = "\t", header = T)
input$expression <- read.csv(input_files$expression)
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns") %>% mutate(PB_associated_gene = word(isoform,c(2), sep = fixed(".")))
input$rnaseq_genotype <- read.csv(input_files$rnaseq_genotype)


input$expression <- input$expression %>% mutate(associated_gene = word(id,c(2), sep = fixed(".")))
input$gene_expression <- aggregate(. ~ associated_gene, input$expression %>% select(-"id"), sum)
rownames(input$gene_expression) <- input$gene_expression$associated_gene


# input phenotpye sample ID has to match the sample ID in expression matrix 
# samples include "ont_" in expression matrix
input$phenotype$sample <- str_remove(as.character(input$phenotype$sample),"FL.")


## ---------- Creating DESeq2 object and analysis -----------------

de <- run_DESeq2(input$gene_expression, input$phenotype)

# annotate results 
de$res_Wald <- merge(de$res_Wald, unique(input$classfiles[,c("associated_gene","PB_associated_gene")]), by.x = 0, by.y = "PB_associated_gene") %>% 
  arrange(padj)

effectsize <- merge(de$res_Wald[,c("log2FoldChange","associated_gene")], input$rnaseq_genotype[,c("GenotypeTG.Age_months8","Gene")], 
                    by.x = "associated_gene", by.y = "Gene") 

density_plot(effectsize, "log2FoldChange", "GenotypeTG.Age_months8", "Iso-Seq log2FoldChange", "RNA-Seq log2FoldChange", "")


