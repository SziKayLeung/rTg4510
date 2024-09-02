## ---------- Script -----------------
##
## Purpose: perform differential analysis on mouse rTg4510 ONT and Iso-Seq targeted datasets using linear regression
## Transcript level separate analysis
## Gene level separate analysis
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

LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))


## ---------- input -----------------

# directory names
dirnames <- list(
  rTg4510 = "/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/",
  output = "/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/01_figures_tables/Targeted_Transcriptome"
)

# read input files
input_files <- list(
  phenotype = paste0(dirnames$rTg4510, "0_metadata/DESeq2SCNPhenotype.txt"), 
  
  # merged data
  expression = paste0(dirnames$rTg4510, "H_Sorted_Nuclei/5_cupcake/6_collapse/demux_fl_count.csv"),
  classfiles = paste0(dirnames$rTg4510, "H_Sorted_Nuclei/5_cupcake/7_sqanti3/rTg4510SCN_collapsed_RulesFilter_result_classification_counts.txt")
)


input <- list()
input$phenotype <- read.table(input_files$phenotype, sep = "\t", header = T) %>% mutate(SampleID = word(sample,c(1),sep=fixed("_")))
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")
input$expression <- data.table::fread(input_files$expression) %>% .[.$id != "0",] %>% filter(id %in% input$classfiles$isoform)


## ---------- ONT: Creating DESeq2 object and analysis -----------------
# remove outliers from phenotype
input$phenotype <- input$phenotype %>% filter(!SampleID %in% c("NeuN72", "NeuN65", "DN3"))

input$expression <- input$expression %>% tibble::column_to_rownames("id") 
input$expression <- input$expression %>% dplyr::select(input$phenotype$sample)
rownames(input$phenotype) <- input$phenotype$sample 
input$phenotype <- input$phenotype %>% dplyr::select(-sample)

if(all(colnames(input$expression) == rownames(input$phenotype))==FALSE){
  print("ERROR: rownames and colnames in expression matrix and phenotype are not in the same order")
  #input$expression %>% select(rownames(input$phenotype))
}

#input$expression$sum <- rowSums(input$expression[,1:15]) 
dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(input$expression)),
                              colData = input$phenotype,
                              design = ~ group + cell)
dds <- estimateSizeFactors(dds)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds_output <- DESeq(dds, test="Wald")
res <- as.data.frame(results(dds_output)) %>% tibble::rownames_to_column("isoform") %>% arrange(padj)
stats <- as.data.frame(mcols(dds_output))
output <- list(dds_output, res, norm, stats)
names(output) <- c("dds_Wald", "res_Wald","norm_counts", "stats_Wald")

#### 


MergedOutputStats <- output$stats_Wald %>% tibble::rownames_to_column(., var = "sorted_isoform") %>% 
  dplyr::select(sorted_isoform, cell_NeuN_vs_DN, WaldPvalue_cell_NeuN_vs_DN, WaldPvalue_group_WT_vs_TG) %>% 
  filter(sorted_isoform %in% mergeAll[mergeAll$dataset %in% c("Both"),"sorted_isoform"])

MergedOutputStats <- merge(MergedOutputStats, mergeAll, all.x = T)
write.csv(MergedOutputStats, "MergedOutputStats.csv")
