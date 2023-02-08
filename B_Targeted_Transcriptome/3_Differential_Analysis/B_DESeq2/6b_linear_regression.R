## ---------- Script -----------------
##
## Purpose: perform differential analysis on mouse rTg4510 ONT targeted datasets using linear regression
## using ONT targeted dataset collapsed from tofu-cupcake
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
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/results/cupcake_collapse/"
)

# read input files
input_files <- list(
  phenotype = paste0(dirnames$rTg4510, "0_metadata/F_ont_targeted/ONT_phenotype.txt"), 
  expression = paste0(dirnames$rTg4510, "F_ONT_Targeted/6b_tofu_sqanti3/demux_fl_count.txt"),
  classfiles = paste0(dirnames$rTg4510, "F_ONT_Targeted/6b_tofu_sqanti3/all_merged_RulesFilter_result_classification.txt")
)

input <- list()
input$phenotype <- read.table(input_files$phenotype, sep = "\t", header = T)
input$expression <- read.csv(input_files$expression) %>% .[.$isoform != "0",]
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")


## ---------- Creating DESeq2 object and analysis -----------------

# results from running DESeq2
res <- list()
res$wald <- run_DESeq2(input$expression,input$phenotype,exprowname="isoform",controlname="CONTROL",interaction="On",test="Wald")
res$lrt <- run_DESeq2(input$expression,input$phenotype,exprowname="isoform",controlname="CONTROL",interaction="On",test="LRT")

# annotate normalised counts
anno_res <- lapply(res, function(x) anno_DESeq2(x, input$classfiles, input$phenotype, controlname="CONTROL"))


# function
time_case_boxplot <- function(normalised_counts, transcript){
  
  df <- normalised_counts %>% filter(isoform == transcript)
  
  p <- ggplot(df, aes(x = group, y = normalised_counts)) + geom_boxplot() +
    geom_point(position="jitter",aes(color = as.factor(time)), size = 3) +
    theme_bw() +
    labs(x = "Genotype", y = "Normalised counts", 
         title = paste0(unique(df$associated_gene),": ", transcript),
         subtitle = df$associated_transcript) +
    scale_colour_manual(name = "Age (months)",
                        values = c(wes_palette("Darjeeling2")[[5]], wes_palette("Zissou1")[[1]], 
                                   wes_palette("Zissou1")[[3]],wes_palette("Zissou1")[[5]]))
  
  return(p)
}




# below 0.1 threshold (p-value)
sig0.1 <- res$wald$res_Wald %>% filter(padj < 0.1) 
sig0.1plots <- lapply(sig0.1$isoform, function(x) time_case_boxplot(anno_res$wald$norm_counts, x))
names(sig0.1plots) <- sig0.1$isoform

# apply and output
pdf(paste0(dirnames$output,"Sig0.1_ONT_DESeq2.pdf"), width = 8, height = 5)
for(i in sig0.1plots){
  print(i)
}
dev.off()
write.csv(res_df, paste0(dirnames$output, "ONT_DESeq2_results.csv"))


# bin1
time_case_boxplot("PB.3915.33_ENSMUST00000234496.1")
