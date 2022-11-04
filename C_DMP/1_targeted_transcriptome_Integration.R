## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
## 
##   
##
##


## ---------- Source function and config files ---------------
# Specific analysis related
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/"
source(paste0(SC_ROOT,"C_DMP/01_source_functions.R"))

# General functions
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/sqanti_general.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/plot_aesthetics.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/plot_tappas_analysis.R")

# Config files 
source(paste0(SC_ROOT,"B_Targeted_Transcriptome/3_Differential_Analysis/rTg4510_differential.config.R"))
source(paste0(SC_ROOT,"C_DMP/rTg4510_methylation.config.R"))


## ---------- Load tappAS files -----------------

loaded <- list(
  iso = input_tappasfiles(TAPPAS_INPUT_DIR$iso),
  ont = input_tappasfiles(TAPPAS_INPUT_DIR$ont)
)


## ---------- Annotate tappAS files -----------------

annotated <- list(
  iso = annotate_tappasfiles(class.files$iso,loaded$iso$input_normalized_matrix,phenotype$iso),
  ont = annotate_tappasfiles(class.files$ont,loaded$ont$input_normalized_matrix,phenotype$ont)
)


## ---------- Identify differentially expressed and methylated genes -----------------

# subset target genes in methylation output 
Diff_genes = list()
categories = c("genotype","interaction","pathology")
for(count in 1:3){
  i = categories[count]
  df <- Whole_DMP[[i]][Whole_DMP[[i]]$SYMBOL %in% TargetGene,]
  cat(i, "DMPs associated to target genes:", as.character(unique(df$SYMBOL)),"\n")
  Diff_genes[[count]] = df
}
names(Diff_genes) = categories

# generate list of DMP 
Diff_genes$pathology$X
dmp_pathology <- data.frame(
  symbol = Diff_genes$pathology$SYMBOL,
  chr = word(Diff_genes$pathology$X,c(1),sep = fixed(":")),
  start = word(Diff_genes$pathology$X,c(2),sep = fixed(":")),
  end = word(Diff_genes$pathology$X,c(2),sep = fixed(":"))
)
write.table(dmp_pathology, paste0(output_dir,"/dmp_pathology.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)

## ---------- Apply functions for plots ----------------

# Methyylation and Expression plots 
norm_in = loaded$ont$input_normalized_matrix
pheno_in = phenotype$ont
Bin1_output <- Methylation_Integration_plots("Bin1","ENSMUST00000234496.1",norm_in,pheno_in)
Clu_output <- Methylation_Integration_plots("Clu","TALONT000440029",norm_in,pheno_in)
Snca_output <- Methylation_Integration_plots("Clu","TALONT001103657",norm_in,pheno_in)


# Transcript Expression Trajection of DMP/DMR Genes
DMP_DMR_Genes = c("Bin1","Snca","Clu")
DMP_DMR_Genes_plots_traj  <- lapply(lapply(DMP_DMR_Genes, function(gene) 
  plot_transexp_overtime_meth(gene,annotated$ont$Norm_transcounts,tappassigtrans$ont$TargetedOnt_Transexp,"ONT Expression") + 
    theme(legend.position = "bottom", plot.title = element_blank())),ggplotGrob)
names(DMP_DMR_Genes_plots_traj) <- DMP_DMR_Genes



## ---------- Output Plots ----------------

pdf(paste0(output_dir,"/TargetedDifferentialAnalysis_DMPDMR.pdf"), width = 15, height = 8)
generate_cowplot(DMP_DMR_Genes_plots_traj$Bin1,Bin1_output[[1]],Bin1_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Clu,Clu_output[[1]],Clu_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Snca,Snca_output[[1]],Snca_output[[2]],num=3,ncol=3,nrow=1)
dev.off()

