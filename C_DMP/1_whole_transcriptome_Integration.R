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
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/sqanti_general.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/plot_aesthetics.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/plot_tappas_analysis.R")

# Config files 
source(paste0(SC_ROOT,"A_Global_Transcriptome/2_Differential_Analysis/rTg4510_differential.config.R"))
source(paste0(SC_ROOT,"C_DMP/rTg4510_whole_methylation.config.R"))

## ---------- Load tappAS files -----------------

loaded <- list(
  iso = input_tappasfiles(TAPPAS_INPUT_DIR$iso),
  rna = input_tappasfiles(TAPPAS_INPUT_DIR$rna)
)


## ---------- Annotate tappAS files -----------------

annotated <- list(
  iso = annotate_tappasfiles(class.files,loaded$iso$input_normalized_matrix,phenotype$iso),
  rna = annotate_tappasfiles(class.files,loaded$rna$input_normalized_matrix,phenotype$rna)
)


## ---------- Identify differentially expressed and methylated genes -----------------

DMPGenotype_Integration <- Methylation_Integration("DMP_Genotype")
DMPInteraction_Integration <- Methylation_Integration("DMP_Interaction")
DMRGenotype_Integration <- Methylation_Integration("DMR_Genotype")

DMP_Combine <- rbind(DMPGenotype_Integration$Positions,DMPInteraction_Integration$Positions) 
DMP_Combine <- DMP_Combine[!duplicated(DMP_Combine), ]
length(unique(DMP_Combine$SYMBOL))


## ---------- Apply functions for plots ----------------

# Methyylation and Expression plots 
norm_in = loaded$rna$input_normalized_matrix
pheno_in = phenotype$rna
Spata13_output <- Methylation_Integration_plots("Spata13","PB.4966.2",norm_in,pheno_in)
Osmr_output <- Methylation_Integration_plots("Osmr","PB.5258.1",norm_in,pheno_in)
Jph1_output <- Methylation_Integration_plots("Jph1","NA",norm_in,pheno_in)
Ncf2_output <- Methylation_Integration_plots("Ncf2","PB.700.1",norm_in,pheno_in)
IRF8_output <- Methylation_Integration_plots("Irf8","PB.15969.1",norm_in,pheno_in)
Susd5_output <- Methylation_Integration_plots("Susd5","PB.16983.1",norm_in,pheno_in)
As3mt_methylation_output <- as3mt_dmr_figure()
Prnp_methylation_output <- prnp_dmr_figure()

# Plot Transcript Expression Trajection of DMP/DMR Genes
DMP_DMR_Genes <- unique(c(DMPGenotype_Integration$Positions$SYMBOL,
                         DMPInteraction_Integration$Positions$SYMBOL,
                         as.character(DMRGenotype_Integration$Positions$SYMBOL)))

DMP_DMR_Genes_plots_traj  <- lapply(lapply(DMP_DMR_Genes, function(gene) 
  plot_transexp_overtime_meth(gene,annotated$rna$Norm_transcounts,tappassigtrans$WholeIso_Transexp,"RNA-Seq Expression") +
    theme(legend.position = "bottom", plot.title = element_blank())),ggplotGrob)
names(DMP_DMR_Genes_plots_traj) <- DMP_DMR_Genes


## ---------- Output Plots ----------------

pdf(paste0(output_dir,"/WholeDifferentialAnalysis_DMPDMR.pdf"), width = 15, height = 8)
generate_cowplot(Prnp_methylation_output[[1]],Prnp_methylation_output[[2]],Prnp_methylation_output[[3]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$As3mt,As3mt_methylation_output[[1]],As3mt_methylation_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Spata13,Spata13_output[[1]],Spata13_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Ncf2,Ncf2_output[[1]],Ncf2_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Irf8,IRF8_output[[1]],IRF8_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Susd5,Susd5_output[[1]],Susd5_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Osmr,Osmr_output[[1]],Osmr_output[[2]],num=3,ncol=3,nrow=1)
dev.off()


## ---------- Output Stats -----------------
Meth_Table = lapply(unique(DMP_Combine$SYMBOL), 
                    function(x) Methylation_Integration_stats(x, tappassigtrans$WholeRNA_Transexp, annotated$rna))
Meth_Table <- do.call("rbind",Meth_Table) 
Meth_Table %>% mutate(loc_simple = word(location,c(1),sep = fixed(" "))) %>% group_by(loc_simple) %>% tally()
binom.test(14, 18, 0.5)
binom.test(5, 18, 0.5)

