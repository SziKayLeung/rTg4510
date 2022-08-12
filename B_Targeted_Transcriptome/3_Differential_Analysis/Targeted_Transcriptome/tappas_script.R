### Differential Analysis 
tappasiso <- input_tappasfiles(tappasiso_input_dir)
tappasrna <- input_tappasfiles(tappasrna_input_dir)
tappas_removediso(tappasiso$tappAS_Transcripts_InputExpressionMatrix.tsv)
tappas_removediso(tappasrna$tappAS_Transcripts_InputExpressionMatrix.tsv)
targetedtappas_isoexp <- tappas_resultsanno(targeted.class.files,tappasiso$input_normalized_matrix.tsv,tappas_phenotype)
targetedtappas_rnaexp <- tappas_resultsanno(targeted.class.files,tappasrna$input_normalized_matrix.tsv,tappas_phenotype)
tappasgenesig_plots <- tappas_genesig()

# Gene Expression Plots 
# using IsoSeq FL read count as expression input
targeted_isogeneexp_plots <- lapply(lapply(unique(targetedtappas_isoexp$GeneExp$associated_gene), function(gene) plot_mergedexp(gene,"NA",targetedtappas_isoexp$GeneExp,targetedtappas_isoexp$Norm_transcounts)),ggplotGrob)
names(targeted_isogeneexp_plots) <- unique(targetedtappas_isoexp$GeneExp$associated_gene)
# using RNASeq abundance as expression input
targeted_rnageneexp_plots <- lapply(lapply(unique(targetedtappas_rnaexp$GeneExp$associated_gene), function(gene) plot_mergedexp(gene,"NA", targetedtappas_rnaexp$GeneExp, targetedtappas_rnaexp$Norm_transcounts)),ggplotGrob)
names(targeted_rnageneexp_plots) <- unique(targetedtappas_rnaexp$GeneExp$associated_gene)

pdf (paste0(output_plot_dir,"/TargetedDifferentialAnalysis.pdf"), width = 10, height = 15)
# IsoSeq FL as expression input 
group_plots(ADReg_Genes,targeted_isogeneexp_plots)
group_plots(GWAS_Genes,targeted_isogeneexp_plots)
group_plots(FTD_Genes,targeted_isogeneexp_plots)
group_plots(EWAS_Genes,targeted_isogeneexp_plots)
# RNASeq as expression input 
group_plots(ADReg_Genes,targeted_rnageneexp_plots)
group_plots(GWAS_Genes,targeted_rnageneexp_plots)
group_plots(FTD_Genes,targeted_rnageneexp_plots)
group_plots(EWAS_Genes,targeted_rnageneexp_plots)
# Trem2 
plot_grid(targeted_isogeneexp_plots$Trem2,plot_transexp_overtime("Trem2",targetedtappas_isoexp$Norm_transcounts), plot_transexp("Trem2",targetedtappas_isoexp$Norm_transcounts,"isoseq_targeted"),labels = c("a","b","c"), nrow = 3, label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
dev.off()

pdf (paste0(output_plot_dir,"/TargetedDifferentialAnalysis_RNAvsIso.pdf"), width = 10, height = 7)
tappasgenesig_plots 
for(gene in c(ADReg_Genes,GWAS_Genes,FTD_Genes,EWAS_Genes)){print(plot_grid(group_plots_rnavsiso(gene,targeted_isogeneexp_plots,targeted_rnageneexp_plots)))}
dev.off()
