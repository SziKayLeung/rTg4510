# Szi Kay Leung
# Call Functions script for Thesis Chapter on Whole Transcriptome IsoSeq

source_root_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Figures_Thesis/"
output_helpfig_dir = paste0(source_root_dir, "Tables4Figures")
output_plot_dir = paste0(source_root_dir, "DiffAnalysis_WholeTranscriptome/Pdf")

# results from Whole transcriptome paper
input_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/Whole_Transcriptome_Paper/Output/Tables"

# source function scripts
source(paste0(source_root_dir,"DiffAnalysis_WholeTranscriptome/IsoSeq_WholeDiffAnalysis_Functions.R"))    ## Global Functions
source(paste0(source_root_dir, "DiffAnalysis_WholeTranscriptome/IsoSeq_WholeDiffAnalysis_Variables.R"))   ## Variables
source(paste0(source_root_dir,"DiffAnalysis_WholeTranscriptome/AS_Functions.R"))                          ## AS Functions

##### Alternative Splicing####################
dASevents_files = read_files_differentialASevents(AS_dir)
group.mono.class.files = dASevents_files$group.mono.class.files
annotated_sqanti_gtf = dASevents_files$annotated_sqanti_gtf
suppa2.output= dASevents_files$suppa2.output
AS_output = AS_events_diff(group.mono.class.files,dASevents_files$annotated_sqanti_gtf,dASevents_files$suppa2.output)

##### Load Files ############################# 
tappasiso <- input_tappasfiles(tappas_dir$glob_iso, "isoseq")
tappasrna <- input_tappasfiles(tappas_dir$glob_rna, "rnaseq")

# Annotate files
IsoExp <- tappas_resultsanno(class.files$glob_iso,tappasiso$input_norm,phenotype$glob_iso)
RNAExp <- tappas_resultsanno(class.files$glob_iso,tappasrna$input_norm,phenotype$glob_rna)

#####  Different models from Tappas output ############################# 
# segregate differential gene and transcript expression results by the different models (using beta coefficient as filters)
# tappassiggene and tappassigtrans input differential expression results
# IsoSeq as Expression Input
gene_sigs_WholeIso_lst = segregate_tappasresults(tappassiggene$WholeIso_Genexp,"IsoSeq")
trans_sigs_WholeIso_lst = segregate_tappasresults(tappassigtrans$WholeIso_Transexp,"IsoSeq")

# RNASeq as Expression Input
gene_sigs_WholeRNA_lst = segregate_tappasresults(tappassiggene$WholeRNA_Genexp,"RNASeq")
trans_sig_WholeRNA_lst = segregate_tappasresults(tappassigtrans$WholeRNA_Transexp,"RNASeq")

# Find examples of the different models using the 3rd gene from each list 
# 3rd gene rather than top 2 for more variation in thesis
genesigs_model = sapply(gene_sigs_WholeIso_lst$models, function(x) x[3,1]) %>% reshape2::melt(value.name = "Gene")
genesigs_model_plots = lapply(lapply(1:length(genesigs_model$Gene), function(i) plot_mergedexp(genesigs_model$Gene[[i]],"NA",IsoExp$GeneExp,IsoExp$Norm_transcounts,paste0("Model ",i))),ggplotGrob)

gene_sigs_seg = rbind(
  gene_sigs_WholeIso_lst$models$`Model 1 Genotype` %>% mutate(Effect = "Genotype"),
  gene_sigs_WholeIso_lst$models$`Model 2 Genotype+Age` %>% mutate(Effect = "Genotype & Age"),    
  gene_sigs_WholeIso_lst$models$`Model 3 Age` %>% mutate(Effect = "Age"),
  gene_sigs_WholeIso_lst$models$`Model 4 Interaction` %>% mutate(Effect = "Interaction"),
  gene_sigs_WholeIso_lst$models$`Model 5 Interaction` %>% mutate(Effect = "Interaction"),
  gene_sigs_WholeIso_lst$models$`Model 6 Interaction` %>% mutate(Effect = "Interaction"),
  gene_sigs_WholeIso_lst$models$`Model 7 Interaction` %>% mutate(Effect = "Interaction"))
write.table(gene_sigs_seg, paste0(output_helpfig_dir,"/DGE_Effect_WholeIsoSeq.csv"), sep = ",", row.names=FALSE)


trans_sigs_seg = rbind(
      trans_sigs_WholeIso_lst$models$`Model 1 Genotype` %>% mutate(Effect = "Genotype"),
      trans_sigs_WholeIso_lst$models$`Model 2 Genotype+Age` %>% mutate(Effect = "Genotype & Age"),    
      trans_sigs_WholeIso_lst$models$`Model 3 Age` %>% mutate(Effect = "Age"),
      trans_sigs_WholeIso_lst$models$`Model 4 Interaction` %>% mutate(Effect = "Interaction"),
      trans_sigs_WholeIso_lst$models$`Model 5 Interaction` %>% mutate(Effect = "Interaction"),
      trans_sigs_WholeIso_lst$models$`Model 6 Interaction` %>% mutate(Effect = "Interaction"),
      trans_sigs_WholeIso_lst$models$`Model 7 Interaction` %>% mutate(Effect = "Interaction"))
write.table(trans_sigs_seg, paste0(output_helpfig_dir,"/DTE_Effect_WholeIsoSeq.csv"), sep = ",", row.names=FALSE)

Genesigs = rbind(gene_sigs_WholeIso_lst$models$`Model 1 Genotype`,
                 gene_sigs_WholeIso_lst$models$`Model 2 Genotype+Age`,
                 gene_sigs_WholeIso_lst$models$`Model 4 Interaction`,
                 gene_sigs_WholeIso_lst$models$`Model 5 Interaction`,
                 gene_sigs_WholeIso_lst$models$`Model 6 Interaction`,
                 gene_sigs_WholeIso_lst$models$`Model 7 Interaction`)

RNAGenesigs = rbind(gene_sigs_WholeRNA_lst$models$`Model 1 Genotype`,
                    gene_sigs_WholeRNA_lst$models$`Model 2 Genotype+Age`,
                    gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`)
##### Differential Gene Expression ############################# 

### List of Interesting Genes 
### 1. Significant genes 
## siggenes = top3 from IsoSeq+IsoSeq except for Mapt = c("Gfap","C4b","Tgfbr1","Slc14a1","Unc93b1","Mapt")
### 2. Novel genes
# novelgenes = c("novelGene_529","novelGene_Fgfr1op_AS","Fgfr1op","novelGene_Htra1_AS","Htra1")
#gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`[grepl("novelGene",gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1),]
## common novel genes identified with significant gene expression changes, using IsoSeq and RNA-Seq as expression
#signovelGenes = intersect(c(tappassiggene$WholeIso_Genexp[grepl("novelGene",tappassiggene$WholeIso_Genexp$...1),1])$`...1`,
#                          c(tappassiggene$WholeRNA_Genexp[grepl("novelGene",tappassiggene$WholeRNA_Genexp$...1),1])$`...1`)
#tappassiggene$WholeRNA_Genexp[grepl("novelGene",tappassiggene$WholeRNA_Genexp$...1),] %>% filter(`R-squared` > 0.5)
### 3. Target genes
# Not all target genes detected in Whole transcriptome, or filtered out due to low count
#TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")
#TargetGene_detected = list() 
#for(gene in TargetGene){if(gene %in% IsoExp$Norm_transcounts$associated_gene ){TargetGene_detected[[gene]] = gene}}
#TargetGene_detected = row.names(do.call(rbind, TargetGene_detected))
siggenes = c("Gfap","C4b","Tgfbr1","Slc14a1","Unc93b1","Mapt", "novelGene_529","novelGene_Fgfr1op_AS","Fgfr1op","novelGene_Htra1_AS","Htra1",
             "Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Fus","Picalm","Snca","Apoe","Trem2","Vgf","Esyt2")

pIsoGeneExp <- generate_diff_plots(siggenes,IsoExp,"Gene","Iso-Seq Expression") # using IsoSeq FL read count as expression input
pRnaGeneExp <- generate_diff_plots(siggenes,RNAExp,"Gene","RNA-Seq Expression") # using RNASeq abundance as expression input

meanexp_output = list()
meanexp_output = lapply(c("Gfap","C4b","Tgfbr1","Slc14a1","Pros1","Unc93b1"), function(x) tabulate_diffgene(x, IsoExp$GeneExp, Interaction))
names(meanexp_output) = c("Gfap","C4b","Tgfbr1","Slc14a1","PRos1","Unc93b1")
meanexp_output <- do.call("rbind",meanexp_output) 

# Concurrence with Isabel's results (IsoSeq significant gene expression results) 
cat("Common significant genes identified with Isabel's results:", length(intersect(Genesigs$...1,Isabel_gene_Tg4510GenotypeDEG$Gene)),"\n")
cat("Common significant genes identified with Isabel's results:", length(intersect(Genesigs$...1,Isabel_gene_Tg4510AgeGenotypeDEG$Gene)),"\n")

# Concurrence with Isabel's results (RNASeq significant gene expression results) 
cat("Common significant genes identified with Isabel's results:", length(intersect(RNAGenesigs$...1,Isabel_gene_Tg4510GenotypeDEG$Gene)),"\n")
cat("Common significant genes identified with Isabel's results:", length(intersect(RNAGenesigs$...1,Isabel_gene_Tg4510AgeGenotypeDEG$Gene)),"\n")

nrow(Isabel_gene_Tg4510GenotypeDEG)
nrow(Isabel_gene_Tg4510AgeGenotypeDEG)
##### Differential Transcript Expression ############################# 
View(rbind(gene_sigs_WholeIso_lst$models$`Model 4 Interaction`,
           gene_sigs_WholeIso_lst$models$`Model 5 Interaction`,
           gene_sigs_WholeIso_lst$models$`Model 6 Interaction`,
           gene_sigs_WholeIso_lst$models$`Model 7 Interaction`)) 

# Heatmap 
pheat = lapply(c("Gfap","C4b","Ctsd","H2-D1","Gatm","Padi2"), function(x) draw_heatmap_gene(x,class.files,IsoExp$Norm_transcounts,"whole"))
names(pheat) = c("Gfap","C4b","Ctsd","H2-D1","Gatm","Padi2")

# Transcript Expression Plots 
sigtrans = c("Gfap","C4b", "Cd68","Osmr","Cd34","Ubqln1","Gjb2","Adam23","Mapt","Ctsd","H2-D1","Gatm","Padi2",
             "Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Fus","Picalm","Snca","Apoe","Trem2","Vgf",
             "Esyt2","Slc1a3","Gja1")

sigtrans = c("Ctsd","H2-D1","Gatm","Padi2","Cd34","Ubqln1")
# using IsoSeq FL read count as expression input
# mean expression per isoform plot  
pIsoTranExp <- generate_diff_plots(sigtrans,IsoExp,"Transcript","Iso-Seq Expression")
pRnaTranExp <- generate_diff_plots(sigtrans,RNAExp,"Transcript_Rnaseq_Targeted","RNA-Seq Expression")
pIsoTranExpTraj <- generate_diff_plots(sigtrans,IsoExp,"Transcript_Isoseq Trajectory","Iso-Seq Expression")
pRnaTranExpTraj <- generate_diff_plots(sigtrans,RNAExp,"Transcript_Rnaseq Trajectory","RNA-Seq Expression")

plot_transexp_overtime("Ubqln1", IsoExp,"Transcript","Iso-Seq Expression")
pdf(paste0(output_plot_dir,"/GfapC4b.pdf"), width = 14, height = 12)
plot_grid(generate_cowplot(pheat$Gfap$gtable,pIsoTranExp$Gfap,num="rwBC",ncol=2,nrow=1),
          generate_cowplot(pIsoTranExpTraj$Gfap,pRnaTranExpTraj$Gfap,num="2DE",ncol=2,nrow = 1), ncol = 1)

plot_grid(generate_cowplot(pheat$C4b$gtable,pIsoTranExp$C4b,num="rwBC",ncol=2,nrow=1),
          generate_cowplot(pIsoTranExpTraj$C4b,pRnaTranExpTraj$C4b,num="2DE",ncol=2,nrow = 1), ncol = 1)
dev.off()

pdf(paste0(output_plot_dir,"/DiffAD.pdf"), width = 19, height = 8)
generate_cowplot(pheat$Ctsd$gtable, pIsoTranExpTraj$Ctsd, pRnaTranExpTraj$Ctsd,num="BCD",ncol=3,nrow=1)
generate_cowplot(pheat$`H2-D1`$gtable, pIsoTranExpTraj$`H2-D1`, pRnaTranExpTraj$`H2-D1`,num="BCD",ncol=3,nrow=1)
generate_cowplot(pheat$Gatm$gtable, pIsoTranExpTraj$Gatm, pRnaTranExpTraj$Gatm,num="BCD",ncol=3,nrow=1)
generate_cowplot(pheat$Padi2$gtable, pIsoTranExpTraj$Padi2, pRnaTranExpTraj$Padi2,num="BCD",ncol=3,nrow=1)
dev.off()

Interaction_Difftrans = rbind(trans_sigs_WholeIso_lst$models$`Model 2 Genotype+Age`,
                              trans_sigs_WholeIso_lst$models$`Model 4 Interaction`,
                              trans_sigs_WholeIso_lst$models$`Model 5 Interaction`,
                              trans_sigs_WholeIso_lst$models$`Model 6 Interaction`,
                              trans_sigs_WholeIso_lst$models$`Model 7 Interaction`)
write.csv(Interaction_Difftrans, paste0(fig_dir,"/Interaction_WholeIsoSeq_DiffTrans.csv"))
All_FDR = bind_rows(trans_sigs_WholeIso_lst$models)

pdf(paste0(output_plot_dir,"/DiffAD2.pdf"), width = 14, height = 12)
generate_cowplot(pIsoTranExpTraj$Ubqln1,pRnaTranExpTraj$Ubqln1,pIsoTranExpTraj$Cd34,pRnaTranExpTraj$Cd34,num=4,ncol=2,nrow=2)
dev.off()

# differentially expressed isoforms
#tlist = c("PB.4255.13","PB.7004.8","PB.2972.16","PB.10959.1","PB.1036.2","PB.15108.6","PB.11607.2","PB.9298.1","PB.7039.1")
tlist = Interaction_Difftrans$isoform

meanexp_output = list()
meanexp_output = lapply(tlist, function(x) tabulate_difftrans(x, IsoExp$Norm_transcounts, All_FDR, "isoseq"))
names(meanexp_output) = tlist
meanexp_output <- do.call("rbind",meanexp_output) 
meanexp_output = meanexp_output %>% mutate(dir = ifelse(log2fc < 0 ,"Down","Up"))
table(meanexp_output$dir)
res = binom.test(510, 673, 0.5)
res$p.value

# Isoforms unique to iso-seq
DEI_genes_unique()


##### Differential Transcript Usage ############################# 
venndiu = DIU_analysis_output_venn()
IFiso_Esyt2 = IF_plot("Esyt2",tappasiso$gene_transcripts.tsv, tappasiso$input_norm, "isoseq")
IFrna_Esyt2 = IF_plot("Esyt2",tappasrna$gene_transcripts.tsv, tappasrna$input_norm, "rnaseq")
DIU_genes_exp_plots = DIU_genes_exp()
DIU_genes = DIU_RNASEQ_results()
write.xlsx(DIU_genes$DIU_DEG_maj, paste0(fig_dir,"/DIU_DEG.xlsx"), sheetName="DIU_DEG_maj")
write.xlsx(DIU_genes$DIU_DEG_nomaj, paste0(fig_dir,"/DIU_DEG.xlsx"), sheetName="DIU_DEG_nomaj",append=TRUE)
write.xlsx(DIU_genes$DIU_notDEG_nomaj, paste0(fig_dir,"/DIU_DEG.xlsx"), sheetName="DIU_notDEG_nomaj",append=TRUE)
write.xlsx(DIU_genes$DIU_notDEG_maj, paste0(fig_dir,"/DIU_DEG.xlsx"), sheetName="DIU_notDEG_maj",append=TRUE)

Cisd3 = diff_across_rnaseq("Cisd3")
Shisa = diff_across_rnaseq("Shisa5")
Fblim1 = diff_across_rnaseq("Fblim1")
Arpc4_Ttll3 = diff_across_rnaseq("Arpc4_Ttll3")

pdf(paste0(output_plot_dir,"/DIU.pdf"), width = 11, height = 12)
generate_cowplot(Cisd3[[1]],Cisd3[[3]],Cisd3[[2]],Cisd3[[4]],num=4,nrow=2,ncol=2)
generate_cowplot(Shisa[[1]],Shisa[[3]],Shisa[[2]],Shisa[[4]],num=4,nrow=2,ncol=2)
generate_cowplot(Fblim1[[1]],Fblim1[[3]],Fblim1[[2]],Fblim1[[4]],num=4,nrow=2,ncol=2)
generate_cowplot(Arpc4_Ttll3[[1]],Arpc4_Ttll3[[3]],Arpc4_Ttll3[[2]],Arpc4_Ttll3[[4]],num=4,nrow=2,ncol=2)
dev.off()
#################################### Generate Plots ############

pdf(paste0(output_plot_dir,"/WholeDifferentialAnalysis.pdf"), width = 10, height = 15)
generate_cowplot(pIsoGeneExp$Gfap,pRnaGeneExp$Gfap,pIsoGeneExp$C4b,pRnaGeneExp$C4b,pIsoGeneExp$Tgfbr1,pRnaGeneExp$Tgfbr1,num=6,ncol=2,nrow=3)
generate_cowplot(pIsoGeneExp$Slc14a1,pRnaGeneExp$Slc14a1,pIsoGeneExp$Unc93b1,pRnaGeneExp$Unc93b1,NULL,NULL,num=4,ncol=2,nrow=3)
generate_cowplot(pIsoTranExpTraj$Gfap,pIsoTranExp$Gfap,num=2,ncol=1,nrow=2)
generate_cowplot(gene_sigs_WholeIso_lst$p,venn_WGCNA(),NULL,NULL,num="rw",ncol=2,nrow=2)
generate_cowplot(plotlist = genesigs_model_plots,num=7,ncol=2,nrow=4)
generate_cowplot(grobTree(venndiu$v1),grobTree(venndiu$v2),grobTree(venndiu$v3),grobTree(venndiu$v4),num=4,ncol=1,nrow=4)
generate_cowplot(pIsoTranExp$Esyt2,pRnaTranExp$Esyt2,IFiso_Esyt2[[3]],IFrna_Esyt2[[3]],num=4,ncol=2,nrow=2)
generate_cowplot(DIU_genes_exp_plots[[2]],DIU_genes_exp_plots[[3]],NULL,NULL,IFrna_Esyt2[[3]],num=2,ncol=2,nrow=2)
plot_grid(generate_cowplot(pRnaGeneExp$novelGene_529,NULL,num="A",ncol=2,nrow=2),
          generate_cowplot(pRnaGeneExp$novelGene_Fgfr1op_AS,pRnaGeneExp$Fgfr1op,num="BC",ncol=2,nrow=2),
          generate_cowplot(pRnaGeneExp$novelGene_Htra1_AS,pRnaGeneExp$Htra1,num="DE",ncol=2,nrow=2), nrow=3, ncol = 1)
generate_cowplot(pIsoTranExpTraj$C4b,pIsoTranExp$C4b,num=2,ncol=1,nrow=2)
generate_cowplot(pRnaTranExpTraj$Gfap,pRnaTranExpTraj$C4b,num=2,ncol=1,nrow=2)
generate_cowplot(pIsoTranExpTraj$Osmr,pRnaTranExpTraj$Osmr,pIsoTranExpTraj$Cd68,pRnaTranExpTraj$Cd68,num=4,ncol=2,nrow=2)
generate_cowplot(pIsoTranExpTraj$Cd34,pRnaTranExpTraj$Cd34,pIsoTranExpTraj$Ubqln1,pRnaTranExpTraj$Ubqln1,num=4,ncol=2,nrow=2)
generate_cowplot(pIsoTranExpTraj$Cd34,pRnaTranExpTraj$Cd34,pIsoTranExpTraj$Ubqln1,pRnaTranExpTraj$Ubqln1,num=4,ncol=2,nrow=2)
generate_cowplot(pIsoTranExpTraj$Slc1a3,pRnaTranExpTraj$Slc1a3,pIsoTranExpTraj$Gja1,pRnaTranExpTraj$Gja1,num=4,ncol=2,nrow=2)
generate_cowplot(pIsoTranExpTraj$Ctsd,pIsoTranExpTraj$`H2-D1`,pIsoTranExpTraj$Gatm,pIsoTranExpTraj$Padi2,num=4,ncol=2,nrow=2)
generate_cowplot(pRnaTranExpTraj$Ctsd,pRnaTranExpTraj$`H2-D1`,pRnaTranExpTraj$Gatm,pRnaTranExpTraj$Padi2,num=4,ncol=2,nrow=2)
dev.off()

pdf(paste0(output_plot_dir,"/NovelGeneExp.pdf"), width = 10, height = 15)
plot_grid(generate_cowplot(pRnaGeneExp$novelGene_529,NULL,num="A",ncol=2,nrow=1),
          generate_cowplot(pRnaGeneExp$novelGene_Fgfr1op_AS,pRnaGeneExp$Fgfr1op,num="BC",ncol=2,nrow=1),
          generate_cowplot(pRnaGeneExp$novelGene_Htra1_AS,pRnaGeneExp$Htra1,num="DE",ncol=2,nrow=1), nrow=3, ncol = 1)
dev.off()

pdf(paste0(output_plot_dir,"/DIU_DEG_major.pdf"), width = 10, height = 15)
for(g in 1:10){print(generate_cowplot(plotlist = diff_across_rnaseq(DIU_genes$DIU_DEG_maj$gene[[g]]),num=4,ncol=2,nrow=2))}
dev.off()

pdf(paste0(output_plot_dir,"/DIU_DEG_nomajor.pdf"), width = 10, height = 15)
for(g in 1:10){print(generate_cowplot(plotlist = diff_across_rnaseq(DIU_genes$DIU_DEG_nomaj$gene[[g]]),num=4,ncol=2,nrow=2))}
dev.off()

pdf(paste0(output_plot_dir,"/DIU_notDEG_nomajor.pdf"), width = 10, height = 15)
for(g in 1:10){print(generate_cowplot(plotlist =  diff_across_rnaseq(DIU_genes$DIU_notDEG_nomaj$gene[[g]]),num=4,ncol=2,nrow=2))}
dev.off()

pdf(paste0(output_plot_dir,"/DIU_notDEG_major.pdf"), width = 10, height = 15)
for(g in 1:10){print(generate_cowplot(plotlist = diff_across_rnaseq(DIU_genes$DIU_notDEG_maj$gene[[g]]),num=4,ncol=2,nrow=2))}
dev.off()

write.xlsx(DIU_genes$DIU_DEG_maj, paste0(output_helpfig_dir,"/DIU_Analysis_Genelist.xlsx"), sheetName="DIU_DEG_maj")
write.xlsx(DIU_genes$DIU_DEG_nomaj, paste0(output_helpfig_dir,"/DIU_Analysis_Genelist.xlsx"), sheetName="DIU_DEG_nomaj", append=TRUE)
write.xlsx(DIU_genes$DIU_notDEG_nomaj, paste0(output_helpfig_dir,"/DIU_Analysis_Genelist.xlsx"), sheetName="DIU_notDEG_nomaj", append=TRUE)
write.xlsx(DIU_genes$DIU_notDEG_maj, paste0(output_helpfig_dir,"/DIU_Analysis_Genelist.xlsx"), sheetName="DIU_notDEG_maj", append=TRUE)

##### DFI ##### 
DFI_plots <- DFI_generateplots()


##### Alternative Splicing Events, FDA ############################# 

#### FDA 
FDA_FP_results = FDA_plot(FDA_FP, "Presence")
FDA_GP_results = FDA_plot(FDA_GP, "Position")
FDA_lengths_results = FDA_lengths_plot()

dASevents_files = read_files_differentialASevents(group_sqanti_dir,suppa_dir)
AS_events_diff_plots = AS_events_diff(dASevents_files$group.mono.class.files, dASevents_files$annotated_sqanti_gtf, dASevents_files$suppa2.output)
#IR_ORF_plots = IR_ORF()

pdf(paste0(output_plot_dir,"/WholeFDA.pdf"), width = 10, height = 15)
#grobs <- ggplotGrob(AS_events_diff_plots[[1]])$grobs
#legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
#plot_grid(AS_events_diff_plots[[1]] + theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1)),AS_events_diff_plots[[3]],labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman",NULL,NULL, nrow = 2, scale = 0.9)
#plot_grid(IR_ORF_plots[[1]],IR_ORF_plots[[2]],labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)
#plot_grid(IR_ORF_plots[[3]],IR_ORF_plots[[4]],labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", nrow = 2, scale = 0.9)

# FDA
top <- plot_grid(FDA_FP_results[[2]] + theme(legend.position = "none"),FDA_GP_results[[2]], align='vh', vjust=1, scale = 1, labels = c("a","b"),label_size = 30, label_fontfamily = "CM Roman")
y.grob <- textGrob("Features", gp=gpar(fontfamily="CM Roman", fontsize=16), rot=90)
x.grob <- textGrob("Number of Genes with Features", gp=gpar(fontfamily="CM Roman", fontsize=16))
bottom <- plot_grid(FDA_lengths_results,NULL,ncol=2,labels = c("c",""), label_size = 30, label_fontfamily = "CM Roman")
plot_grid(grid.arrange(arrangeGrob(top, left = y.grob, bottom = x.grob)),bottom, nrow = 3, scale = 0.9,rel_heights = c(0.5,0.4,0.1))
# DFI
plot_grid(DFI_plots[[2]],DFI_plots[[3]],NULL,NULL, labels = c("a","b"),label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
dev.off()

##### Summary Dataset ############################# 
summaryoutput = summary_dataset() # Global table 

#### Target Genes 
pdf(paste0(output_plot_dir,"/WholeDifferentialAnalysis_TargetGenes.pdf"), width = 10, height = 15)
for(gene in TargetGene){
  print(plot_grid(pIsoGeneExp[[gene]],pRnaGeneExp[[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
  print(plot_grid(pIsoTranExpTraj [[gene]],pRnaTranExpTraj [[gene]],ncol = 1,labels = "auto", label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
}
dev.off()

#### WGCNA
venn_WGCNA()

#### Methylation Integration #### 
# List of Genes identified with DMR/DMP (from Isabel's analysis)
DMPGenotype_Integration = Methylation_Integration("DMP_Genotype")
DMPInteraction_Integration = Methylation_Integration("DMP_Interaction")
DMRGenotype_Integration = Methylation_Integration("DMR_Genotype")

DMP_Combine = rbind(DMPGenotype_Integration$Positions,DMPInteraction_Integration$Positions) 
DMP_Combine = DMP_Combine[!duplicated(DMP_Combine), ]
length(unique(DMP_Combine$Chipseeker_SYMBOL))

# Stats 
Meth_Table = lapply(unique(DMP_Combine$Chipseeker_SYMBOL), function(x) Methylation_Integration_stats(x))
Meth_Table <- do.call("rbind",Meth_Table) 
Meth_Table %>% mutate(loc_simple = word(location,c(1),sep = fixed(" "))) %>% group_by(loc_simple) %>% tally()
binom.test(14, 18, 0.5)
binom.test(5, 18, 0.5)

plot_transexp_overtime("Osmr",RNAExp$Norm_transcounts,"isoseq","NA")
# Methyylation and Expression plots 
Spata13_output = Methylation_Integration_plots("Spata13","PB.4966.2")
Osmr_output = Methylation_Integration_plots("Osmr","PB.5258.1")
Jph1_output = Methylation_Integration_plots("Jph1","NA")
Ncf2_output = Methylation_Integration_plots("Ncf2","PB.700.1")
IRF8_output = Methylation_Integration_plots("Irf8","PB.15969.1")
Susd5_output = Methylation_Integration_plots("Susd5","PB.16983.1")
As3mt_methylation_output = as3mt_dmr_figure()
Prnp_methylation_output = prnp_dmr_figure()

# Plot Transcript Expression Trajection of DMP/DMR Genes
DMP_DMR_Genes = unique(c(DMPGenotype_Integration$Positions$Chipseeker_SYMBOL,
         DMPInteraction_Integration$Positions$Chipseeker_SYMBOL,
         as.character(DMRGenotype_Integration$Positions$Chipseeker_SYMBOL)))

DMP_DMR_Genes_plots_traj  <- lapply(lapply(DMP_DMR_Genes, function(gene) 
  plot_transexp_overtime(gene,RNAExp$Norm_transcounts,"isoseq","RNA-Seq Isoform Expression") +
    theme(legend.position = "bottom", plot.title = element_blank())),ggplotGrob)
names(DMP_DMR_Genes_plots_traj) <- DMP_DMR_Genes

#pdf(paste0(output_plot_dir,"/WholeDifferentialAnalysis_DMPDMR.pdf"), width = 10, height = 5)
#for(i in DMP_DMR_Genes){print(plot_grid(DMP_DMR_Genes_plots_traj[[i]]))}
#dev.off()

pdf(paste0(output_plot_dir,"/WholeDifferentialAnalysis_DMPDMR.pdf"), width = 15, height = 8)
generate_cowplot(Prnp_methylation_output[[1]],Prnp_methylation_output[[2]],Prnp_methylation_output[[3]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$As3mt,As3mt_methylation_output[[1]],As3mt_methylation_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Spata13,Spata13_output[[1]],Spata13_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Ncf2,Ncf2_output[[1]],Ncf2_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Irf8,IRF8_output[[1]],IRF8_output[[2]],num=3,ncol=3,nrow=1)
generate_cowplot(DMP_DMR_Genes_plots_traj$Susd5,Susd5_output[[1]],Susd5_output[[2]],num=3,ncol=3,nrow=1)

generate_cowplot(DMP_DMR_Genes_plots_traj$Osmr,Osmr_output[[1]],Osmr_output[[2]],num=3,ncol=3,nrow=1)
plot_grid(DMP_DMR_Genes_plots_traj$Spata13,Spata13_output[[1]],Spata13_output[[2]], ncol = 3,
          labels = "auto",label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(DMP_DMR_Genes_plots_traj$Ncf2,Ncf2_output[[1]],Ncf2_output[[2]], ncol = 3,
          labels = "auto",label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(DMP_DMR_Genes_plots_traj$Susd5,Susd5_output[[1]],Susd5_output[[2]], ncol = 3,
          labels = "auto",label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(DMP_DMR_Genes_plots_traj$Irf8,IRF8_output[[1]],IRF8_output[[2]], ncol = 3,
          labels = "auto",label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
dev.off()

# differentially expressed isoforms
tlist = c("PB.700.1","PB.15969.1")
meanexp_output = list()
meanexp_output = lapply(tlist, function(x) tabulate_difftrans(x, IsoExp$Norm_transcounts, All_FDR,"isoseq"))
names(meanexp_output) = tlist
meanexp_output <- do.call("rbind",meanexp_output) 
meanexp_output

View(DMP_Combine[DMP_Combine$Chipseeker_SYMBOL == "Irf8",])
