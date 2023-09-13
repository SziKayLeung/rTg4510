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



## ---------- Source function and config files -----------------

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(SC_ROOT, "rTg4510_config.R"))
#source(paste0(SC_ROOT,"bin/draw_heatmap_gene_level.R"))


## ---------- Global Iso-Seq dataset -----------------

# mean number of genes
cat("Mean number of genes:", mean(sapply(sub_class.files, function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("StDev number of genes:", sd(sapply(sub_class.files, function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("Mean number of genes in WT samples:", mean(sapply(sub_class.files[wholeWT], function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("StDev number of genes in WT samples:", sd(sapply(sub_class.files[wholeWT], function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("Mean number of genes in TG samples:", mean(sapply(sub_class.files[wholeTG], function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("StDev number of genes in TG samples:", sd(sapply(sub_class.files[wholeTG], function(x) length(unique(x[["associated_gene"]])))),"\n")

cat("Mean number of genes in 2 months samples:", mean(sapply(sub_class.files[whole2mos], function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("StDev number of genes in 2 months samples:", sd(sapply(sub_class.files[whole2mos], function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("Mean number of genes in 8 months samples:", mean(sapply(sub_class.files[whole8mos], function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("StDev number of genes in 8 months samples:", sd(sapply(sub_class.files[whole8mos], function(x) length(unique(x[["associated_gene"]])))),"\n")

# mean number of isoforms
cat("Mean number of isoforms:", mean(sapply(sub_class.files, function(x) nrow(x))),"\n")
cat("StDev number of isoforms:", sd(sapply(sub_class.files, function(x) nrow(x))),"\n")
cat("Mean number of isoforms in WT samples:", mean(sapply(sub_class.files[wholeWT], function(x) nrow(x))),"\n")
cat("StDev number of isoforms in WT samples:", sd(sapply(sub_class.files[wholeWT], function(x) nrow(x))),"\n")
cat("Mean number of isoforms in TG samples:", mean(sapply(sub_class.files[wholeTG], function(x) nrow(x))),"\n")
cat("StDev number of isoforms in TG samples:", sd(sapply(sub_class.files[wholeTG], function(x) nrow(x))),"\n")

cat("Mean number of isoforms in 2 months samples:", mean(sapply(sub_class.files[whole2mos], function(x) nrow(x))),"\n")
cat("StDev number of isoforms in 2 months samples:", sd(sapply(sub_class.files[whole2mos], function(x) nrow(x))),"\n")
cat("Mean number of isoforms in 8 months samples:", mean(sapply(sub_class.files[whole8mos], function(x) nrow(x))),"\n")
cat("StDev number of isoforms in 8 months samples:", sd(sapply(sub_class.files[whole8mos], function(x) nrow(x))),"\n")


# mean descriptives
cat("Mean isoform length in all samples:", mean(class.files$glob_iso$length), "\n")
cat("StDev isoform length in all samples:", sd(class.files$glob_iso$length), "\n")
cat("Mean number of exons in all samples:", mean(class.files$glob_iso$exons), "\n")
cat("StDev number of exons in all samples:", sd(class.files$glob_iso$exons), "\n")
cat("Mean number of exons in all samples:", mean(as.data.frame(class.files$glob_iso %>% group_by(associated_gene) %>% tally())$n), "\n")

# Global data 
for(i in c("Gfap","C4b")){
  cat("***************** Gene level\n")
  cat("RNA-Seq gene level statistics for", i, "\n")
  print(GlobalDESeq$RresGeneAnno$wald$anno_res %>% filter(associated_gene == i))
  print(GlobalDESeq$RresGeneAnno$wald$stats_Wald %>% filter(associated_gene == i) %>% select(associated_gene, WaldStatistic_groupCASE.time8))
  cat("Iso-Seq gene level statistics for", i, "\n")
  print(GlobalDESeq$resGeneAnno$wald$anno_res %>% filter(associated_gene == i))
  print(GlobalDESeq$resGeneAnno$wald$stats_Wald %>% filter(associated_gene == i) %>% select(associated_gene, WaldStatistic_groupCASE.time8))
}

for(i in c("PB.2973.16","PB.7022.9")){
  cat("***************** Transcript level\n")
  cat("RNA-Seq gene level statistics for", i, "\n")
  #print(GlobalDESeq$RresTranAnno$wald$anno_res %>% filter(isoform == i))
  #print(GlobalDESeq$RresTranAnno$wald$stats_LRT %>% filter(rownames(.) == i) %>% select(groupCASE.time8, LRTPvalue, LRTStatistic))
  cat("Iso-Seq gene level statistics for", i, "\n")
  print(GlobalDESeq$resTranAnno$wald$anno_res %>% filter(isoform == i))
  print(GlobalDESeq$resTranAnno$wald$stats_Wald %>% filter(rownames(.) == i) %>% select(WaldPvalue_groupCASE.time8))
}

cat("Total number of reads in filterd dataset (ONT and Iso-Seq combined):", round(sum(class.files$targ_filtered$nreads)/1000000,2),"million \n")
TargetedMergedCouts <- list(
  WT = class.files$targ_filtered %>% select(contains(targetedWT)) %>% reshape2::melt() %>% mutate(dataset = word(variable,c(1),sep=fixed("_"))),
  TG = class.files$targ_filtered %>% select(contains(targetedTG)) %>% reshape2::melt() %>% mutate(dataset = word(variable,c(1),sep=fixed("_")))
)
lapply(TargetedMergedCouts, function(x) x %>% group_by(dataset, variable) %>% tally(value) %>% group_by(dataset) %>% dplyr::summarize(mean_value = mean(n)))
lapply(TargetedMergedCouts, function(x) x %>% group_by(dataset) %>% dplyr::summarize(mean_value = mean(value)))
lapply(TargetedMergedCouts, function(x) x %>% group_by(variable) %>% tally(value) %>% dplyr::summarize(mean_value = mean(n)))

TargetedMergedCoutsc <- lapply(TargetedMergedCouts, function(x) x %>% group_by(dataset, variable) %>% tally(value) %>% group_by(dataset)) %>% bind_rows(.id = "group")
TargetedMergedCoutsc %>% group_by(group) %>% dplyr::summarize(mean_value = mean(n))
wilcox.test(n~ group, data = TargetedMergedCoutsc[TargetedMergedCoutsc$dataset == "Iso.Seq",], exact = FALSE)
wilcox.test(n~ group, data = TargetedMergedCoutsc[TargetedMergedCoutsc$dataset == "ONT",], exact = FALSE)
wilcox.test(n~ group, data = TargetedMergedCoutsc, exact = FALSE)

cat("Number of matched samples in Whole vs Targeted Transcriptome:", length(wholesamples), "\n")
wholevsTargeted <- whole_vs_targeted_plots(class.files$iso_match,paste0("FL.WholeIso", wholesamples), paste0("FL.TargetedIso", wholesamples), TargetGene)[[3]]
cat("Number of isoforms detected by both whole and targeted, whole only and unique only\n:")
wholevsTargetedTally <- wholevsTargeted %>% group_by(dataset) %>% tally()
wholevsTargetedTally
cat("Ratio of targeted to whole transcriptome of detecting transcripts associated to target genes:", 
    as.numeric(sum(wholevsTargetedTally[wholevsTargetedTally$dataset != "Whole","n"])/
                 sum(wholevsTargetedTally[wholevsTargetedTally$dataset != "Targeted","n"])))
# wholevsTargeted %>% group_by(dataset) %>% tally(sumWhole)
# wholevsTargeted %>% group_by(dataset) %>% tally(sumTargeted)

# fisher's exact test
data <- matrix(c(815, 440, 15, 440), nrow = 2)
rownames(data) <- c("Unique_Transcripts", "Not_Unique_Transcripts")
colnames(data) <- c("Targeted_Sequencing", "Whole_Transcriptome_Sequencing")
res <- fisher.test(data)
res$p.value

recapitulate_gene_level()

## ---------- Merged Targeted Iso-Seq and ONT datasets -----------------

cat("Total number of transcripts before filtering:", nrow(class.files$targ_all),"\n")
cat("Total number of unique ONT transcripts before filtering:", class.files$targ_all %>% filter(Dataset == "ONT") %>% nrow(),"\n")
cat("% of unique ONT transcripts before filtering:", class.files$targ_all %>% filter(Dataset == "ONT") %>% nrow()/nrow(class.files$targ_all),"\n")
cat("Total number of Iso-Seq transcripts before filtering:", class.files$targ_all %>% filter(Dataset != "ONT") %>% nrow(), "\n")
cat("Total number of Iso-Seq transcripts also detected in ONT before filtering:", class.files$targ_all %>% filter(Dataset == "Both") %>% nrow(), "\n")
cat("% of Iso-Seq transcripts detected in ONT", class.files$targ_all %>% filter(Dataset == "Both") %>% nrow()/class.files$targ_all %>% filter(Dataset != "ONT") %>% nrow())
cat("% of ONT transcripts detected in Iso-Seq", class.files$targ_all %>% filter(Dataset == "Both") %>% nrow()/class.files$targ_all %>% filter(Dataset != "Iso-seq") %>% nrow())
    
nrow(class.files$targ_all)
nrow(class.files$targ_all %>% filter(associated_transcript != "novel"))
nrow(class.files$targ_all %>% filter(associated_transcript == "novel"))

nrow(class.files$targ_filtered)
nrow(class.files$targ_filtered %>% filter(associated_transcript != "novel"))
nrow(class.files$targ_filtered %>% filter(associated_transcript == "novel"))

class.files$targ_filtered %>% filter(Dataset == "ONT") %>% nrow()
class.files$targ_filtered %>% filter(Dataset == "Iso-Seq") %>% nrow()
class.files$targ_filtered %>% filter(Dataset == "Both") %>% nrow()

nrow(class.files$targ_all %>% filter(nreads < 10))

# summary
all_summarise_gene_stats(Targeted$Gene_class, class.files$targ_filtered, Targeted$cpat, Targeted$noORF, Targeted$Genes)

# median number of isoforms 
median(as.data.frame(class.files$targ_filtered %>% group_by(associated_gene) %>% tally())$n)
class.files$targ_filtered %>% group_by(associated_gene) %>% tally() %>% arrange(n)
class.files$targ_filtered %>% mutate(noveltrans = ifelse(associated_transcript == "novel","novel","known")) %>% 
  group_by(associated_gene, noveltrans) %>% tally() %>% arrange(n)

# isoform fraction 
statsIFAll <- plotIFAll(Exp=Exp$targ_ont$normAll %>% select(-associated_gene),
          classf=class.files$targ_all,
          pheno=phenotype$targeted_rTg4510_ont,
          majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso))[2][[1]]

statsIFAllMajor <- statsIFAll %>% filter(structural_category != "minor") %>% group_by(associated_gene) %>% tally()
mean(statsIFAllMajor$n)
sd(statsIFAllMajor$n)

# Targeted stats 
totaln = sum(Merged_gene_class_df["Total.Number.of.Transcripts",])
cat("Number of A5A3 events across 20 target genes:", sum(Merged_gene_class_df["A5A3",]), "\n")
cat("Number of trancripts with A5A3' events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.A5A3.Events",]), "\n")
cat("% of trancripts with A5A3' events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.A5A3.Events",])/totaln * 100, "\n")
cat("Number of ES events across 20 target genes:", sum(Merged_gene_class_df["ES",]), "\n")
cat("Number of trancripts with ES events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.ES.Events",]), "\n")
cat("Proportion of trancripts with ES events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.ES.Events",])/totaln * 100, "\n")
cat("Number of ES events across 20 target genes:", sum(Merged_gene_class_df["IR",]), "\n")
cat("Number of IR events across 20 target genes:", sum(Merged_gene_class_df["IR",]), "\n")
cat("Number of trancripts with IR events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.IR.Events",]), "\n")
cat("Proportion of trancripts with IR events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.IR.Events",])/totaln * 100, "\n")

IR_tab_exonoverlap <- input_FICLE_splicing_results(dirnames$targ_anno, "IntronRetentionExonOverlap")
cat("Number of IR events overlapping 2 or more exons:", nrow(IR_tab_exonoverlap[IR_tab_exonoverlap$IRNumExonsOverlaps >= 2,]), "\n") 
cat("Number of IR events overlapping 5 exons:", nrow(IR_tab_exonoverlap[IR_tab_exonoverlap$IRNumExonsOverlaps == 5,]), "\n") 
cat("Proportion of IR events overlapping 2 or more exons:", nrow(IR_tab_exonoverlap[IR_tab_exonoverlap$IRNumExonsOverlaps >= 2,])/sum(Merged_gene_class_df["IR",]) * 100, "\n") 


# Apoe
subset(Targeted$ref_gencode, associated_gene == "Apoe")["MaxGeneLength"]
subset(Targeted$ref_gencode, associated_gene == "Apoe")["Maxexons"]

## ---------- Differential isoform expression analyis - Targeted -----------------

# Differentially expressed across pathology
reportStats(res=TargetedDESeq$ontResTranAnno$wald$anno_res,stats=TargetedDESeq$ontResTranAnno$wald$stats_Wald, isoList=c("PB.14646.139","PB.20818.54"))
unique(TargetedDESeq$ontResTranAnno$wald$anno_res$associated_gene)

# isoform detected in wald across age, and not at 8 months
setdiff(TargetedDESeq$ontResTranAnno$wald$anno_res$isoform,TargetedDESeq$ontResTranAnno$wald8mos$anno_res$isoform)

# isoform detected in wald at 8 months and not across age
# further detection
addedWald8mosIso = setdiff(TargetedDESeq$ontResTranAnno$wald8mos$anno_res$isoform, TargetedDESeq$ontResTranAnno$wald$anno_res$isoform)
cat("Further number of transcripts identified at 8 months:", length(addedWald8mosIso), "\n")
cat("Annotated to:", unique(TargetedDESeq$ontResTranAnno$wald8mos$anno_res %>% filter(isoform %in% addedWald8mosIso) %>% .[,c("associated_gene")],"\n"))
reportStats(res=TargetedDESeq$ontResTranAnno$wald8mos$res_Wald,stats=TargetedDESeq$ontResTranAnno$wald8mos$stats_Wald, isoList=c("PB.22007.99"))
reportStats(res=TargetedDESeq$ontResTranAnno$wald8mos$anno_res,stats=TargetedDESeq$ontResTranAnno$wald8mos$stats_Wald, 
            isoList=c("PB.8675.37810","PB.22007.224","PB.38419.87","PB.19309.7497","PB.40586.875"))

# Differentially express across genotype 
TargetedDESeq$ontResTranAnno$waldgenotype$anno_res <- TargetedDESeq$ontResTranAnno$waldgenotype$anno_res %>% filter(padj < 0.05)
nrow(TargetedDESeq$ontResTranAnno$waldgenotype$anno_res)

# Iso-Seq Trem2 results 
# genotype
reportStats(TargetedDESeq$isoResTranAnno$waldgenotype$anno_res,TargetedDESeq$isoResTranAnno$waldgenotype$stats_Wald,c("PB.20818.62","PB.20818.54"))
# pathology
reportStats(TargetedDESeq$isoResTranAnno$wald$anno_res,TargetedDESeq$isoResTranAnno$wald$stats_Wald,c("PB.20818.62","PB.20818.54"))
TargetedDESeq$isoResTranAnno$wald$stats_Wald %>% filter(row.names(.) == "PB.20818.62") %>% select(groupCASE.time8)


# remove transcripts already differentially expressed across genotype
genotypeOnlyIso = setdiff(TargetedDESeq$ontResTranAnno$waldgenotype$anno_res$isoform,TargetedDESeq$ontResTranAnno$wald$anno_res$isoform)
intersect(TargetedDESeq$ontResTranAnno$waldgenotype$anno_res$isoform,TargetedDESeq$ontResTranAnno$wald$anno_res$isoform)
genotypeOnly = TargetedDESeq$ontResTranAnno$waldgenotype$anno_res %>% filter(isoform %in% genotypeOnlyIso)

nrow(genotypeOnly)
length(unique(genotypeOnly$associated_gene))
sort(unique(genotypeOnly$associated_gene))

TargetedMergedDESeq$waldGenotype %>% filter(isoform %in% genotypeOnlyIso) %>% group_by(ont_direction) %>% tally()
TargetedMergedDESeq$waldGenotype %>% filter(isoform %in% genotypeOnlyIso) %>% group_by(associated_gene) %>% tally()
TargetedMergedDESeq$waldGenotype %>% filter(isoform %in% genotypeOnlyIso) %>% group_by(structural_category) %>% tally()
nrow(TargetedMergedDESeq$waldGenotype %>% filter(isoform %in% genotypeOnlyIso) %>% filter(!is.na(isoseq_direction)))



# Replication with Iso-Seq
# Calculate the adjusted threshold 
num_tests = nrow(TargetedDESeq$ontResTranAnno$wald$anno_res)
threshold <- 0.05 / num_tests
TargetedDESeq$isoResTranAnno$wald$res$significant <- TargetedDESeq$isoResTranAnno$wald$res$pvalue <= threshold
replicatedIso = intersect(TargetedDESeq$isoResTranAnno$wald$res[TargetedDESeq$isoResTranAnno$wald$res$significant == TRUE,"isoform"], 
          TargetedDESeq$ontResTranAnno$wald$anno_res$isoform)
setdiff(TargetedDESeq$ontResTranAnno$wald$anno_res$isoform,TargetedDESeq$isoResTranAnno$wald$res[TargetedDESeq$isoResTranAnno$wald$res$significant == TRUE,"isoform"])

# effect size = group_CASE_vs_CONTROL 
replicatedVals = rbind(TargetedDESeq$isoResTranAnno$wald$stats_Wald %>% rownames_to_column(var = "isoform") 
                       %>% filter(isoform %in% replicatedIso) %>% .[,c("isoform","group_CASE_vs_CONTROL")] %>% mutate(dataset = "IsoSeq"),
      TargetedDESeq$ontResTranAnno$wald$stats_Wald %>% rownames_to_column(var = "isoform") %>% 
        filter(isoform %in% replicatedIso) %>% .[,c("isoform","group_CASE_vs_CONTROL")] %>% mutate(dataset = "ONT")) %>% 
  spread(., dataset,group_CASE_vs_CONTROL) 
ggplot(replicatedVals, aes(x = ONT, y = IsoSeq)) + geom_point()
cor.test(replicatedVals$IsoSeq, replicatedVals$ONT, method = "spearman")

SigPathologyVals = rbind(TargetedDESeq$isoResTranAnno$wald$stats_Wald %>% rownames_to_column(var = "isoform") 
                       %>% filter(isoform %in% TargetedDESeq$ontResTranAnno$wald$anno_res$isoform) %>% .[,c("isoform","group_CASE_vs_CONTROL")] %>% mutate(dataset = "IsoSeq"),
                       TargetedDESeq$ontResTranAnno$wald$stats_Wald %>% rownames_to_column(var = "isoform") %>% 
                         filter(isoform %in% TargetedDESeq$ontResTranAnno$wald$anno_res$isoform) %>%
                         .[,c("isoform","group_CASE_vs_CONTROL")] %>% mutate(dataset = "ONT")) %>% 
  spread(., dataset,group_CASE_vs_CONTROL) 
ggplot(SigPathologyVals, aes(x = ONT, y = IsoSeq)) + geom_point()
cor.test(SigPathologyVals$IsoSeq, SigPathologyVals$ONT, method = "spearman")

SigGenotypeVals = rbind(TargetedDESeq$isoResTranAnno$waldgenotype$stats_Wald %>% rownames_to_column(var = "isoform") %>% 
                          filter(isoform %in% TargetedDESeq$ontResTranAnno$waldgenotype$anno_res$isoform) %>% 
                          .[,c("isoform","group_CASE_vs_CONTROL")] %>% mutate(dataset = "IsoSeq"),
                         TargetedDESeq$ontResTranAnno$waldgenotype$stats_Wald %>% rownames_to_column(var = "isoform") %>% 
                          filter(isoform %in% TargetedDESeq$ontResTranAnno$waldgenotype$anno_res$isoform) %>%
                           .[,c("isoform","group_CASE_vs_CONTROL")] %>% mutate(dataset = "ONT")) %>% 
  spread(., dataset,group_CASE_vs_CONTROL) 
ggplot(SigGenotypeVals, aes(x = ONT, y = IsoSeq)) + geom_point()
cor.test(SigGenotypeVals$IsoSeq, SigGenotypeVals$ONT, method = "spearman")


res <- binom.test(consistentNum, totalNum, alternative = c("two.sided"))
res$p.value

# Bin1 224
Bin1224_IF <- plotIF("Bin1",ExpInput=Exp$targ_ont$normAll,pheno=phenotype$targeted_rTg4510_ont,
                     cfiles=class.files$targ_all,design="time_series",majorIso=NULL,isoSpecific=c("PB.22007.224"),rank=0,stat=TRUE) %>% filter(isoform == "PB.22007.224")

aggregate(perc ~ group, Bin1224_IF, mean)
aggregate(perc ~ group + time, Bin1224_IF, mean)
wilcox.test(perc ~ group, data = Bin1224_IF, exact = FALSE)
wilcox.test(perc ~ time, data = Bin1224_IF[Bin1224_IF$group == "TG" & Bin1224_IF$Age %in% c(2,8),], exact = FALSE)

## ---------- Differential isoform usage analysis - Targeted -----------------

isoUsage = c("PB.42931.201","PB.22007.224","PB.40586.875","PB.40586.872","PB.19309.7368","PB.19309.7374","PB.19309.7389")
TargetedDESeq$ontResTranAnno$wald$res_Wald %>% filter(isoform %in% c("PB.42931.201","PB.22007.224","PB.40586.875","PB.40586.872"))
TargetedDESeq$ontResTranAnno$wald8mos$res_Wald %>% filter(isoform %in% isoUsage)
class.files$targ_all[class.files$targ_all$isoform %in% isoUsage,c("isoform","associated_gene","structural_category","associated_transcript","subcategory")]
TargetedDIU$ontDIUGeno$resultDIU %>% filter(Gene %in% c("Fus","Bin1","Trpa1","Apoe","App"))
TargetedDIU$ontDIUPath$resultDIU %>% filter(Gene %in% c("Fus","Bin1","Apoe"))

## ---------- Trem2 ----------
# Read in FICLE files
Trem2FICLE <- input_FICLE_all_perGene(dirnames$targ_anno,"Trem2")

# Numbers
nTrem2 <- Merged_gene_class_df["Total.Number.of.Transcripts","Trem2"]
nTrem2NoExon3 <- nrow(Trem2FICLE$Exonskipping_tab %>% filter(ES == "Gencode_3"))
nTrem2NoExon2 <- nrow(Trem2FICLE$Exonskipping_tab %>% filter(ES == "Gencode_2"))
nTrem2A5A3 <- length(unique(Trem2FICLE$A5A3_tab$transcript_id))
nTrem2NovelExon <- length(unique(Trem2FICLE$NE$transcript_id))
nTrem2NovelExonBeyond <- Trem2FICLE$NE_counts[Trem2FICLE$NE_counts$novelexons == "Beyond_First","Cat"]

cat("Total number of Trem2 transcripts:",nTrem2,"\n")
cat("Number of Trem2 transcripts with no exon 3:",nTrem2NoExon3,"(",round(nTrem2NoExon3/nTrem2 * 100,2),"%)\n")
cat("Number of Trem2 transcripts with A5A3:",nTrem2A5A3,"(",round(nTrem2A5A3/nTrem2 * 100,2),"%)\n")
cat("Number of Trem2 transcripts with no exon 2:",nTrem2NoExon2,"(",round(nTrem2NoExon2/nTrem2 * 100,2),"%)\n")
cat("Number of Trem2 transcripts with A5A3:",nTrem2A5A3,"(",round(nTrem2A5A3/nTrem2 * 100,2),"%)\n")
cat("Number of Trem2 transcripts with novel exons:",nTrem2NovelExon,"(",round(nTrem2NovelExon/nTrem2 * 100,2),"%)\n")
cat("Number of Trem2 transcripts with novel exons:",nTrem2NovelExonBeyond,"(",round(nTrem2NovelExonBeyond/nTrem2 * 100,2),"%)\n")

# Differential expression gene level 
reportStats(res=TargetedDESeq$ontResGeneAnno$wald$res_Wald, stats=TargetedDESeq$ontResGeneAnno$wald$stats_Wald, isoList=c("20818"))

## ---------- Bin1 ----------
# Read in FICLE files
Bin1FICLE <- input_FICLE_all_perGene(dirnames$targ_anno,"Bin1")

# Numbers
nBin1 <- Merged_gene_class_df["Total.Number.of.Transcripts","Bin1"]
ES <- input_FICLE_splicing_results(dirnames$targ_anno,"Exonskipping_tab")
Bin1ESMax <- ES %>% filter(associated_gene == "Bin1") %>% group_by(transcript_id) %>% tally() %>% filter(n > 10)
cat("Number of Bin1 transcripts with more than 10 exons skipped:",length(Bin1ESMax$transcript_id),"(",round(length(Bin1ESMax$transcript_id)/nBin1*100,2),"%)\n")


cat("Total number of Bin1 transcripts:",nBin1,"\n")
for(i in c(14,15,16)){cat("Number of Bin1 transcripts with no exon",i,":", nrow(Bin1FICLE$Exonskipping_tab %>% filter(ES == paste0("Gencode_",i))),"\n")}
cat("Number of transcripts with IR in CLAP domain:",nrow(unique(Bin1FICLE$IntronRetention_tab %>% mutate(exon = word(IR,c(3),sep=fixed("_"))) %>% filter(exon%in%c(13,14,15,16)) %>% select("transcript_id"))),"\n")

# dominant isoform = novel isoform
#View(class.files$targ_all[class.files$targ_all$associated_gene == "Bin1",c("isoform","nreads","ONT_sum_FL","Iso.Seq_sum_FL")])
#sum(class.files$targ_all[class.files$targ_all$associated_gene == "Bin1","ONT_sum_FL"])
36593/sum(class.files$targ_all[class.files$targ_all$associated_gene == "Bin1","ONT_sum_FL"]) * 100
16624/sum(class.files$targ_all[class.files$targ_all$associated_gene == "Bin1","ONT_sum_FL"]) * 100
8739/sum(class.files$targ_all[class.files$targ_all$associated_gene == "Bin1","Iso.Seq_sum_FL"]) * 100
4105/sum(class.files$targ_all[class.files$targ_all$associated_gene == "Bin1","Iso.Seq_sum_FL"]) * 100

Bin1Isoforms = c("PB.22007.101","PB.22007.224","PB.22007.99","PB.22007.861","PB.22007.45098","PB.22007.176")
class.files$targ_all[class.files$targ_all$isoform == "PB.22007.101",c("ONT_sum_FL","Iso.Seq_sum_FL")]
class.files$targ_all[class.files$targ_all$isoform %in% Bin1Isoforms, c("structural_category","associated_gene","associated_transcript","subcategory")]

# Differential transcript and gene expression statistics
reportStats(res=TargetedDESeq$ontResTranAnno$wald8mos$anno_res, stats=TargetedDESeq$ontResTranAnno$wald8mos$stats_Wald, isoList=c("PB.22007.224"))
reportStats(res=TargetedDESeq$isoResTranAnno$wald8mos$anno_res, stats=TargetedDESeq$isoResTranAnno$wald8mos$stats_Wald, isoList=c("PB.22007.224"))
reportStats(res=TargetedDESeq$ontResGeneAnno$wald$res_Wald, stats=TargetedDESeq$ontResGeneAnno$wald$stats_Wald, isoList=c("22007"))

## ---------- Clu ----------
# Read in FICLE files
CluFICLE <- input_FICLE_all_perGene(dirnames$targ_anno,"Clu")

# Numbers
nClu <- Merged_gene_class_df["Total.Number.of.Transcripts","Clu"]
CluFirstExon <- CluFICLE$Exonskipping_generaltab %>% filter(Gencode_1 == "FirstExon")
nCluFirstExon_SkippedAF <- nrow(CluFICLE$Exonskipping_tab %>% filter(ES %in% c("Gencode_2","Gencode_3","Gencode_4","Gencode_5","Gencode_6")) 
                               %>% group_by(transcript_id) %>% tally() %>% filter(n == 5))
for(i in c(9,10,11)){cat("Number of Clu transcripts with no exon",i,":", nrow(CluFICLE$Exonskipping_tab %>% filter(ES == paste0("Gencode_",i))),"\n")}
cat("Total number of Clu transcripts:",nClu,"\n")
cat("Number of Clu transcripts with Gencode 1 as first exon:",nrow(CluFirstExon),"(",round(nrow(CluFirstExon)/nClu*100,2),"%)\n")

# Differential transcript and gene expression statistics
CluIso = c("PB.14646.139","PB.14646.60837","PB.14646.68849","PB.14646.39341")
class.files$targ_all[class.files$targ_all$isoform %in% CluIso, c("structural_category","associated_gene","associated_transcript","subcategory")]
reportStats(res=TargetedDESeq$ontResTranAnno$wald8mos$res_Wald, stats=TargetedDESeq$ontResTran$wald8mos$stats_Wald, isoList=CluIso)
reportStats(res=TargetedDESeq$ontResTranAnno$wald$res_Wald, stats=TargetedDESeq$ontResTran$wald$stats_Wald, isoList=CluIso)
reportStats(res=TargetedDESeq$ontResGeneAnno$wald$res_Wald, stats=TargetedDESeq$ontResGeneAnno$wald$stats_Wald, isoList=c("14646"))
reportStats(res=TargetedDESeq$ontResGeneAnno$waldgenotype$res_Wald, stats=TargetedDESeq$ontResGeneAnno$waldgenotype$stats_Wald, isoList=c("14646"))

## ---------- App ----------
# Read in FICLE files
AppFICLE <- input_FICLE_all_perGene(dirnames$targ_anno,"App")

# Numbers
ES <- input_FICLE_splicing_results(dirnames$targ_anno,"Exonskipping_tab")
nApp <- Merged_gene_class_df["Total.Number.of.Transcripts","App"]
nAppES <- length(unique(AppFICLE$Exonskipping_tab$transcript_id))
AppESMax <- ES %>% filter(associated_gene == "App") %>% group_by(transcript_id) %>% tally() %>% filter(n > 10)
cat("Number of App transcripts with more than 10 exons skipped:",length(AppESMax$transcript_id),"(",round(length(AppESMax$transcript_id)/nApp*100,2),"%)\n")
cat("Number of App transcripts with exon skipping:",nAppES,"(",round(nAppES/nApp*100,2),"%)\n")
for(i in c(7,8,14,15,17)){cat("Number of App transcripts with no exon",i,":", nrow(AppFICLE$Exonskipping_tab %>% filter(ES == paste0("Gencode_",i))),"\n")}
AppIso = group_class.files.diff.targeted[group_class.files.diff.targeted$associated_gene == "App",]
nrow(AppIso[AppIso$WTFL == 0,])
nrow(AppIso[AppIso$TGFL == 0,])

# Differential transcript expression statistics
AppIso = c("PB.19309.7368","PB.19309.7374","PB.19309.7389")
class.files$targ_all[class.files$targ_all$isoform %in% AppIso, c("structural_category","associated_gene","associated_transcript","subcategory")]
reportStats(res=TargetedDESeq$ontResTranAnno$wald8mos$res_Wald, stats=TargetedDESeq$ontResTran$wald8mos$stats_Wald, isoList=AppIso)
reportStats(res=TargetedDESeq$ontResTranAnno$wald$res_Wald, stats=TargetedDESeq$ontResTran$wald$stats_Wald, isoList=AppIso)
reportStats(res=TargetedDESeq$ontResTranAnno$waldgenotype$res_Wald, stats=TargetedDESeq$ontResTran$waldgenotype$stats_Wald, isoList=AppIso)
reportStats(res=TargetedDESeq$ontResGeneAnno$wald$res_Wald, stats=TargetedDESeq$ontResGeneAnno$wald$stats_Wald, isoList=c("19309"))
reportStats(res=TargetedDESeq$ontResGeneAnno$waldgenotype$res_Wald, stats=TargetedDESeq$ontResGeneAnno$waldgenotype$stats_Wald, isoList=c("19309"))

# Differential transcript usage
TargetedDIU$ontDIUGeno$resultDIU %>% filter(Gene == "App")


## ---------- Apoe ----------

# Differential transcript expression statistics
ApoeIso = c("PB.40586.871","PB.40586.872","PB.40586.875")
class.files$targ_all[class.files$targ_all$isoform %in% ApoeIso, c("structural_category","associated_gene","associated_transcript","subcategory")]
reportStats(res=TargetedDESeq$ontResTranAnno$wald8mos$res_Wald, stats=TargetedDESeq$ontResTran$wald8mos$stats_Wald, isoList=ApoeIso)
reportStats(res=TargetedDESeq$ontResTranAnno$wald$res_Wald, stats=TargetedDESeq$ontResTran$wald$stats_Wald, isoList=ApoeIso)
reportStats(res=TargetedDESeq$ontResGeneAnno$wald$res_Wald, stats=TargetedDESeq$ontResGeneAnno$stats_Wald, isoList=c("40586"))

# number of isoforms with 4 known exons
Apoe4Ex <- ApoeExon %>% select(paste0("Gencode_",c(5,6,8)))
table(rowSums(Apoe4Ex  == "No" | Apoe4Ex  == "FirstExon" | Apoe4Ex  == "IR"))


# Differential transcript usage
TargetedDIU$ontDIUGeno$resultDIU %>% filter(Gene == "Apoe")
