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

#SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Paper_Figures/"
#source(paste0(SC_ROOT, "rTg4510_config.R"))
#source(paste0(SC_ROOT, "0_source_functions.R"))
#source(paste0(SC_ROOT,"bin/draw_heatmap_gene_level.R"))


## ---------- Global Iso-Seq dataset -----------------

# mean number of genes
cat("Mean number of genes:", mean(sapply(sub_class.files, function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("Mean number of genes in WT samples:", mean(sapply(sub_class.files[wholeWT], function(x) length(unique(x[["associated_gene"]])))),"\n")
cat("Mean number of genes in TG samples:", mean(sapply(sub_class.files[wholeTG], function(x) length(unique(x[["associated_gene"]])))),"\n")


# mean number of isoforms
cat("Mean number of isoforms:", mean(sapply(sub_class.files, function(x) nrow(x))),"\n")
cat("Mean number of isoforms in WT samples:", mean(sapply(sub_class.files[wholeWT], function(x) nrow(x))),"\n")
cat("Mean number of isoforms in TG samples:", mean(sapply(sub_class.files[wholeTG], function(x) nrow(x))),"\n")


# mean descriptives
cat("Mean isoform length in all samples:", mean(class.files$glob_iso$length), "\n")
cat("Mean number of exons in all samples:", mean(class.files$glob_iso$exons), "\n")
cat("Mean number of exons in all samples:", mean(class.files$glob_iso$exons), "\n")
cat("Mean number of exons in all samples:", mean(as.data.frame(class.files$glob_iso %>% group_by(associated_gene) %>% tally())$n), "\n")


# Global data 
for(i in c("Gfap","C4b")){
  cat("***************** Gene level\n")
  cat("RNA-Seq gene level statistics for", i, "\n")
  print(GlobalDESeq$RresGeneAnno$lrt$anno_res %>% filter(associated_gene == i))
  print(GlobalDESeq$RresGeneAnno$lrt$stats_LRT %>% filter(associated_gene == i) %>% select(associated_gene, LRTStatistic))
  cat("Iso-Seq gene level statistics for", i, "\n")
  print(GlobalDESeq$resGeneAnno$lrt$anno_res %>% filter(associated_gene == i))
  print(GlobalDESeq$resGeneAnno$lrt$stats_LRT %>% filter(associated_gene == i) %>% select(associated_gene, LRTStatistic))
}

for(i in c("PB.2973.16","PB.7022.9")){
  cat("***************** Transcript level\n")
  cat("RNA-Seq gene level statistics for", i, "\n")
  print(GlobalDESeq$RresTranAnno$lrt$anno_res %>% filter(isoform == i))
  print(GlobalDESeq$RresTranAnno$lrt$stats_LRT %>% filter(rownames(.) == i) %>% select(LRTPvalue, LRTStatistic))
  cat("Iso-Seq gene level statistics for", i, "\n")
  print(GlobalDESeq$resTranAnno$lrt$anno_res %>% filter(isoform == i))
  print(GlobalDESeq$resTranAnno$lrt$stats_LRT %>% filter(rownames(.) == i) %>% select(LRTPvalue, LRTStatistic))
}


## ---------- Merged Targeted Iso-Seq and ONT datasets -----------------

class.files$targ_filtered %>% filter(dataset == "ONT") %>% nrow()
class.files$targ_all %>% filter(dataset == "Both") %>% nrow()
class.files$targ_all %>% filter(dataset == "Iso-Seq") %>% nrow()

nrow(class.files$targ_all)
nrow(class.files$targ_all %>% filter(associated_transcript != "novel"))
nrow(class.files$targ_all %>% filter(associated_transcript == "novel"))

nrow(class.files$targ_filtered)
nrow(class.files$targ_filtered %>% filter(associated_transcript != "novel"))
nrow(class.files$targ_filtered %>% filter(associated_transcript == "novel"))

class.files$targ_filtered %>% filter(dataset == "ONT") %>% nrow()
class.files$targ_filtered %>% filter(dataset == "Iso-Seq") %>% nrow()

nrow(class.files$targ_all %>% filter(nreads < 10))

# summary
all_summarise_gene_stats(Targeted$Gene_class, class.files$targ_filtered, Targeted$cpat, Targeted$noORF, Targeted$Genes)

# mean number of isoforms 
mean(as.data.frame(class.files$targ_filtered %>% group_by(associated_gene) %>% tally())$n)
class.files$targ_filtered %>% group_by(associated_gene) %>% tally() %>% arrange(n)
class.files$targ_filtered %>% mutate(noveltrans = ifelse(associated_transcript == "novel","novel","known")) %>% 
  group_by(associated_gene, noveltrans) %>% tally() %>% arrange(n)


# Targeted stats 
cat("Number of ES events across 20 target genes:", sum(Merged_gene_class_df["ES",]), "\n")
cat("Number of trancripts with ES events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.ES.Events",]), "\n")
cat("Proportion of trancripts with ES events:", sum(Merged_gene_class_df["Number.of.Transcripts.with.ES.Events",]), "\n")
sum(Merged_gene_class_df["Total.Number.of.Transcripts",])






