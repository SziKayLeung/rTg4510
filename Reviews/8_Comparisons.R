refMap <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged/bulkSorted.rTg4510SCN_collapsed.filtered_targetgenes_filtered.gtf.refmap", header = T)
srefMap <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged/Sortedbulk.all_iso_ont_collapsed.filtered_counts_filtered.gtf.refmap", header = T)
sTMap <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged/Sortedbulk.all_iso_ont_collapsed.filtered_counts_filtered.gtf.tmap", header = T)

length(unique(refMap$ref_id))
length(unique(srefMap$ref_id))
table(sTMap$class_code)


class.files$targ_filtered[class.files$targ_filtered$associated_gene == "Trem2",]

MatchedDiff <- sTMap[sTMap$qry_id %in% TargetedDESeq$ontResTranAnno$wald$anno_res$isoform,]
MatchedDiff <- merge(class.files$targ_sorted[,c("isoform","associated_transcript")], 
      MatchedDiff[,c("ref_id","qry_id")], by.x = "isoform", by.y = "ref_id")
colnames(MatchedDiff) <- c("SCN_isoform","SCN_associated_transcipt","bulk_isoform")

MaxMatchedDiff <- class.files$targ_sorted[class.files$targ_sorted$associated_transcript %in% 
                          MatchedDiff$SCN_associated_transcipt[!grepl("novel", MatchedDiff$SCN_associated_transcipt)],] %>% 
  mutate(AllSum = NeuNSum + DNSum) %>%
  group_by(associated_transcript) %>%
  top_n(1, AllSum) %>% as.data.frame()

MatchedDiff <- merge(MatchedDiff, MaxMatchedDiff[,c("isoform","associated_transcript")], by.x = "SCN_associated_transcipt", by.y = "associated_transcript", all.x = T)

SortedDiffplots <- list()
for(i in 1:nrow(MatchedDiff)){
  isoform <- MatchedDiff$isoform[[i]]
  gene <- class.files$targ_sorted[class.files$targ_sorted$isoform == isoform, "associated_gene"]
  transcript <- class.files$targ_sorted[class.files$targ_sorted$isoform == isoform, "associated_transcript"]
  print(isoform)
  SortedDiffplots[[i]] <- plot_boxplot_SCN(Exp$targ_sorted_all, isoform, ageDiv = FALSE) + 
    theme(legend.title=element_blank()) + labs(title = "", subtitle = paste0(gene,"\n", transcript)) + theme(legend.position = "None")
}


plot_boxplot_SCN(Exp$targ_sorted_all, "PB.90374.201")
plot_boxplot_SCN(Exp$targ_sorted_all, "PB.67362.11")
plot_boxplot_SCN(Exp$targ_sorted_all, "PB.67362.1063")
PB.90374.201
PB.67362.11
PB.67362.1063


### human
library("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

getLDS("ensembl_transcript_id_version", "ensembl_transcript_id", "ENSMUST00000022616", mouse, "ensembl_transcript_id_version", martL = human)
getLDS("ensembl_transcript_id_version", "ensembl_transcript_id", "ENSMUST00000159265", mouse, "ensembl_transcript_id_version", martL = human)
