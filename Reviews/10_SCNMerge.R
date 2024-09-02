read_stat <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged/merged_collapse.read_stat.txt", header = T)
sampleID <- read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged/original_fasta/merged_sample_id.csv")

read_stat <- merge(read_stat, sampleID,by = "id")

read_stat <- read_stat %>% mutate(dataset = ifelse(grepl("SCN",primer), "sorted","bulk"))
head(read_stat)
spread(read_stat, key = "pbid", value = "dataset")
bulk <- read_stat %>% filter(dataset == "bulk")
sorted <- read_stat %>% filter(dataset == "sorted")

mergeAll <- merge(sorted, bulk, by = "pbid", all = T) %>% dplyr::select(pbid, id.x, id.y)
colnames(mergeAll) <- c("newcollapsed_isoform","sorted_isoform","bulk_isoform")

identify_dataset_by_na <- function(col1,col2,name1,name2){

  
  if(!is.na(col1) & !is.na(col2)){return("Both")
  }else if(is.na(col1) & !is.na(col2)){return(name2)
  }else if(!is.na(col1) & is.na(col2)){return(name1)
  }else{return("NA")}
  
}

mergeAll$dataset <- apply(mergeAll, 1, function(x) identify_dataset_by_na (x[["sorted_isoform"]], x[["bulk_isoform"]], "Sorted","Bulk"))

mergeAll <- merge(class.files$targ_filtered[,c("isoform","associated_gene","associated_transcript","structural_category")], mergeAll, by.x = "isoform", by.y = "bulk_isoform", all = T) %>% mutate(sorted_reIsoform = ifelse(is.na(sorted_isoform),"PB", isoform))
vennMerged <- twovenndiagrams(mergeAll[mergeAll$dataset %in% c("Bulk","Both"),"isoform"], mergeAll[mergeAll$dataset %in% c("Sorted","Both"),"sorted_reIsoform"],"Bulk","Sorted")
plot_grid(vennMerged)

merge(class.files$targ_filtered, mergeAll[,c("isoform","dataset")], all.x = T, by = "isoform") %>% 
  group_by(associated_gene, structural_category, dataset) %>% 
  tally() %>% 
  ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + geom_bar(stat = "identity") + facet_grid(~dataset) 

class.files$targ_filtered %>% filter(isoform %in% mergeAll[mergeAll$dataset == "Bulk", ])


## ----- 12 DTE
newCollasedID <- read_stat %>% filter(dataset == "bulk", id %in% TargetedDESeq$ontResTranAnno$wald$anno_res$isoform) %>% .[,c("pbid")]

read_stat <- merge(read_stat %>% filter(pbid %in% newCollasedID), 
      TargetedDESeq$ontResTranAnno$wald$anno_res[,c("isoform","associated_transcript","associated_gene","pvalue")], by.x = "id", by.y = "isoform", all.x = T)

merged <- merge(read_stat %>% filter(dataset == "bulk"), read_stat %>% filter(dataset == "sorted") %>% dplyr::select(id, pbid), by = "pbid", all = T)
colnames(merged) <- c("newcollapsed_isoform","bulk_isform","primer","dataset","associated_transcript","associated_gene","pvalue","sorted_isoform")
merged <- merged %>% dplyr::select(-dataset, -primer) %>% mutate(Detection = ifelse(is.na(sorted_isoform),"Bulk","Both"))
mergedBoth <- merged %>% filter(Detection == "Both")

SortedDiffplots <- list()
for(i in 1:nrow(mergedBoth)){
  isoform <- mergedBoth$sorted_isoform[[i]]
  gene <- class.files$targ_sorted[class.files$targ_sorted$isoform == isoform, "associated_gene"]
  transcript <- class.files$targ_sorted[class.files$targ_sorted$isoform == isoform, "associated_transcript"]
  print(isoform)
  SortedDiffplots[[i]] <- plot_boxplot_SCN(Exp$targ_sorted_all, isoform, ageDiv = FALSE, genotypeDiv = FALSE) + 
    theme(legend.title=element_blank()) + labs(title = "", subtitle = paste0(gene,"\n", transcript)) + theme(legend.position = "None")
}
plot_grid(plotlist = SortedDiffplots)



## ----- > 30 DTE by genotype 

newCollasedID <- read_stat %>% filter(dataset == "bulk", id %in% TargetedDESeq$ontResTranAnno$waldgenotype$anno_res$isoform[1:53]) %>% .[,c("pbid")]

read_stat <- merge(read_stat %>% filter(pbid %in% newCollasedID), 
                   TargetedDESeq$ontResTranAnno$wald$anno_res[,c("isoform","associated_transcript","associated_gene","pvalue")], by.x = "id", by.y = "isoform", all.x = T)

merged <- merge(read_stat %>% filter(dataset == "bulk"), read_stat %>% filter(dataset == "sorted") %>% dplyr::select(id, pbid), by = "pbid", all = T)
colnames(merged) <- c("newcollapsed_isoform","bulk_isform","primer","dataset","associated_transcript","associated_gene","pvalue","sorted_isoform")
merged <- merged %>% dplyr::select(-dataset, -primer) %>% mutate(Detection = ifelse(is.na(sorted_isoform),"Bulk","Both"))
mergedBoth <- merged %>% filter(Detection == "Both")
nrow(mergedBoth)

BSortedDiffplots <- list()
for(i in 1:nrow(mergedBoth)){
  isoform <- mergedBoth$sorted_isoform[[i]]
  gene <- class.files$targ_sorted[class.files$targ_sorted$isoform == isoform, "associated_gene"]
  transcript <- class.files$targ_sorted[class.files$targ_sorted$isoform == isoform, "associated_transcript"]
  print(isoform)
  BSortedDiffplots[[i]] <- plot_boxplot_SCN(Exp$targ_sorted_all, isoform, ageDiv = FALSE) + 
    theme(legend.title=element_blank()) + labs(title = "", subtitle = paste0(gene,"\n", transcript)) + theme(legend.position = "None")
}
plot_grid(plotlist = BSortedDiffplots)



