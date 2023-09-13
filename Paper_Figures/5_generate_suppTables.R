
formulateStats <- function(dat, gene){
  
  df <- dat[["norm_counts"]] %>%
    filter(isoform %in%  dat[["anno_res"]][["isoform"]]) %>%
    group_by(isoform, associated_gene, time, group) %>%
    dplyr::summarize(mean_norm_counts = mean(normalised_counts)) %>% 
    mutate(group = ifelse(group == "CONTROL","WT","TG"), isoform = as.character(isoform)) %>%
    mutate(groups = paste0(group," ",time,"mos Counts")) %>% as.data.frame() %>%
    dplyr::select(isoform, groups, mean_norm_counts)  %>%
    dcast(., isoform  ~ groups, value.var = "mean_norm_counts") %>%
    full_join(.,dat[["anno_res"]], by = "isoform") %>% 
    mutate_if(is.numeric, ~ signif(., digits = 3)) %>%
    dplyr::rename("FDR" = "padj", "P Value" = "pvalue", "Associated Gene" = "associated_gene")
  
  if(missing(gene)){
    df <- df %>%
      dplyr::select(`Associated Gene`, isoform, associated_transcript, log2FoldChange, lfcSE, `P Value`, FDR, ends_with("Counts")) %>%
      dplyr::rename("Associated Transcript" = "associated_transcript", "Isoform" = "isoform")
  }else{
    df <- df %>% dplyr::select(`Associated Gene`, log2FoldChange, lfcSE, `P Value`, FDR, ends_with("Counts")) 
  }
  
  return(df)
}



### ------------- Supplementary Table 3 -------------------
# DESeq2 gene level expression analyses output for whole transcriptome data

TSupp3 <- formulateStats(GlobalDESeq$resGeneAnno$wald,"gene")

### ------------- Supplementary Table 4 -------------------
# DESeq2 transcript level expression analyses output for whole transcriptome data

TSupp4 <- formulateStats(GlobalDESeq$resTranAnno$lrt)
TSupp4B <- formulateStats(GlobalDESeq$resTranAnno$wald)


### ------------- Supplementary Table 5 -------------------
# DESeq2 transcript level expression analyses output for targeted transcriptome data

# Additional transcripts that are significant at 8 months in ONT dataset
TargetedDESeq$ontResTranAnno$wald8mosUnique <- list(
  norm_counts = TargetedDESeq$ontResTranAnno$wald8mos$norm_counts,
  anno_res = TargetedDESeq$ontResTranAnno$wald8mos$anno_res %>% filter(!isoform %in% TargetedDESeq$ontResTranAnno$wald$anno_res$isoform)
)


TSupp5 <- dplyr::bind_rows(
  formulateStats(TargetedDESeq$ontResTranAnno$wald) %>% dplyr::mutate(Dataset = "ONT", Model = "WaldAllAges"),
  formulateStats(TargetedDESeq$isoResTranAnno$wald) %>% dplyr::mutate(Dataset = "Iso-Seq", Model = "WaldAllAges"),
  formulateStats(TargetedDESeq$ontResTranAnno$wald8mosUnique) %>% dplyr::mutate(Dataset = "ONT", Model = "Wald8mos")) %>%
  dplyr::select(Dataset, Model, everything())

TargetedMergedDESeq$waldGenotypeAgeInteraction <- merge(TargetedMergedDESeq$waldGenotypeAgeInteraction,class.files$targ_all[,c("isoform","subcategory")], by = "isoform", all.x=T)
TargetedMergedDESeq$waldGenotype <- merge(TargetedMergedDESeq$waldGenotype,class.files$targ_all[,c("isoform","subcategory")], by = "isoform", all.x=T)

write.table(TargetedMergedDESeq$waldGenotypeAgeInteraction, paste0(dirnames$targ_output,"/TargetedDESeq_ontResTranAnno_resWaldGenotypeAge.txt"), quote = F, row.names = F, sep = "\t")
write.table(TargetedMergedDESeq$waldGenotype, paste0(dirnames$targ_output,"/TargetedDESeq_ontResTranAnno_resWaldGenotype.txt"), quote = F, row.names = F, sep = "\t")

top10Genotype <- TargetedMergedDESeq$waldGenotype %>% arrange(padj_ont) %>% .[1:10,"isoform"]
bottom10Genotype <- TargetedMergedDESeq$waldGenotype %>% filter(padj_ont < 0.05) %>% arrange(-padj_ont) %>% .[1:10,"isoform"]

pdf(paste0(output_dir,"/GenotypeTop10Bottom10.pdf"), width = 10, height = 8)
for(i in c(top10Genotype, bottom10Genotype)){
  g = class.files$targ_all[class.files$targ_all$isoform == i,"associated_gene"]
  print(plot_transexp_overtime(g,TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=2,isoSpecific=c(i),setorder=c("CONTROL","CASE")))
}
dev.off()

nrow(TargetedMergedDESeqSig$waldGenotypeAge)
nrow(TargetedMergedDESeqSig$waldGenotypeAgeInteraction) 
nrow(TargetedMergedDESeqSig$waldGenotype)
nrow(TargetedMergedDESeqSig$lrt)
InteractionOnly = setdiff(TargetedMergedDESeqSig$waldGenotypeAgeInteraction$isoform,TargetedMergedDESeqSig$waldGenotype$isoform)
GenotypeAlone = setdiff(TargetedMergedDESeqSig$waldGenotype$isoform,TargetedMergedDESeqSig$waldGenotypeAgeInteraction$isoform)

pdf(paste0(output_dir,"/GenotypeOnly.pdf"), width = 10, height = 8)
for(i in GenotypeAlone){
  g = class.files$targ_all[class.files$targ_all$isoform == i,"associated_gene"]
  print(plot_transexp_overtime(g,TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=2,isoSpecific=c(i),setorder=c("CONTROL","CASE")))
}
dev.off()

 

pdf(paste0(output_dir,"/InteractionOnly.pdf"), width = 10, height = 8)
for(i in InteractionOnly){
  g = class.files$targ_all[class.files$targ_all$isoform == i,"associated_gene"]
  print(plot_transexp_overtime(g,TargetedDESeq$ontResTranAnno$wald$norm_counts,show="specific",rank=2,isoSpecific=c(i),setorder=c("CONTROL","CASE")))
}
dev.off()


### ------------- Supplementary Table 6 -------------------
# DIU analyses output for targeted transcriptome data

TSupp6 <- dplyr::bind_rows(
  TargetedDIU$ontDIUGeno$resultDIU %>% dplyr::mutate(Dataset = "ONT", Model = "Genotype"),
  TargetedDIU$ontDIUPath$resultDIU %>% dplyr::mutate(Dataset = "ONT", Model = "Pathology")) %>%
  dplyr::select(Dataset, Model, everything()) %>%
  dplyr::rename("P Value" = "p.value") %>%
  mutate_if(is.numeric, ~ signif(., digits = 3)) 

