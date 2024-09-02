dirnames$protein <- paste0("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/B_cupcake_pipeline/4_proteogenomics/")
protein = list(
  cpat = read.table(paste0(dirnames$protein,"5_calledOrfs/all_iso_ont_best_orf.tsv"), sep ="\t", header = T),
  t2p.collapse = read.table(paste0(dirnames$protein,"6_refined_database/all_iso_ont_orf_refined.tsv"), sep = "\t", header = T)
)

dirnames$protein <- paste0(dirnames$targ_root,"/4_proteogenomics/7_classified_protein/")

class.files$protein <- read.table(paste0(dirnames$protein,"all_iso_ont.sqanti_protein_classification.tsv"), sep = "\t", header = T)
nmd <- read.table(paste0(dirnames$protein, "all_iso_ont.classification_filtered.tsv"), sep = "\t", header = T)
idx <- match(nmd$pb,class.files$ptarg_filtered$base_acc)
nmd <- transform(nmd, corrected_acc = ifelse(!is.na(idx), as.character(class.files$ptarg_filtered$corrected_acc[idx]), NA))

## ---------- cryptic exons -----------------

InternalNovelExons <- lapply(FICLENE, function(x) x[x$classification == "Internal_NovelExon",])
InternalNovelExons <- lapply(InternalNovelExons, function(x) merge(x, 
                                                                   protein$cpat[,c("pb_acc","coding_score","orf_calling_confidence")], 
                                                                   by.x = "transcriptID", by.y = "pb_acc"))
InternalNovelExons <- lapply(InternalNovelExons, function(x) x %>% mutate(coding = ifelse(coding_score < 0.44, "noncoding","coding")))
InternalNovelExons <- lapply(InternalNovelExons, function(x) merge(x, 
                                                                   nmd[,c("pb","is_nmd","has_stop_codon")], 
                                                                   by.x = "transcriptID", by.y = "pb",all.x=TRUE)) 

plotCE <- function(gene, type){
  representative <- as.data.frame(gtf$ref_target %>% filter(type == "transcript" & transcript_type == "protein_coding") %>% group_by(gene_name) %>% top_n(1, width))
  allRep <- as.data.frame(gtf$ref_target)
  # manual drop gene
  print(gene)
  df <- InternalNovelExons[[gene]]
  if(gene == "Apoe"){
    # PB.40586.2026 = reference match ("PB.40586.2026")
    df <- df[df$transcriptID %in% c("PB.40586.26305"),]
    ref <-data.frame("PB.40586.2026",NA,NA,NA,NA,NA,"coding","False","True")
    colnames(ref) <- colnames(df)
    df <- rbind(df, ref)
  }else if(gene == "Bin1"){
    df <- df[df$transcriptID %in% paste0("PB.22007.",c("925","1554","1033","4967","1470","14222","1014")),]
  }else if(gene == "Clu"){
    df <-  df[df$transcriptID %in% paste0("PB.14646.",c("6992")),]
    # PB.40586.2026 = reference match ("PB.14646.139", "PB.14646.483")
    ref <- data.frame("PB.14646.139",NA,NA,NA,NA,NA,"coding","False","True")
    colnames(ref) <- colnames(df)
    df <- rbind(df, ref)
  }else if(gene == "Snca"){
    df <-  df[df$transcriptID %in% paste0("PB.38419.",c("3145")),]
  }else if(gene == "Trem2"){
    df <- df[df$transcript %in% paste0("PB.20818.", c("80","192","573","1074","493","362","1096")),]
  }else{
    return(NULL)
  }
  transcriptList <- list(
    #Reference = c("ENSMUST00000024791.14","ENSMUST00000113237.3"),
    Reference =  unique(allRep[allRep$gene_name == gene, "transcript_id"]),
    #Reference =  representative[ representative$gene_name == gene, "transcript_id"],
    `not NMD` = as.character(unique(df %>% filter(coding == "coding" & is_nmd == "False") %>% .[,c("transcriptID")])),
    `NMD` = as.character(unique(df %>% filter(coding == "coding" & is_nmd == "True") %>% .[,c("transcriptID")])),
    `non-coding` = as.character(unique(df %>% filter(coding == "noncoding") %>% .[,c("transcriptID")]))
  )
  proteinList <- list(
    #Reference =  unique(allRep[allRep$gene_name == gene, "transcript_id"]),
    `not NMD` =  unique(gtf$ptarg_merged %>% filter(transcript %in% transcriptList$`not NMD`) %>% .[["gene_id"]]),
    `NMD` = unique(gtf$ptarg_merged %>% filter(transcript %in% transcriptList$`NMD`) %>% .[["gene_id"]]),
    `non-coding` = unique(gtf$ptarg_merged %>% filter(transcript %in% transcriptList$`non-coding`) %>% .[["gene_id"]]))
  
  if(type == "transcript"){
    p <-  generalggTranPlots(transcriptList, gtf$targ_merged, class.files$targ_filtered, gene) 
  }else if(type == "protein"){
    p <- generalggTranPlots(proteinList, gtf$targ_merged, class.files$targ_filtered, gene) 
  }else{
    bothList <- c(transcriptList, proteinList)
    p <- generalggTranPlots(bothList, gtf$targ_merged, class.files$targ_filtered, gene) 
  }
  
  return(p)
}
