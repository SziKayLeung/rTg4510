plot_expression_summed <- function(TList, TName = NULL, AgeDiv = FALSE){
  
  dat <- subset(Exp$targ_ont$normAll, row.names(Exp$targ_ont$normAll) %in% TList) %>% dplyr::select(-associated_gene) %>%
    apply(.,2,sum) %>% 
    reshape2::melt(value.name = "sumReads") %>% 
    tibble::rownames_to_column(., var = "sample") %>% 
    mutate(sample = word(sample, c(2), sep = fixed("_"))) %>% 
    left_join(., phenotype$targ_ont, by = "sample") %>%
    mutate(group = factor(ifelse(group == "CASE","TG","WT"), levels = c("WT","TG")))
  
  
  if(isFALSE(AgeDiv)){
    p <- ggplot(dat, aes(x = group, y = sumReads)) + geom_boxplot(outlier.shape = NA) + 
      geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), aes(colour = factor(time)), size = 3) + 
      mytheme + labs(x = "", y = "Normalized counts", subtitle = TName) + 
      scale_colour_manual(values = c("grey","azure4","black","red"), name = "Age (months)") +
      theme(legend.position = "top") 
  }else{
    p <- ggplot(dat, aes(x = group, y = sumReads, colour = as.factor(time))) + geom_boxplot() +
      mytheme + labs(x = "", y = "Normalized counts", subtitle = TName) + 
      scale_colour_manual(values = c("grey","azure4","black","red"), name = "Age (months)") +
      theme(legend.position = "top") 
  }

  
  return(p)
}

plot_expression_summed(CE, "Cryptic exons", AgeDiv = TRUE)

CE <- c(paste0("PB.22007.",c("925","1554","1033","4967","1470","14222","1014")),
        "PB.40586.26305",
        paste0("PB.14646.",c("6992")),
        paste0("PB.38419.",c("3145")),
        paste0("PB.20818.", c("80","192","573","1074","493","362","1096")))
plot_expression_summed(CE)

NMD <- lapply(TargetGene,function(x) class.files$protein_filtered_final[class.files$protein_filtered_final$is_nmd == "TRUE" & class.files$protein_filtered_final$associated_gene == x,"corrected_acc"])
names(NMD) <- TargetGene

NMDp <- lapply(NMD, function(x) plot_expression_summed(x, AgeDiv = TRUE))
for(i in 1:length(TargetGene)){NMDp[[i]] <- NMDp[[i]] + labs(subtitle = TargetGene[[i]])}
NMDGenop <- lapply(NMD, function(x) plot_expression_summed(x, AgeDiv = FALSE))
for(i in 1:length(TargetGene)){NMDGenop[[i]] <- NMDGenop[[i]] + labs(subtitle = TargetGene[[i]])}

names(NMDp) <- TargetGene
plot_grid(plotlist = NMDp)
plot_grid(plotlist = NMDGenop)


IRList <- lapply(TargetGene,function(x) IR[IR$associated_gene == x,"transcript_id"])
names(IRList) <- TargetGene
IRGenop <- lapply(IRList, function(x) plot_expression_summed(x, AgeDiv = FALSE))
for(i in 1:length(TargetGene)){IRGenop[[i]] <- IRGenop[[i]] + labs(subtitle = TargetGene[[i]])}
plot_grid(plotlist = IRGenop)
