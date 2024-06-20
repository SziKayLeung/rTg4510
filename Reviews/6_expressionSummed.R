



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
