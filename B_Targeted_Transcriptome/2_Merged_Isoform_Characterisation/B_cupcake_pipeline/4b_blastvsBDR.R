library("dplyr")

blastOutput <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/4_characterise/BDRBlast/allIsoOntBDRHigh.txt")
blastOutput <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/4_characterise/BDRBlast/allIsoOntBDRLow.txt")
colnames(blastOutput) <- c("human","mouse","Perc_Identity","alignment_length","mismatches", "gap","q.start","q.end","s.start","s.end","evalue","bit_score")

classfile <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/2_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt", sep = "\t", as.is = T)

humanclassfile <- read.table("/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/E_MergedTargeted/3_sqanti3/all_iso_ont_scn_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt", sep = "\t", as.is = T)

## Filtering requirements
# > 200bp length 
# > 90% blast identity 

# filter blast hits > 200bp 
filtered <- blastOutput[blastOutput$alignment_length > 200 & blastOutput$Perc_Identity > 90, ]

topranked <- blastOutput %>% group_by(human) %>% filter(Perc_Identity == max(Perc_Identity) & bit_score == max(bit_score))

topranked <- merge(topranked, classfile[,c("isoform","structural_category","associated_gene","associated_transcript")], by.x = "mouse",by.y = "isoform", all.x = T)
filtered[filtered$mouse == "PB.20818.54",]

length(unique(filtered$human))



