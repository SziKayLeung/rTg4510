trans = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/G_Merged_Targeted/B_cupcake_pipeline/4_characterise/Transdecoder/pfam.domtblout.read")
colnames(trans) = c("targetName","targetAccession","tlen","queryName","accession","querylen","E-value","score","bias","domainNum","numDomains","cEvalue","iEvalue","scoreDomain",
                    "scoreBias","fromHmm","toHmm","fromAli","toAli","fromEnv","toEnv","acc","description","len","type","orf")

trans = trans %>% mutate(isoform = word(targetName,c(1),sep=fixed("_")))

transAnno = merge(class.files$targ_filtered[,c("isoform","associated_transcript","associated_gene","structural_category")], trans, by = "isoform")                 
transAnno[transAnno$isoform == "PB.20818.54",]
transAnno[transAnno$isoform == "PB.20818.1074",]

Trem2TransAnno <- transAnno[transAnno$associated_gene == "Trem2",c("isoform","associated_transcript","associated_gene","structural_category","targetName","queryName","iEvalue")] 

#Replace missing values with Inf
Trem2TransAnno[is.na(Trem2TransAnno)] <- Inf

# Group the Trem2TransAnno frame by isoform and queryName
groups <- split(Trem2TransAnno, list(Trem2TransAnno$isoform, Trem2TransAnno$queryName))

# Select the row with the minimum non-Inf iEvalue for each group
result <- do.call(rbind, lapply(groups, function(group) {
  min_iEvalue <- min(group$iEvalue[group$iEvalue < Inf], na.rm = TRUE)
  group[group$iEvalue == min_iEvalue, ]
}))

result <- result %>% tidyr::spread(queryName, iEvalue)

# note the query value is very high for some of the domains

loopDomains <- function(result){
  cols <- as.data.frame(t(result)) %>% select(-c("isoform","associated_transcript","targetName","associated_gene","structural_category"))
  result_list <- list()
  for(i in cols){
    if(is.na(i)){
      result_list <- append(result_list, 0)
    }else{
      result_list <- append(result_list, 1)
    }
  }
  return(paste0(unlist(result_list), collapse = ""))
}
result$combo <- apply(result,1,loopDomains)

result %>% group_by(structural_category, combo) %>% tally() %>% 
  ggplot(., aes(x = combo, y = n, fill = structural_category)) + geom_bar(stat = "identity") +
  labs(x = "Combination of domains", y = "Number of transcripts") +
  theme_classic()
