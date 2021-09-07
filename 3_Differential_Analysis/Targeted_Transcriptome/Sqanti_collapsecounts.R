library("dplyr")
library("stringr")
library(tidyverse)
library("ggplot2")

output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_INPUT/IsoSeq_Expression"

# find the isoforms that are differentially expressed
targeted.class.files = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt", sep = "\t",as.is = T, header = T)

sqanti = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI2/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt", header = T)

tama_discard = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAMA/AllMouseTargeted_discarded.txt") %>% 
  mutate(isoform = word(V4,c(2), sep = ";"))

# removed transcript structural category
sqanti %>% filter(isoform %in% tama_discard$isoform) %>% group_by(structural_category) %>% tally() %>% ggplot(., aes(x = structural_category, y = n)) + geom_bar(stat = "identity")

TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")
allcounts_known = sqanti %>% filter(associated_transcript != "novel") %>% filter(associated_gene %in% TargetGene) %>% select(associated_transcript,starts_with("FL.")) %>% 
  group_by(associated_transcript) %>%
  summarise_each(funs(sum))

# allcounts_known[allcounts_known$associated_transcript == "ENSMUST00000114268.4","FL.K17"] == sum(sqanti[sqanti$associated_transcript == "ENSMUST00000114268.4","FL.K17"])

retained_known = targeted.class.files[targeted.class.files$associated_transcript %in% allcounts_known$associated_transcript,c("isoform","structural_category","associated_transcript")]
retained_known_annotate = merge(retained_known,allcounts_known, all.x = T)

# determine the number of transcripts named more than once in sqanti due to incomplete splice match
morethanone = as.data.frame(retained_known_annotate %>% group_by(associated_transcript) %>% tally() %>% .[.$n > 1,"associated_transcript"])
oneFSM = list()
moreFSM = list()
noFSM = list()
for(trans in morethanone$associated_transcript){
  cate = retained_known_annotate[retained_known_annotate$associated_transcript == trans,"structural_category"]
  if('full-splice_match' %in% cate){
    if(length(cate[cate == "full-splice_match"]) == 1){
     oneFSM[[trans]] = trans 
    }else{
     moreFSM[[trans]] = trans  
    }
  }else{
    noFSM[[trans]] = trans
  }
}

# oneFSM = keep transcripts with only one FSM
# noFSM = keep ISM transcripts as no FSM 
# moreFSM = empty()
retained_known_annotate_oneFSM = targeted.class.files[targeted.class.files$associated_transcript %in% row.names(do.call("rbind", oneFSM)) & targeted.class.files$structural_category == "full-splice_match",]

targeted.class.files[targeted.class.files$associated_transcript %in% row.names(do.call("rbind", noFSM)),]

#nrow(retained_known_annotate_oneFSM) == length(oneFSM)

# 
#targeted.class.files[targeted.class.files$associated_transcript == "ENSMUST00000220858.1",]
#targeted.class.files[targeted.class.files$associated_gene == "Vipas39",]
#sqanti[sqanti$associated_transcript == "ENSMUST00000220858.1",]
#sqanti[sqanti$associated_gene == "Vipas39",]
#tama_discard[tama_discard$isoform == "PB.1998.1",]

allcounts_novel = targeted.class.files[targeted.class.files$associated_transcript == "novel",] %>% select(isoform,starts_with("FL."))


allcounts = rbind(retained_known_annotate %>% select(-associated_transcript),allcounts_novel)
#nrow(allcounts) == nrow(targeted.class.files)

# only retain AD transcripts 
AD_transcripts = targeted.class.files[toupper(targeted.class.files$associated_gene) %in% TargetGene,"isoform"]

allcounts_AD = allcounts %>% filter(isoform %in% AD_transcripts)
nrow(allcounts_AD) == length(AD_transcripts)

colnames(allcounts_AD)[1] = ""
write.table(allcounts_AD,paste0(output,"/AllMouseTargeted_sqantitamafiltered_collapsed.expression.txt"), quote = F, row.names = F, sep = "\t")
