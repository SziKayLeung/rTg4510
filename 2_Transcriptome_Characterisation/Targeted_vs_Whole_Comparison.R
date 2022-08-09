suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("wesanderson"))
suppressMessages(library("tidyr"))

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")

# plot theme
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=18,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 18,family="CM Roman"),
                 axis.text.x= element_text(size=16,family="CM Roman"),
                 axis.text.y= element_text(size=16,family="CM Roman"))


### Read files #################################
# output file from TAMA merge 
TAMA_transfile <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Whole_vs_Targeted/merged_whole_targeted_trans_report.txt", header = T) %>% mutate(gene_id = word(transcript_id, c(1), sep = fixed("."))) %>% mutate(gene_name = word(word(all_source_trans, c(3), sep = fixed("_")),c(1), sep = ","))

# output gene report from TAMA merge
genereport <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Whole_vs_Targeted/merged_whole_targeted_gene_report.txt", header = T)

# gtf from TAMA
tama_gtf <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Whole_vs_Targeted/merged_whole_targeted.gtf", as.is = T, sep = "\t")

# sqanti2 classification file from whole transcriptome 
whole.class.names.files <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.txt"
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

# sqanti2 classification file from targeted transcriptome 
targeted.class.names.files <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt"

### File preparation #################################
# Read sqanti files
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")
targeted.class.files <- SQANTI_class_preparation(targeted.class.names.files,"standard")

# subset isoforms annotated to target genes
TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")

wholeAD_Pbid <- whole.class.files[toupper(whole.class.files$associated_gene) %in% TargetGene,"isoform"]
targetedAD_Pbid <- targeted.class.files[toupper(targeted.class.files$associated_gene) %in% TargetGene,"isoform"]

# Number of transcripts from tama merge by source (independent of target genes)
trans %>% group_by(sources) %>% tally()

# Subset from Tama merge file only transcripts from AD target genes 
# https://stackoverflow.com/questions/42646626/searching-for-a-list-of-string-in-a-dataframe-in-r
wholeADtrans <- trans[unlist(lapply(wholeAD_Pbid,function(x) grep(x, trans$all_source_trans, fixed = TRUE))),]
targetedADtrans <- trans[unlist(lapply(targetedAD_Pbid,function(x) grep(x, trans$all_source_trans, fixed = TRUE))),]

setdiff(wholeADtrans$transcript_id,targetedADtrans$transcript_id)
setdiff(targetedADtrans$transcript_id,wholeADtrans$transcript_id)

# Common transcripts in both approach 
intersect(wholeADtrans$transcript_id,targetedADtrans$transcript_id)

### 

whole_vs_targeted_plots <- function(){
  # from TAMA merge transcript file, filter transcripts associated with target gene 
  ADtrans <- TAMA_transfile %>% filter(toupper(gene_name) %in% TargetGene)
  # from TAMA merge transcript file, filter transcripts not associated with target gene 
  nonADtrans <- TAMA_transfile %>% filter(!toupper(gene_name) %in% TargetGene)
  
  #### plots ###
  #p1 = number of isoforms per target gene identified in whole transcriptome, targeted transcriptome, and in both
  #p2 = FL read counts in targeted transcriptome of isoforms annotated to AD genes (target) that are detected uniquely in targeted transcriptome and also in whole
  #p3 = FL read counts in whole transcriptome of isoforms not annotated to AD genes (off target) that are detected uniquely in whole and also in targeted
  
  ##### p1
  # calculate the total number of merged transcripts for plot in reorder 
  ADtrans_ttgene <- TAMA_transfile %>% group_by(gene_name) %>% tally()
  
  # Tally the merged transcripts by where it is from and also by gene name 
  # source from tama merge transcript file would either be "Targeted" or "Targeted,Whole" or "Whole"
  p1 <- ADtrans %>% group_by(sources,gene_name) %>% tally() %>%
    mutate(sources = recode(sources,"Targeted,Whole" = "Both")) %>%
    mutate(sources = factor(sources, levels = c("Whole","Both","Targeted"))) %>%
    left_join(.,ADtrans_ttgene, by = "gene_name") %>%
    ggplot(.,aes(x = reorder(gene_name,-n.y), y = n.x, fill = sources)) + geom_bar(stat = "identity") + mytheme + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "Number of Isoforms") +
    scale_fill_manual(name = "", values = c(label_colour("whole"),label_colour("whole+targeted"),label_colour("targeted")))
  
  ##### p2, p3
  # iso_filter
  # filter the transcripts in whole or targeted transcriptome to determine the difference in expression 
  # if type == ADGenes, then subset the transcripts that are in tama merge file associated with ADGenes and are either in targeted, or in targeted and equivalent whole, then find the targeted PacBio Id of that transcript, and find the FL; 
  # if type == NonADGenes, then subset the transcripts that are in tama merge file not associated with ADGenes and are either in whole, or in targeted and equivalent whole, then find the Whole PacBio Id of that transcript, and find the FL; 
  iso_filter <- function(input_dat, input_class,type){
    if(type == "ADGenes"){
      fil_group <- c("Targeted","Targeted,Whole")
      name <- "Targeted"
      colour <- scale_fill_manual(name = "", values = c(label_colour("targeted"),label_colour("whole+targeted")))
    }else if(type == "NonADGenes"){
      fil_group <- c("Whole","Targeted,Whole")
      name <- "Whole"
      colour <-  scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("whole")))
    }else{
      print("1")
    }
    
    dat <- input_dat %>% filter(sources %in% fil_group) %>% 
      # separate column by the different sources if detected in both whole and targeted
      # the order from TAMA is not consistent i.e one row would be wholeXXX,targetedXXX and another would be targetedXXX,wholeXXX
      # therefore use ifelse to grep the correct PBID reference to whole and targeted transcriptome
      separate(all_source_trans, c("match1", "match2"), ",") %>% 
      mutate(Whole = word(ifelse(word(match1, c(1), sep = fixed("_")) == "Whole",
                                 word(match1, c(2), sep = fixed("_")),
                                 word(match2, c(2), sep = fixed("_"))))) %>% 
      mutate(Targeted = word(ifelse(word(match1, c(1), sep = fixed("_")) == "Targeted",
                                    word(match1, c(2), sep = fixed("_")),
                                    word(match2, c(2), sep = fixed("_"))))) 
    
    # merge the pacbio ID from the whole or targeted transcriptome to the classification file to extract the FL reads
    dat <- merge(dat[,c("sources",name)],input_class[,c("isoform","FL")], by.x = name, by.y = "isoform", all.x = T) %>% 
      mutate(sources = recode(sources, "Targeted,Whole" = "Whole & Targeted", name = paste0(name," Only")))
    
    p <- ggplot(dat, aes(x = sources, y = log10(FL), fill = sources)) + geom_boxplot() + mytheme + 
      labs(x = "Transcriptome Approach", y = paste0("FL Read Counts \n",name," Transcriptome (Log10)")) +
      colour + theme(legend.position = "none")
    
    return(p)
  }
  
  p2 <- iso_filter(ADtrans,targeted.class.files,"ADGenes")
  p3 <- iso_filter(nonADtrans,whole.class.files,"NonADGenes")
  
  return(list(p1,p2,p3))
}

