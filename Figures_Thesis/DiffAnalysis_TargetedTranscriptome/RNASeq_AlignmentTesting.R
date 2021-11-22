library("dplyr")
library("ggplot2")
library("cowplot")
library("stringr")

root_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/RNASeq_SQANTI3/TESTING/"
diff_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/"
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7","Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")
output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_TargetedTranscriptome"

files = list(paste0(root_dir,"K24_AllReadsCollapsed_FSMLength/abundance.tsv"),
             #paste0(root_dir,"K24_AllReadsCollapsed_FSMExp_BTS/abundance.tsv"),
             paste0(root_dir,"K24_AllReadsCollapsed_FSMExp/abundance.tsv"),
             paste0(root_dir,"K24_TargetedReadsCollapsed_FSMLength/abundance.tsv"),
             paste0(root_dir,"K24_TargetedReadsCollapsed_FSMExp/abundance.tsv"),
             paste0(root_dir,"K24_AllReadsNotCollapsed/abundance.tsv"))

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.70, 0.7),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 10,family="ArialMT"),
                 axis.text.x= element_text(size=12,  family="ArialMT"),
                 axis.text.y= element_text(size=12,  family="ArialMT"))


input_files = lapply(files, function(x) read.table(x, header = T))
names(input_files) = c("AllReadsCollapsed_FSMLength","AllReadsCollapsed_FSMExp",
                       "TargetedReadsCollapsed_FSMLength","TargetedReadsCollapsed_FSMExp","AllReadsNotCollapsed")

class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"
class.files = read.table(class.names.files, sep = "\t", as.is = T, header = T)

input_files = lapply(input_files, function(x) merge(x, class.files[,c("isoform","associated_transcript","associated_gene")], by.x = "target_id", by.y = "isoform"))

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
# SQANTI filtered file 
sqanti.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt"
sqanti.class.files <- SQANTI_class_preparation(sqanti.names.files,"standard")

# SQANTI file of merged whole and targeted transcriptome 
merged.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/Whole_Targeted/merged_whole_targeted_classification.txt"
merged.class.files <- read.table(merged.names.files, header = T)

# Tama-merged corresponding isoforms
merged_targetedid <- read.table(paste0(diff_dir,"Whole_Targeted/merged_targetedPBID.txt"), header = T)
merged_wholeid <- read.table(paste0(diff_dir,"Whole_Targeted/merged_wholePBID.txt"), header = T)

reannotate_tamamerged_output <- function(tama_merged_transcript_file){
  # merge the Norm_transcounts from RNA-Seq Expression with tama merge annotation file, with the corresponding Pbid from targeted transcriptome 
  colnames(tama_merged_transcript_file)[[1]] = "tama_id"  # initially called "isoform" 
  dat = merge(tama_merged_transcript_file, merged_targetedid[,c("transcript_id","isoform")], by.x = "tama_id", by.y = "transcript_id", all.x = T)
  
  # note there are some transcripts that were not detected in targeted transcriptome but were in whole transcriptome and were aligned with RNA-Seq
  targeted_detected = dat[!is.na(dat$isoform),]  # these are the isoforms that were detected 
  missing_isoforms = unique(dat[is.na(dat$isoform),"tama_id"]) # these are the isoforms that were not detected
  
  # find the corresponding pbid from the whole transcriptome 
  wholeco_missing_isoforms = merged_wholeid[merged_wholeid$transcript_id %in% missing_isoforms,c("transcript_id","isoform")]
  
  # remerge the Norm_transcounts with the list of reannotated pbid from whole transcriptome
  targeted_not_detected = merge(dat %>% select(-c("isoform")),wholeco_missing_isoforms,by.x = "tama_id", by.y = "transcript_id", all.y = T) %>% 
    mutate(isoform = paste0(isoform,"_Whole"))
  
  combined = rbind(targeted_detected,targeted_not_detected)    
  return(combined)
}

tama_targetgene = merged.class.files[,c("isoform","associated_gene","associated_transcript")] %>% filter(toupper(associated_gene) %in% TargetGene)
merged_WholeTargeted = read.table(paste0(root_dir,"K24_Merged_WholeTargeted/abundance.tsv"),sep = "\t", as.is = T, header = T)
input_files$Merged_WholeTargeted = reannotate_tamamerged_output(merged_WholeTargeted %>% filter(target_id %in% tama_targetgene$isoform))
input_files$Merged_WholeTargeted = merge(input_files$Merged_WholeTargeted ,tama_targetgene[,c("isoform","associated_gene","associated_transcript")], by.x = "tama_id", by.y = "isoform")
colnames(input_files$Merged_WholeTargeted)[1] <- "target_id"

sqanti_plot_counts <- function(dat, input_title){
  p = ggplot(dat, aes(x = reorder(isoform, -FL), y = FL,  fill = associated_transcript)) + geom_bar(stat = "identity") + labs(x = "PB.ID", y = "Total Iso-Seq FL Counts", title = input_title) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + mytheme + scale_fill_discrete(name = "Associated Transcript")
  
  return(p)
}

plot_counts <- function(dat, input_title){
  p = ggplot(dat, aes(x = reorder(target_id, -est_counts), y = est_counts,  fill = associated_transcript)) + geom_bar(stat = "identity") + labs(x = "PB.ID", y = "RNA-Seq Estimated Counts", title = input_title) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + mytheme + scale_fill_discrete(name = "Associated Transcript") 
  
  
  return(p)
}

plot_tpm <- function(dat, input_title){
  p = ggplot(dat, aes(x = reorder(target_id, -tpm), y = tpm,  fill = associated_transcript)) + geom_bar(stat = "identity") + labs(x = "PB.ID", y = "TPM", title = input_title) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p)
}

plot_pergene <- function(input_gene){
  dat = lapply(input_files,function(x) x %>% filter(associated_gene == input_gene) %>% arrange(-est_counts) %>% .[1:10,]) 
  p <- plot_grid(plot_counts(dat$AllReadsNotCollapsed,paste0(input_gene, ": All Reads Not Collapsed")),
            plot_counts(dat$AllReadsCollapsed_FSMLength,"All Reads Collapsed by longest FSM"),
            plot_counts(dat$AllReadsCollapsed_FSMExp,"All Reads Collapsed by most abundant FSM"),
            plot_counts(dat$Merged_WholeTargeted,"Merged Whole & Targeted Reads \n Collapsed by most abundant FSM"),
            plot_counts(dat$TargetedReadsCollapsed_FSMLength,"Targeted Reads Collapsed by longest FSM"),
            plot_counts(dat$TargetedReadsCollapsed_FSMExp,"Targeted Reads Collapsed by most abundant FSM"))
  
  return(p)
}

plottpm_pergene <- function(input_gene){
  dat = lapply(input_files,function(x) x %>% filter(associated_gene == input_gene) %>% arrange(-est_counts) %>% .[1:10,]) 
  p <- plot_grid(plot_tpm(dat$AllReadsNotCollapsed,paste0(input_gene, ": All Reads Not Collapsed")),
                 plot_tpm(dat$AllReadsCollapsed_FSMLength,"All Reads Collapsed by longest FSM"),
                 plot_tpm(dat$AllReadsCollapsed_FSMExp,"All Reads Collapsed by most abundant FSM"),
                 plot_tpm(dat$Merged_WholeTargeted,"Merged Whole & Targeted Reads \n Collapsed by most abundant FSM"),
                 plot_tpm(dat$TargetedReadsCollapsed_FSMLength,"Targeted Reads Collapsed by longest FSM"),
                 plot_tpm(dat$TargetedReadsCollapsed_FSMExp,"Targeted Reads Collapsed by most abundant FSM"))
  
  return(p)
}

# novel vs known isoforms first isoform 
input_files$AllReadsCollapsed_FSMLength %>% filter(associated_gene %in% TargetGene) 

output = data.frame()
for(g in 1:length(TargetGene)){
  output[g,1] = TargetGene[g]
  output[g,2] = input_files$AllReadsCollapsed_FSMLength %>% filter(associated_gene == TargetGene[g]) %>% arrange(-est_counts) %>% .[1,c("associated_transcript")]
  output[g,3] = input_files$AllReadsCollapsed_FSMExp %>% filter(associated_gene == TargetGene[g]) %>% arrange(-est_counts) %>% .[1,c("associated_transcript")]
  output[g,4] = input_files$TargetedReadsCollapsed_FSMLength %>% filter(associated_gene == TargetGene[g]) %>% arrange(-est_counts) %>% .[1,c("associated_transcript")]
  output[g,5] = input_files$TargetedReadsCollapsed_FSMExp %>% filter(associated_gene == TargetGene[g]) %>% arrange(-est_counts) %>% .[1,c("associated_transcript")]
  output[g,6] = input_files$AllReadsNotCollapsed %>% filter(associated_gene == TargetGene[g]) %>% arrange(-est_counts) %>% .[1,c("associated_transcript")]
}
colnames(output) = c("Gene","AllReadsCollapsed_FSMLength","AllReadsCollapsed_FSMExp","TargetedReadsColapsed_FSMLength","TargetedReadsCollapsed_FSMExp","AllReadsNotCollapsed")
reshape2::melt(output, id = "Gene") %>% mutate(Transcript_type = ifelse(value == "novel","novel","known")) %>% group_by(variable,Transcript_type) %>% tally() %>% 
  ggplot(., aes(x = variable, y = n, fill = Transcript_type)) + geom_bar(stat = "identity") + mytheme + labs(y = "Number of Target Genes", x = "Method") + theme(legend.position = "top")

pdf(paste0(output_plot_dir,"/TargetedDifferentialAnalysis_testing_panelled.pdf"), width = 15, height = 10)
for(i in TargetGene){print(plot_pergene(i))}
dev.off()

pdf(paste0(output_plot_dir,"/TargetedDifferentialAnalysis_testing_panelledTPM.pdf"), width = 15, height = 10)
for(i in TargetGene){print(plottpm_pergene(i))}
dev.off()

dat = input_files$AllReadsCollapsed_FSMExpBTS %>% filter(associated_gene == "Snca") %>% arrange(-est_counts) %>% .[1:10,]
#plot_counts(dat,"Targeted Reads Collapsed by most abundant FSM")


pdf(paste0(output_plot_dir,"/TargetedDifferentialAnalysis_testing_finaloutput.pdf"), width = 20, height = 5)
genes = c("Bin1","Clu")
for(g in genes){
  print(plot_grid(sqanti_plot_counts(sqanti.class.files %>% filter(associated_gene == g) %>% arrange(-FL) %>% .[1:10,],paste(g,"Iso-Seq Expression")),
            plot_counts(input_files$AllReadsNotCollapsed %>% filter(associated_gene == g) %>% 
                          arrange(-est_counts) %>%.[1:10,],paste(g,"RNA-Seq Expression \nAll Reads, Not Collapsed")),
            plot_counts(input_files$AllReadsCollapsed_FSMLength %>% filter(associated_gene == g) %>% 
                          arrange(-est_counts) %>%  .[1:10,],paste(g,"RNA-Seq Expression \nAll Reads, longest FSM Collapse")), labels = "auto", ncol = 3, label_size = 30,scale = 0.9))
}
dev.off()


