# 04/11/2020: Using Aaron's script from Iso-Seq paper to find differential transcript usage between WT and TG

##### Analysis #####
# 1. Input SQANTI2 classification file from chained data FL count of all samples (n = 12)
# 2. Convert FL count per sample to TPM using FL_count/sum(total_counts) * million
# 3. Wilcoxin-rank sum test on TPM between WT and TG
# 4. Calculte mean expression of TPM per transcript between WT and TG
# 5. Convert gene name from ensembl id to gene symbol using reference mm10 gtf
# 6. Extract transcripts with significant wilcoxin p value (u < 0.05)
# 6. Convert gene name from ensembl id to gene symbol using reference mm10 gtf
# 7. Only keep genes if differential transcript usage is observed in both WT and TG
# 8. Subset only genes with significant differential expression and genes observed in both genotype
# 9. Output Plots
# 10. Output Tables
###################

# Packages
library("dplyr")
library("ggplot2")
library("magicfor") # save in list
library("tidyr") # wide to long

# functions for themes
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq3_Tg4510/Isoseq3/Isoseq3_QC/Rmarkdown_Input.R")
# output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Rmarkdown"
##### 1. Input SQANTI2 classification file
input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/SQANTI/"
class_file <- read.table(paste0(input_dir,"WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.txt"),header=T)

# Sample Phenotype classification
TG <- c("O18","K18","S18","L22","Q20","K24")
WT <- c("Q21","K17","M21","O23","S23","K23")

##### 2. Convert FL count per sample to TPM
colnames(class_file)
# columns 46 to 56 (last column) of df contain FL counts
DTU <- apply(class_file[,46:length(class_file)], 2, function(x) x/sum(x)*1000000) %>%
  as.data.frame() %>%
  # rename columns with phenotype + TPM
  dplyr::rename_all(paste0, "_TPM")


##### 3. Wilcoxin-rank sum test on TPM

# smallest dected significance assuming completely transcript expression = 0.002165
# df <- data.frame("WT" = c(1:6), "TG" = c(51:56))
# wilcox.test(df$WT,df$TG).pvalue


# Name columns as WT and TG for downstream analysis
colnames(DTU) <- lapply(colnames(DTU), function(x)
  if(grepl(paste(WT, collapse="|"), x)){
    paste0(x, "_WT_TPM")
  } else if (grepl(paste(TG, collapse="|"), x)){
    paste0(x, "_TG_TPM")
  } else {
    paste0("NA")})

# order dataframe from columns 1:6 as WT and 7:12 as TG
DTU <- cbind(DTU[,grepl("WT",colnames(DTU))],DTU[,grepl( "TG",colnames(DTU))])

# colnames(DTU) # check column names

# Wilcoxin Test: WT vs TG
DTU <- DTU %>% mutate(u = apply(DTU, 1, function(df) wilcox.test(df[1:6],df[7:12])$p.value))
class_relevant_cols <- c("isoform","associated_gene","associated_transcript","structural_category","subcategory")
DTU <- cbind(class_file[,class_relevant_cols],
             # reorder so final column "u" from mutate comes first before other sample counts
             DTU[,c(length(DTU),1:length(DTU)-1)])

##### 4. Calculate mean & median expression per transcript for WT vs TG
FL_count_colnames <- c("FL.K17_TPM_WT_TPM", "FL.K23_TPM_WT_TPM","FL.M21_TPM_WT_TPM",
                       "FL.O23_TPM_WT_TPM", "FL.Q21_TPM_WT_TPM", "FL.S23_TPM_WT_TPM",
                       "FL.K18_TPM_TG_TPM", "FL.K24_TPM_TG_TPM", "FL.L22_TPM_TG_TPM",
                       "FL.O18_TPM_TG_TPM", "FL.Q20_TPM_TG_TPM", "FL.S18_TPM_TG_TPM"  )
DTU <- DTU %>%
  mutate(WT_mean_exp = apply(.[,grepl( "WT",colnames(.))], 1, function(x) mean(x))) %>%
  mutate(TG_mean_exp = apply(.[,grepl( "TG",colnames(.))], 1, function(x) mean(x))) %>%
  mutate(WT_median_exp = apply(.[,grepl( "WT",colnames(.))], 1, function(x) median(x))) %>%
  mutate(TG_median_exp = apply(.[,grepl( "TG",colnames(.))], 1, function(x) median(x))) %>%
  mutate(mean_exp_diff = abs(WT_mean_exp - TG_mean_exp)) %>%
  mutate(Direction = ifelse(WT_mean_exp > TG_mean_exp, "WT", "TG")) %>%
  .[,c("isoform","associated_gene","associated_transcript","structural_category","subcategory",
       "u","WT_median_exp","TG_median_exp","WT_mean_exp","TG_mean_exp","mean_exp_diff","Direction",FL_count_colnames)]


##### 5. Convert gene name from ensembl id to gene symbol using reference mm10 gtf
# classification file associated_gene column in Ensembl id rather than gene name
# convert ensembl id using reference
mm10_reference_file <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.vM22_gene_annotation_table.txt"
mm10_reference <- read.table(mm10_reference_file,as.is = T, header=T, sep = "\t")
DTU <- merge(DTU, mm10_reference[c("gene_id","GeneSymbol")], by.x = ("associated_gene"), by.y= "gene_id", all.x = T) %>%
  .[,c("isoform","GeneSymbol","associated_gene","associated_transcript","structural_category","subcategory",
       "u","WT_mean_exp","TG_mean_exp","WT_median_exp","TG_median_exp","mean_exp_diff","Direction",FL_count_colnames)]

##### 6. Extract transcripts with significant wilcoxin p value (u < 0.05)
DTU_sg <- subset(DTU, DTU$u < 0.05)

##### 7. Tabulate the number of samples with transcripts observed
# columns 12:17 - WT, columns 18:23 TG
DTU_sg <- DTU_sg %>%
  mutate("WT_nonzero" = apply(.[,12:17],1,function(x)sum(x != 0))) %>%
  mutate("TG_nonzero" = apply(.[,18:23],1,function(x)sum(x != 0)))

DTU_sg_filter <- DTU_sg[,c(1:12,26,27)]

# plot of the number of samples of WT and TG for differentially expressed transcripts
p1 <- DTU_sg_filter %>% group_by(WT_nonzero, TG_nonzero) %>% count() %>%
  mutate(combo = paste0("WT",WT_nonzero,"_TG",TG_nonzero)) %>%
  mutate(WT_nonzero = factor(WT_nonzero, levels = c(0:6))) %>%
  mutate(TG_nonzero = factor(TG_nonzero, levels = c(0:6)))  %>%
  mutate(perc = n/sum(.$n) * 100) %>%
  ggplot(data = ., aes(x = WT_nonzero, y = TG_nonzero)) +
  geom_tile(aes(fill = perc)) +
  geom_text(aes(label = signif(perc,3)), color = "black", fontface = "bold", size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "NA", name = "Percentage of DE transcripts") +
  labs(x = "Number of WT Samples", y = "Number of TG samples")


# remove differentially expressed transcripts from 1 sample of WT or TG
DTU_sg_filter <- DTU_sg[DTU_sg$WT_nonzero != "1" & DTU_sg$TG_nonzero != "1",]

##### 7. Only keep genes if differential transcript usage is observed in both WT and TG
genes <- unique(DTU_sg$GeneSymbol)
#genes <- unique(DTU_sg_filter$GeneSymbol)

magic_for(print, silent = TRUE)
for(i in genes){
  genetoexamine <- DTU %>% filter(GeneSymbol == i)
  genetoexamine_Direction <- unique(genetoexamine$Direction)
  if(length(genetoexamine_Direction) == 2){print(i)}
}

genelist <- magic_result_as_vector() %>% .[!is.na(.)]

##### 8. Subset only genes with significant differential expression and genes observed in both genotype
dtu_list <- subset(DTU,DTU$GeneSymbol %in% genelist) %>%
  # log10 expression differece for cut-off
  mutate(log10_mean_exp_diff = log10(mean_exp_diff))
dtu_list[["transcript_name_id"]] <- paste0(dtu_list$associated_transcript,"_", dtu_list$isoform)

# Apply mean difference threshold and rank by biggest difference
# histogram of log10_exp_diff (absolute mean difference between TG and WT) --> cut off at 0 (i.e TPM diff = 1)
# hist(dtu_list$log10_exp_diff)
dtu_list_expcutoff <- subset(dtu_list, dtu_list$log10_mean_exp_diff > 0) %>% arrange(-mean_exp_diff)


##### 9. Output plots
# function to plot differential expression changes
plot_transcript_expression <- function(input_df, gene){

  # subset df of gene of interest and prepare for plot
  df <- input_df %>%
    # subset df to gene of interest
    filter(GeneSymbol == gene) %>%
    # only include relevant columns for plotting
    select(GeneSymbol,isoform, associated_transcript, structural_category, u, WT_mean_exp, TG_mean_exp) %>%
    # classify transcripts that are significant in plot
    mutate(sig = ifelse(u < 0.05, "yes","no")) %>%
    mutate(sig = factor(sig, levels = c("yes","no"))) %>%
    # wide to long for plotting
    gather(., Genotype, Expression, WT_mean_exp:TG_mean_exp, factor_key=TRUE)

  # plot
  p <- ggplot(df, aes(x = Genotype, y = Expression, group = isoform, colour = structural_category, linetype = sig)) +
    geom_line() + geom_point() +
    scale_linetype_manual(values=c("solid", "dotted"), name = "Wilcoxon Test", labels = c("P < 0.05", "P > 0.05")) +
    # only includ name of transcript from the WT end if differentially transcribed between WT and TG
    geom_text(data=filter(df,Genotype =="WT_mean_exp"),
              aes(label = ifelse(sig == "yes", associated_transcript, "")),
              nudge_x = 0.25) +
    labs(x = "", y = "Mean Transcript Expression (TPM)", title = paste(gene, "\n\n")) +
    mytheme + theme(legend.position = c(0.9,0.8), plot.title = element_text(face="italic")) +
    scale_x_discrete(labels= c("WT", "TG")) +
    scale_colour_discrete(name = "Structural Category", labels = c("FSM", "NIC", "NNC"))

  return(p)
}

## plot expression for all genes that have differential transcript expression between WT and TG
#pdf(paste0(output_plot_dir,"/All_Differential_Transcript_Expression.pdf"), width = 11, height = 8.5)
#for(gene in unique(dtu_list_expcutoff$GeneSymbol)){
  #print(plot_transcript_expression(dtu_list_expcutoff, gene))
#}
#dev.off()

## plot expression fot top 50 genes with differentially expressed transcripts
#pdf(paste0(output_plot_dir,"/Top50_Differential_Transcript_Expression.pdf"), width = 11, height = 8.5)
#for(gene in unique(dtu_list_expcutoff$GeneSymbol)[1:50]){
#  print(plot_transcript_expression(dtu_list_expcutoff, gene))
#}
#dev.off()


##### 10. Output Tables
write.csv(DTU, file="Iso_Seq_Aaron_DTU_All.csv")  # All transcripts with wilcoxin rank sum test results
write.csv(dtu_list, file="Iso_Seq_Aaron_DTU_GeneExpressionBothGenotype.csv") # Only genes with significant differential expression and genes observed in both genotype
