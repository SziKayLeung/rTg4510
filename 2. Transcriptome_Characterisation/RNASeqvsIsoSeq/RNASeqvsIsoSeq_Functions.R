# Szi Kay Leung 
# Aim: Functions Script to be called for RNAseqvsIsoSeq.Rmd 
# 19/03/2020: Cleaned up script

library(stringr) # Extract strings
library(dplyr)   # dataframe splitting
library("ggpubr") # correlation plot
library(plotly) # interactive graphs
library(ggplot2)
library(grid) # ggplot add corr value
library(gridExtra) # ggplot add corr value
library(ggthemes) # ggplot minimum theme

## Define all output directory
output_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/Correlations/"
  
SQANTI_input <- function(sample, sqanti_input_dir){
  ## Define directory input
  sqanti_dir <- sqanti_input_dir
  
  #print("All Classification files from SQANTI")
  all_classification_files <- list.files(path = sqanti_dir , pattern = ".collapsed.filtered.rep_classification.filtered_lite_classification.txt", full.names = TRUE)
  #print(all_classification_files)
  #return(classification.files)
  
  ## Grep specific file based on sample name
  print(paste("Input SQANTI Filter output file for Sample", sample))
  classification_file <- grep(sample, all_classification_files, fixed = TRUE, value = TRUE)
  #print(classification_file)
  
  ## Read specific sample file
  classification_dat <- read.delim2(file = classification_file, header = TRUE)
  classification_dat <<- data.frame(classification_dat)
  
  print(paste("SQANTI Classification file of Sample", sample))
  return(head(classification_dat))
}

TOFU_input <- function(sample, tofu_input_dir){
  ## Define directory input
  tofu_dir <- tofu_input_dir
  
  ## List all files in TOFU directory
  #print("All Abundance files from TOFU")
  all_abundance_files <- list.files(path = tofu_dir , pattern = ".collapsed.filtered.abundance.txt", full.names = TRUE)
  #print(all_abundance_files)
  
  ## Grep specific file based on sample name
  print(paste("Input SQANTI Filter output file for Sample", sample))
  abundance_file <- grep(sample, all_abundance_files, fixed = TRUE, value = TRUE)
  print(abundance_file)
  
  ## Read specific sample file
  abundance_dat <- read.table(file = abundance_file, header = TRUE)[-c(1:14),] #remove rows 1:14 of descriptor
  abundance_dat <<- data.frame(abundance_dat)
  colnames(abundance_dat) <- c("PacBio_Id","count_fl","count_nfl","count_nfl_amb","norm_fl","norm_nfl","norm_nfl_amb")
  
  # View specific sample file
  print(paste("Abundance file of Sample", sample))
  abundance_dat <<- abundance_dat 
  head(abundance_dat)
}


Annotate2Abundance_IsoSeq <- function(sample,classification_dat,abundance_dat){
  # Merge SQANTI_Classssification file for genename with TOFU_Abundance file for FL counts
  Merge_IsoSeq <- merge(classification_dat,abundance_dat,by.x = "isoform", by.y = "PacBio_Id", all.x = TRUE)
  #write.csv(Merge, file = print(paste0(sample,"Merged.csv")))
  
  # save Merge file to workspace: https://stackoverflow.com/questions/32563153/function-save-returned-data-frame-to-workspace
  print(paste("Merged file of SQANTI Classification and Abundance File of Sample", sample))
  Merge_IsoSeq <<- Merge_IsoSeq
  head(Merge_IsoSeq)
}

SumFLCounts <- function(dat){
  # input dat dataframe, group by "associated_gene" column (i.e. group by genes), and include 2 colums:
  # column 1 = labelled "PacBio Isoform" == isoform column; note [1] for referring to entire column
  # column 2 = labelled "PacBio FL Counts" == sum of "count_Fl" column
  
  Merge_IsoSeq_SumFL <- dat %>% group_by(associated_gene) %>%
    summarize(PacBio_Isoform = isoform[1], PacBio_FL_Counts = sum(count_fl))
  Merge_IsoSeq_SumFL$associated_gene <- as.character(Merge_IsoSeq_SumFL$associated_gene)
  Merge_IsoSeq_SumFL <<- Merge_IsoSeq_SumFL
  return(head(Merge_IsoSeq_SumFL))
}

Validation_SumFLCounts <- function(gene,dat){
  # Validation that the function works 
  print(paste("Validation of summing PacBio FL"))
  print(paste("Original input data from ToFU Abundance files for the Gene", gene))
  df <- dat[which(dat$associated_gene == gene),c(7,40)] #only print column 7 and 40 = "associated_gene" and "fl_count"
  
  print(paste("Summed PacBio FL count for the Gene", gene , "saved as new dataframe for downstream analysis"))
  # Merge_IsoSeq_SumFL is the output dataframe from SumFLCounts Function
  df2<- Merge_IsoSeq_SumFL[which(Merge_IsoSeq_SumFL$associated_gene == gene),] 
  
  validation <- list(df,df2)
  validation <<- validation
  validation
}

Input_RNAseq <- function(sample){
  # Define directory input 
  featurecount_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/WT8_ISOSEQ/Post/FeatureCounts"
  
  print("Input FeatureCount for All Samples")
  # Import featurecount file with annotated mgisymbol names
  featurecount_file <- list.files(path = featurecount_dir, pattern = "Final_WT8_genenames.mx", full.names = TRUE) 
  featurecount_dat <- read.delim2(featurecount_file, row.names = 1, header = TRUE)
  featurecount_dat <<- featurecount_dat #save to workspace
  print(head(featurecount_dat))
  
  # Grep specific sample column with counts and ensembl_name (column 1)
  #sample_with_quotes <- c(print(sample, sep=""))
  RNASeq <- data.frame(featurecount_dat[,sample],row.names = (rownames(featurecount_dat)))
  
  colnames(RNASeq)[1] <- paste("RNASeq",sample,"Raw Counts")
  print(paste("Input FeatureCount for Sample", sample))
  #print(head(RNASeq))

  # To extract the mgi_symbol separated by _ per row (carried out with a for loop)
  Genes <- data.frame(lapply(rownames(RNASeq), function(x){gene<-word(x,2, sep = fixed ('_'))}))
  Genes <- data.frame(t(Genes))
  colnames(Genes)[1] <- "Mgi_Symbol"
  RNASeq <- cbind(RNASeq,Genes)
  
  RNASeq <<- RNASeq

  # Set a threshold for RNA-Seq FeatureCounts 
  RNASeq <- RNASeq[which(RNASeq[,1] > 5),] #refers to column 1, do not specific name of column as different for each sample
  RNASeq$Mgi_Symbol <- as.character(RNASeq$Mgi_Symbol) #essential for correct merging downstream 
  
  RNASeq <<- RNASeq
  head(RNASeq)
}

AD_Counts <- function(AD_Genes){
  
  # Create new dataframe for appending output from loop
  AD_Genes_Merged <- data.frame() 
  
  for(i in AD_Genes){
    df <- Full_Merge[which(Full_Merge$associated_gene == i),]
    AD_Genes_Merged <- rbind(AD_Genes_Merged, df)
  }
  
  AD_Genes_Merged
}

Missing_Reads_Review <- function(){
  print(paste("Total Number of Genes in Full_Merge of IsoSeq and RNASeq:", dim(Full_Merge)[1]))
  
  # Create "Detection column" in Full_Merge depending on whether reads present/absent in technology
  colnames(Full_Merge)
  colnames(Full_Merge) <- c("associated_gene","PacBio_Isoform_ID","Isoseq_FL_Counts","RNASeq_Raw_Counts") #double check col names
  
  Full_Merge$Detection_Isoseq = ifelse(is.na(Full_Merge$Isoseq_FL_Counts),"noIsoSeq","IsoSeq")
  Full_Merge$Detection_RNAseq = ifelse(is.na(Full_Merge$RNASeq_Raw_Counts),"noRNASeq","RNASeq") 
  Full_Merge$Detection = paste(Full_Merge$Detection_Isoseq,Full_Merge$Detection_RNAseq,sep="_") #concatenate columns
  
  # Check factor level for detection
  Full_Merge$Detection <- as.factor(Full_Merge$Detection)
  # print(levels(Full_Merge$Detection))
  
  print(paste("Total Number of Genes Detected in IsoSeq AND RNASeq:",length(which(Full_Merge$Detection == "IsoSeq_RNASeq"))))
  print(paste("Total Number of Genes Detected in IsoSeq but not RNASeq:",length(which(Full_Merge$Detection == "IsoSeq_noRNASeq"))))
  print(paste("Total Number of Genes Detected in RNASeq but not IsoSeq:",length(which(Full_Merge$Detection == "noIsoSeq_RNASeq"))))
  
  Full_Merge <<- Full_Merge #important to upFull_Mergee df 
  head(Full_Merge)
  
  write.csv(Full_Merge, file = print(paste0(output_dir, sample,"_Full_Merge.csv"))) # for downstream
}

Log_Counts <- function(){
  Log_Full_Merge <- Full_Merge
  
  # Correlation of Raw Data to plot 0 values
  # Substitute NA values for 0 to plot correlation
  Full_Merge$Isoseq_FL_Counts[is.na(Full_Merge$Isoseq_FL_Counts)] <- 0
  Full_Merge$RNASeq_Raw_Counts[is.na(Full_Merge$RNASeq_Raw_Counts)] <- 0
  
  # Correlation of Log(Raw_Data)
  # Unable to log 0 counts as log0 = infinity, therfore keep values as "NA"
  Log_Full_Merge <- na.omit(Log_Full_Merge)
  
  # Log Counts
  Log_Full_Merge$Log_Isoseq_FL_Counts <- log10(Log_Full_Merge$Isoseq_FL_Counts)
  Log_Full_Merge$Log_RNASeq_Raw_Counts <- log10(Log_Full_Merge$RNASeq_Raw_Counts)
  
  Log_Full_Merge <<- Log_Full_Merge
  Full_Merge <<- Full_Merge
  print(head(Log_Full_Merge))
}


# https://stackoverflow.com/questions/50960339/create-ggplot2-function-and-specify-arguments-as-variables-in-data-as-per-ggplot
Run_Corplot <- function(dat,x.var,y.var){
  corr <- grobTree(textGrob(paste("R : ", round(cor(dat[,x.var], dat[,y.var]), 4) ), x = 0.05, y = 0.97, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "italic")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  Corplot <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) + 
    geom_point(size = 0.4) + geom_smooth(method=lm) + theme_tufte() + scale_color_manual(values=c("royalblue2", "grey40", "firebrick2")) + annotation_custom(corr) 

  print(Corplot)
}

# https://stackoverflow.com/questions/50960339/create-ggplot2-function-and-specify-arguments-as-variables-in-data-as-per-ggplot
Interactive_Log <- function(dat){
  Corplot <- ggplot(dat, aes_string(Log_IsoSeq_FL_Counts,Log_RNASeq_Raw_Counts)) + 
    geom_point() + geom_smooth(method=lm) + theme_tufte() 
  
  # plotly:https://plot.ly/ggplot2/user-guide/
  Corplot <- ggplotly(Corplot)

  Corplot <- style(Corplot, text= dat$associated_gene, hoverinfo="text+x+y")
  print(Corplot)
}



Missing_Reads_Plot <- function(dat){
  # Only genes with either IsoSeq/RNASeq Reads
  Missing <- dat[which(dat$Detection == "noIsoSeq_RNASeq" | dat$Detection == "IsoSeq_noRNASeq" ),]
  
  
  Missing_Plot <- ggplot(Missing, aes(x = Isoseq_FL_Counts, y = RNASeq_Raw_Counts, color = Detection)) + geom_point() + geom_smooth(method=lm) + 
    scale_color_manual(values=c("royalblue2", "firebrick2")) 
  
  Missing_Plot <- Missing_Plot + theme_tufte()
  
  # plotly:https://plot.ly/ggplot2/user-guide/
  #Missing_Plot <- ggplotly(Missing_Plot)
  
  #Missing_Plot <- style(Missing_Plot, text= dat$associated_gene, hoverinfo="text+x+y")
  
  
  Missing <<- Missing # for downstream use
  Missing_Plot
}


Missing_Genes <- function(dat,value){
  print(paste("Genes with no IsoSeq Reads but RNASeq RawCounts >", value))
  print(dat[which(dat$RNASeq_Raw_Counts > value),1])
  print("Genes with only IsoSeq Reads, and no RNASeq Reads")
  print(dat[which(dat$Detection == "IsoSeq_noRNASeq"),1])
}


density_plot <- function(dat,x.var,y.var){
  
  corr <- grobTree(textGrob(paste("R : ", round(cor(dat[,x.var], dat[,y.var]), 4) ), x = 0.05, y = 0.97, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "italic")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1) +
    theme_tufte() +
    labs(x = "Log (Isoseq Full Length Counts)", y = "Log (RNASeq Raw Counts)") + 
    geom_smooth(method=lm, colour = "black") 
  
  print(p)
}

#hist(Log_Full_Merge$Log_Isoseq_FL_Counts[Log_Full_Merge$Log_Isoseq_FL_Counts<10])

Run_FSM <- function(sample,dat){
  
  # Merge SQANTI_Classssification file for genename with TOFU_Abundance file for FL counts
  FSM_Merge_IsoSeq <- merge(dat,abundance_dat,by.x = "isoform", by.y = "PacBio_Id", all.x = TRUE)
  #write.csv(Merge, file = print(paste0(sample,"Merged.csv")))
  
  # save Merge file to workspace: https://stackoverflow.com/questions/32563153/function-save-returned-data-frame-to-workspace
  print(paste("Merged file of SQANTI Classification and Abundance File of Sample", sample))
  FSM_Merge_IsoSeq <<- FSM_Merge_IsoSeq
  
  # sum FL Counts
  FSM_Merge_IsoSeq_SumFL <- FSM_Merge_IsoSeq %>% group_by(associated_gene) %>%
    summarize(PacBio_Isoform = isoform[1], PacBio_FL_Counts = sum(count_fl))
  FSM_Merge_IsoSeq_SumFL$associated_gene <- as.character(FSM_Merge_IsoSeq_SumFL$associated_gene)
  FSM_Merge_IsoSeq_SumFL <<- FSM_Merge_IsoSeq_SumFL
  
  # MERGE ISOSEQ with RNASEQ
  FSM_Full_Merge <- merge(FSM_Merge_IsoSeq_SumFL,RNASeq,by.x = "associated_gene",by.y = "Mgi_Symbol",all = TRUE)
  
  print(paste("Total Number of Genes in FSM_Full_Merge of IsoSeq and RNASeq:", dim(FSM_Full_Merge)[1]))
  
  # Create "Detection column" in FSM_Full_Merge depending on whether reads present/absent in technology
  colnames(FSM_Full_Merge)
  colnames(FSM_Full_Merge) <- c("associated_gene","PacBio_Isoform_ID","Isoseq_FL_Counts","RNASeq_Raw_Counts") #double check col names
  
  FSM_Full_Merge$Detection_Isoseq = ifelse(is.na(FSM_Full_Merge$Isoseq_FL_Counts),"noIsoSeq","IsoSeq")
  FSM_Full_Merge$Detection_RNAseq = ifelse(is.na(FSM_Full_Merge$RNASeq_Raw_Counts),"noRNASeq","RNASeq") 
  FSM_Full_Merge$Detection = paste(FSM_Full_Merge$Detection_Isoseq,FSM_Full_Merge$Detection_RNAseq,sep="_") #concatenate columns
  
  # Check factor level for detection
  FSM_Full_Merge$Detection <- as.factor(FSM_Full_Merge$Detection)
  # print(levels(FSM_Full_Merge$Detection))
  
  print(paste("Total Number of Genes Detected in IsoSeq AND RNASeq:",length(which(FSM_Full_Merge$Detection == "IsoSeq_RNASeq"))))
  print(paste("Total Number of Genes Detected in IsoSeq but not RNASeq:",length(which(FSM_Full_Merge$Detection == "IsoSeq_noRNASeq"))))
  print(paste("Total Number of Genes Detected in RNASeq but not IsoSeq:",length(which(FSM_Full_Merge$Detection == "noIsoSeq_RNASeq"))))
  
  FSM_Full_Merge <<- FSM_Full_Merge #important to update FSM_Full_Merge df 
  head(FSM_Full_Merge)
  
  #write.csv(FSM_Full_Merge, file = print(paste0(output_dir, sample,"_FSM_Full_Merge.csv"))) # for downstream
  
  # Log Counts
  Log_FSM_Full_Merge <- FSM_Full_Merge
  
  # Correlation of Raw Data to plot 0 values
  # Substitute NA values for 0 to plot correlation
  FSM_Full_Merge$Isoseq_FL_Counts[is.na(FSM_Full_Merge$Isoseq_FL_Counts)] <- 0
  FSM_Full_Merge$RNASeq_Raw_Counts[is.na(FSM_Full_Merge$RNASeq_Raw_Counts)] <- 0
  
  # Correlation of Log(Raw_Data)
  # Unable to log 0 counts as log0 = infinity, therfore keep values as "NA"
  Log_FSM_Full_Merge <- na.omit(Log_FSM_Full_Merge)
  
  # Log Counts
  Log_FSM_Full_Merge$Log_Isoseq_FL_Counts <- log10(Log_FSM_Full_Merge$Isoseq_FL_Counts)
  Log_FSM_Full_Merge$Log_RNASeq_Raw_Counts <- log10(Log_FSM_Full_Merge$RNASeq_Raw_Counts)
  
  Log_FSM_Full_Merge <<- Log_FSM_Full_Merge
  FSM_Full_Merge <<- FSM_Full_Merge
  print(head(Log_FSM_Full_Merge))
  
  density_plot(Log_FSM_Full_Merge,"Log_Isoseq_FL_Counts","Log_RNASeq_Raw_Counts")
  return(density_plot)
}

Prepare_Gencode <- function(gencode){
  # replace the first row as column name
  colnames(gencode) <- as.character(unlist(gencode[1,]))
  gencode <- gencode[-1, ]
  # removed white space in the column for grep
  gencode$GeneSymbol <- gsub('\\s+', '', gencode$GeneSymbol)
  
  print("Gencode vM20 Gene Annotation")
  datatable(gencode)
  
  # Sum up the number of transcripts in gencode
  # tabulate the freqency of unique input in GeneSymbol column  
  gencode$GeneSymbol = as.character(gencode$GeneSymbol)
  gencode_count <- data.frame(table(gencode$GeneSymbol))
  colnames(gencode_count) <- c("GeneName", "Gencode_Number")
  
  # convert to character for merging
  gencode_count$GeneName <- as.character(gencode_count$GeneName)
  #removed white space in the column for later merging
  gencode_count$GeneName <- gsub('\\s+', '', gencode_count$GeneName)
  
  gencode_count <<- gencode_count #save for later downstream function
}

Sum_Merge4TranscriptAbundance <- function(dat){
  # tabulate the freqency of unique input in associated_gene column 
  dat$associated_gene = as.character(dat$associated_gene)
  sqanti_count = data.frame(table(dat$associated_gene))
  
  colnames(sqanti_count) <- c("GeneName","Isoseq_Number")
  sqanti_count$GeneName <- as.character(sqanti_count$GeneName)
  datatable(sqanti_count)
  
  Transcript_Abundance <- merge(sqanti_count, gencode_count, by="GeneName",all.x = TRUE)
  Transcript_Abundance <- Transcript_Abundance[order(-Transcript_Abundance$Isoseq_Number),]
  Transcript_Abundance <<- Transcript_Abundance
  datatable(Transcript_Abundance)
}

disease_plot <- function(dat_disease){
  
  # gather the data for plotting; only column 2 and 3, to be sorted into Class and Number of Isoform 
  x_dat <- gather(dat_disease,"Class","Numbers_Of_Isoform",2:3)
  
  p <- ggplot(x_dat, aes(x=GeneName, y=Numbers_Of_Isoform, fill=Class)) + geom_bar(stat="identity", position = "dodge")
  
  p + theme(axis.text.x = element_text(size = 12, angle = 45,hjust = 1),) + scale_fill_discrete(name = "Class", labels = c("Gencode", "Isoseq")) + scale_fill_brewer(palette="Spectral") + labs(x = "Genes", y = "Counts" )
}

