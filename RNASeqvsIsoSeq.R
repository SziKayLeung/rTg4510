library(stringr) # Extract strings
library(dplyr)   # dataframe splitting
library("ggpubr") # correlation plot


SQANTI_input <- function(sample){
  ## Define directory input
  sqanti_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/SQANTI2"
  
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

TOFU_input <- function(sample){
  ## Define directory input
  tofu_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/ToFU"
  
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
  return(head(abundance_dat))
}


Annotate2Abundance_IsoSeq <- function(sample,classification_dat,abundance_dat){
  # Merge SQANTI_Classssification file for genename with TOFU_Abundance file for FL counts
  Merge_IsoSeq <- merge(classification_dat,abundance_dat,by.x = "isoform", by.y = "PacBio_Id")
  #write.csv(Merge, file = print(paste0(sample,"Merged.csv")))
  
  # save Merge file to workspace: https://stackoverflow.com/questions/32563153/function-save-returned-data-frame-to-workspace
  print(paste("Merged file of SQANTI Classification and Abundance File of Sample", sample))
  Merge_IsoSeq <<- Merge_IsoSeq
  return(Merge_IsoSeq)
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
  return(validation)
}

Input_RNAseq <- function(sample){
  # Define directory input 
  featurecount_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/WT8_ISOSEQ/Post/FeatureCounts"
  
  print("Input FeatureCount for All Samples")
  # Import featurecount file with annotated mgisymbol names
  featurecount_file <- list.files(path = featurecount_dir, pattern = "Final_WT8_genenames.mx", full.names = TRUE) 
  featurecount_dat <- read.delim2(featurecount_file, row.names = 1, header = TRUE)
  # Featurecount_dat <<- featurecount_dat #save to workspace
  print(head(featurecount_dat))
  
  # Grep specific sample column with counts and ensembl_name (column 1)
  RNASeq <- data.frame(featurecount_dat[,sample],row.names = (rownames(featurecount_dat)))
  colnames(RNASeq)[1] <- paste("RNASeq",sample,"Raw Counts")
  print(paste("Input FeatureCount for Sample", sample))
  #print(head(RNASeq))
  
  # To extract the mgi_symbol separated by _ per row (carried out with a for loop)
  Genes <- data.frame(lapply(rownames(RNASeq), function(x){gene<-word(x,2, sep = fixed ('_'))}))
  Genes <- data.frame(t(Genes))
  colnames(Genes)[1] <- "Mgi_Symbol"
  RNASeq <- cbind(RNASeq,Genes)
  
  # Set a threshold for RNA-Seq FeatureCounts 
  RNASeq <- RNASeq[which(RNASeq$`RNASeq K17 Raw Counts` > 5),] 
  RNASeq$Mgi_Symbol <- as.character(RNASeq$Mgi_Symbol) #essential for correct merging downstream 
  
  RNASeq <<- RNASeq
  return(head(RNASeq))
}

Merge_RNASeq_Isoseq <- function(IsoSeq,RNASeq){
  Full_Merge <- merge(IsoSeq,RNASeq,by.x = "associated_gene",by.y = "Mgi_Symbol",all = TRUE)
  Full_Merge <<- Full_Merge
  return(head(Full_Merge))
}

AD_Counts <- function(AD_Genes){
  
  # Create new dataframe for appending output from loop
  AD_Genes_Merged <- data.frame() 
  
  for(i in AD_Genes){
    df <- Full_Merge[which(Full_Merge$associated_gene == i),]
    AD_Genes_Merged <- rbind(AD_Genes_Merged, df)
  }
  
  return(AD_Genes_Merged)
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
  
  # Substitute NA values for 0 to plot correlation
  Full_Merge$Isoseq_FL_Counts[is.na(Full_Merge$Isoseq_FL_Counts)] <- 0
  Full_Merge$RNASeq_Raw_Counts[is.na(Full_Merge$RNASeq_Raw_Counts)] <- 0
  
  Full_Merge <<- Full_Merge #important to update df 
  return(head(Full_Merge))
}

Run_Corplot <- function(dat){
  Corplot <- ggscatter(dat, y = "RNASeq_Raw_Counts", x = "Isoseq_FL_Counts", 
                       color = "Detection",
                       palette = c("red", "black"," blue"),
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson", 
                       ylab = "RNASeq_Raw_Counts", xlab = "PacBio_FL_Counts")
  
  return(Corplot)
}

Missing_Reads_Plot <- function(dat){
  # Only genes with either IsoSeq/RNASeq Reads
  Missing <- dat[which(dat$Detection == "noIsoSeq_RNASeq" | dat$Detection == "IsoSeq_noRNASeq" ),]
  
  Missing_Plot <- ggscatter(Missing, y = "RNASeq_Raw_Counts", x = "Isoseq_FL_Counts", 
                            color = "Detection",
                            palette = c("red", " blue"),
                            ylab = "RNASeq_Raw_Counts", xlab = "PacBio_FL_Counts")
  
  Missing <<- Missing # for downstream use
  return(Missing_Plot)
}

Missing_Genes <- function(value){
  print(paste("Genes with no IsoSeq Reads but RNASeq RawCounts >", value))
  print(Missing[which(Missing$RNASeq_Raw_Counts > value),1])
  print("Genes with only IsoSeq Reads, and no RNASeq Reads")
  print(Missing[which(Missing$Detection == "IsoSeq_noRNASeq"),1])
}