library(knitr)
library(readr)
library(dplyr)


samplelist <- c("B21", "C21", "K17", "K23", "M21", "O23", "Q21", "S23")

for (sample in samplelist) {
  
  rmarkdown::render(input = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq3_Tg4510/RNASeqvsIsoSeq_Sample.Rmd", 
                    output_file = paste0(sample,"_IsoseqRNASeq.html"),
                    output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213sl693/RNASeq/Correlations")
}


output_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/Correlations/"
all_Full_Merge <- list.files(path = output_dir , 
                                       pattern = "_Full_Merge.csv", full.names = TRUE)

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq3_Tg4510/RNASeqvsIsoSeq.R")

for (file in all_Full_Merge){
  print(file)
  
  df <- read.csv(file, header = TRUE)[,-1] #remove first column of index from import
  
  print(cor.test(df$Isoseq_FL_Counts, df$RNASeq_Raw_Counts, method = "pearson"))
  
  #pdf(paste0(output_dir,"AllSamples.pdf"))
  Run_Corplot(df)
  print(Corplot)
  #print(Missing_Reads_Plot(df))
  #dev.off()
}
  
  Missing_Genes(5000)
  
}

df <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/Correlations//B21_Full_Merge.csv",
              header = TRUE)





# Running individual samples through all the functions 
# for (sample in samplelist) {
#   sample
#   source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq3_Tg4510/RNASeqvsIsoSeq.R")
#   SQANTI_input(sample)
#   TOFU_input(sample)
#   Annotate2Abundance_IsoSeq(sample,classification_dat,abundance_dat)
#   Input_RNAseq(sample)
#   Merge_RNASeq_Isoseq(Merge_IsoSeq_SumFL,RNASeq) # output file name = Full_Merge
#   Missing_Reads_Review()
# }

