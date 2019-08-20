# Szi Kay Leung
# 20/08/2019: 

#********************** Packages (install if not found)
list_of_packages <- c("DT","stringr")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

suppressMessages(library(DT))
suppressMessages(library(stringr))
#********************** Tabualate details for sequenced rTg4510 sample
input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Samples"
tg4510 <- read.csv(print(paste0(input_dir,"/","Tg4510_fullsample.csv")))[,-1]

sequenced <- c("Q21","O18", "L22","K18","O23","S23","S18","K17","M21","K23","Q20","K24")
sequenced_info <- subset(tg4510, tg4510$Sample.ID %in% sequenced)
#sequenced_info[order(-sequenced_info$ng.ul),]

Info <- datatable(sequenced_info) %>%
            formatStyle(
              'Genotype', 
              backgroundColor = styleEqual(
                unique(sequenced_info$Genotype), c('lightblue', 'lightgreen')
                )
              ) %>%
            formatStyle(
              'Age_in_months', 
              backgroundColor = styleEqual(
                unique(sequenced_info$Age_in_months), c('orange', 'yellow')
              )
            )

print(Info)

#********************** Read Sequel output from sequenced rTg4510 samples 
sequel_output <- read.csv(print(paste0(input_dir,"/","Sequenced_mouse_output.csv")), header = TRUE, stringsAsFactors=FALSE)[-1]

#********************** Merge Details and Sequel output
final <- merge(sequenced_info,sequel_output, by.x ="Sample.ID", by.y = "Sample")

headers <- c("Sample.ID","Genotype","Age_in_months", "RIN", "ng.ul","Mouse","Total.Bases..GB.","Unique.Molecular.Yield..GB.","P1_Counts","P1_Perc")
datatable(final[c(headers)]) %>%
  formatStyle(
    'Genotype', 
    backgroundColor = styleEqual(
      unique(final$Genotype), c('lightblue', 'lightgreen')
    )
  ) %>%
  formatStyle(
    'Age_in_months', 
    backgroundColor = styleEqual(
      unique(final$Age_in_months), c('orange', 'yellow')
    )
  )

#********************** Sample Selection 
all <- merge(tg4510, sequel_output, by.x ="Sample.ID", by.y = "Sample", all.x = TRUE)
datatable(all[c(headers)]) %>%
  formatStyle(
    'Genotype', 
    backgroundColor = styleEqual(
      unique(all$Genotype), c('lightblue', 'lightgreen')
    )
  ) %>%
  formatStyle(
    'Age_in_months', 
    backgroundColor = styleEqual(
      unique(all$Age_in_months), c('orange', 'yellow',"pink","red")
    )
  )
  


# # prepare P0, P1 and P2 to separate out percentages and counts
# Ps <- c("P0", "P1", "P2")
# for (i in Ps){
#   input <- sequel_output[,c(i)]
#   perc <- word(input, c(1), sep = fixed ('('))
#   perc <- word(perc, c(1), sep = fixed ('%') )
#   perc <- paste0(perc,"%")
#   
#   counts <- word(input, c(2), sep = fixed ('('))
#   counts <- word(counts, c(1), sep = fixed (')'))
#   counts <- as.numeric(gsub(",","",counts))
#   
#   sequel_output[[paste0(i,"_Perc")]] <- perc
#   sequel_output[[paste0(i,"_Counts")]] <- counts
# }
# # remove initial P0, P1, and P2 containing percentages and counts
# sequel_output <- sequel_output[,-c(12:14)]
# 
# Mouse <- sequel_output[c(1:8),c(-1)]
# write.csv(Mouse,print(paste0(input_dir,"/","Sequenced_mouse_output.csv")))
