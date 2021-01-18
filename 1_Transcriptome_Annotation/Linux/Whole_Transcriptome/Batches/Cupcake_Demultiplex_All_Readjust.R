#!/usr/bin/env Rscript
# Szi Kay Leung
# 04/01/2021: Read in collapsed file from merged whole transcriptome mouse data to demultiplex using sample ID (ERCC)
# only 10 samples had ERCCs (first 2 samples did not, Q21, O18)

library("stringr")

# read in collapsed file from merged data
stats <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/TOFU_ADJUST/All_Merged.collapsed.read_stat.txt", header = TRUE)

# extract the run_id from the first part of the id
stats$run_id <- word(stats$id,c(1), sep = fixed("/"))
# unique(stats$run_id)

# description of sample to run_id
sample <- data.frame(sample <- c("Q21_WT","O18_TG","L22_TG","K18_TG","O23_WT","S23_WT","S18_TG","K17_WT","M21_WT","K23_WT","Q20_TG","K24_TG"),
                     run_id <- c("m54082_180607_173058","m54082_180605_141944",
                                 "m54082_190306_083150","m54082_190307_045507",
                                 "m54082_190401_165425","m54082_190403_135102",
                                 "m54082_190404_101400","m54082_190405_063832",
                                 "m54082_190430_163756","m54082_190524_145911",
                                 "m54082_190527_173356","m54082_190529_082942"))
colnames(sample) <- c("sample","run_id")

# merge the data to find the run_id to the sample_id
merge_final <- merge(stats, sample, by = "run_id")

# Filter only uniquely mapped
merge_final <- merge_final %>% filter(stat == "unique")

# tally the number of reads per sample and final output
final <- merge_final[,c("sample","pbid")]
final_counts <- final %>% group_by(sample,pbid) %>% tally()
final_output <- spread(final_counts, sample,n)
colnames(final_output)[1] <- "id"
final_output[is.na(final_output)] <- 0

write.csv(final_output,"/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/TOFU_ADJUST/Mouse.Demultiplexed_Abundance.txt",
          row.names = FALSE, quote = FALSE)
