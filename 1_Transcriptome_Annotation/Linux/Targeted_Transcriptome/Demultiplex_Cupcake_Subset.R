#!/usr/bin/env Rscript

library("stringr")
library("dplyr")
library("tidyr")

refine_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq/REFINE"
flnc_report <- list.files(path = refine_dir, pattern = "*flnc.report.csv", full.names = T)
flnc_report_names <- list.files(path = refine_dir, pattern = "*flnc.report.csv")

# do not include batched report i.e. Targeted_Seq_1.flnc.report.csv
flnc_report <- flnc_report[!grepl("Targeted_Seq", flnc_report)]

report <- lapply(flnc_report, function(x) read.csv(x))
names(report) <- flnc_report_names[!grepl("Targeted_Seq", flnc_report_names)]

# check if there are any overlapping ids across the samples
all <- bind_rows(report, .id = "column_label")
all$sample <- word(all$column_label,c(1), sep = fixed("."))
n_occur <- data.frame(table(all$id))
n_occur[n_occur$Freq > 1,]

# check ids from flnc report is within the clustered.csv
cluster_report <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged_Subset/CLUSTER/All_Targeted_Merged_Subset.clustered.cluster_report.csv")

# clusterd report contains a subset if id from flnc report
# not included ids are lowly clustered
setdiff(all$id, cluster_report$read_id)
setdiff(cluster_report$read_id,all$id)

# check ids from flnc report same as the read_stats.txt from All merged data
read_stats <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged_Subset/TOFU/All_Targeted_Merged_Subset.collapsed.read_stat.txt", header = TRUE)

# id from clustered report same as the read_stats
setdiff(cluster_report$read_id, read_stats$id)
setdiff(read_stats$id,cluster_report$read_id)

# merged stats
merged_stats <- merge(read_stats, all[,c("id","sample")], by = "id")

# Unique Mapped Reads and FL
unique <- merged_stats[merged_stats$stat == "unique",]
final <- unique[,c("sample","pbid")]
final_counts <- final %>% group_by(sample,pbid) %>% tally()
final_output <- spread(final_counts, sample,n)
colnames(final_output)[1] <- "id"
final_output[is.na(final_output)] <- 0

write.csv(final_output,"/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged_Subset/TOFU/All_Targeted_Merged_Subset.Demultipled_Abundance.txt", row.names = FALSE, quote = FALSE)
