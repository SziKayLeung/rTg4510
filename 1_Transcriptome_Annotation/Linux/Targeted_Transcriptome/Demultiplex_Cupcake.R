#!/usr/bin/env Rscript
<<<<<<< HEAD
# script.R <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
# script.R /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq/REFINE /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/TOFU/All_Targeted_Merged.collapsed.read_stat.txt

args = commandArgs(trailingOnly=TRUE)
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))

refine_dir <- args[1]       # refine directory containing flnc.report.csv
input_cluster_report <- args[2]   # all merged clustered report.csv
input_tofu_readstat <- args[3]
output_path_file <- args[4]


# input flnc report
flnc_report <- list.files(path = refine_dir, pattern = "*flnc.report.csv", full.names = T)
flnc_report_names <- list.files(path = refine_dir, pattern = "*flnc.report.csv")
cat("Processing:", flnc_report, "\n")
=======

library("stringr")
library("dplyr")
library("tidyr")

refine_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq/REFINE"
flnc_report <- list.files(path = refine_dir, pattern = "*flnc.report.csv", full.names = T)
flnc_report_names <- list.files(path = refine_dir, pattern = "*flnc.report.csv")
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5

report <- lapply(flnc_report, function(x) read.csv(x))
names(report) <- flnc_report_names

# check if there are any overlapping ids across the samples
all <- bind_rows(report, .id = "column_label")
all$sample <- word(all$column_label,c(1), sep = fixed("."))
n_occur <- data.frame(table(all$id))
n_occur[n_occur$Freq > 1,]

# check ids from flnc report is within the clustered.csv
<<<<<<< HEAD
cat("Using:", input_cluster_report, "\n")
cluster_report <- read.csv(input_cluster_report)
=======
cluster_report <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv")
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5

# clusterd report contains a subset if id from flnc report
# not included ids are lowly clustered
setdiff(all$id, cluster_report$read_id)
setdiff(cluster_report$read_id,all$id)

# check ids from flnc report same as the read_stats.txt from All merged data
<<<<<<< HEAD
cat("Using:", input_tofu_readstat, "\n")
read_stats <- read.table(input_tofu_readstat, header = TRUE)
=======
read_stats <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/TOFU/All_Targeted_Merged.collapsed.read_stat.txt",
                         header = TRUE)
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5

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

<<<<<<< HEAD
cat("Output:", output_path_file, "\n")
write.csv(final_output,output_path_file,row.names = FALSE, quote = FALSE)
=======
write.csv(final_output,"/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt",
          row.names = FALSE, quote = FALSE)
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5
