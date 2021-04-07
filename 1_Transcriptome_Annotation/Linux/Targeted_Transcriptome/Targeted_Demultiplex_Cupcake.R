#!/usr/bin/env Rscript
# Szi Kay Leung
# 21/01/2020: script to demultiplex targeted data

### Rationale ###
# note script is different from cupcake_demultiplex.R as not able to use run_id for demultiplexing given that samples were barcoded in the same run (i.e samples from the same run would have the same id)
# instead use the id from flnc report generated from REFINE of the individual samples
# Rscript Targeted_Cupcake_Demultiplex.R <path/refine_dir> <cluster_report> <path/collapsed_read_stat file> <output_file>
###

### Arguments info ###
# <path/refine_dir> = refine directory containing all the individual sample repots from targeted data (flnc.report.csv)
# <cluster_report> = cluster_report.csv from merging all individual targeted samples (merged data, n = 24)
# <path/collapsed_read_stat file> = path of cupcake collapsed read stat file generated from merged targeted data
# <output_file> = output of abundance file
###

### Example Run in Linux ###
# REFINE_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq/REFINE
# Cluster_report=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv
# read_stat=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/TOFU/All_Targeted_Merged.collapsed.read_stat.txt
# output_file=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
# Rscript $REFINE $Cluster_report $read_stat $output_file

args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

refine_dir <- args[1]
input_cluster_report <- args[2]
input_cupcake_readstat <- args[3]
output_abuncance_file <- args[4]

cat("processing files from refine directory:", refine_dir, "\n")
flnc_report <- list.files(path = refine_dir, pattern = "*flnc.report.csv", full.names = T)
flnc_report_names <- list.files(path = refine_dir, pattern = "*flnc.report.csv")

# do not include batched report i.e. Targeted_Seq_1.flnc.report.csv
flnc_report <- flnc_report[!grepl("Targeted_Seq", flnc_report)]
flnc_report

report <- lapply(flnc_report, function(x) read.csv(x))
names(report) <- flnc_report_names[!grepl("Targeted_Seq", flnc_report_names)]

# check if there are any overlapping ids across the samples
cat("Overlapping ids across targeted samples \n")
all <- bind_rows(report, .id = "column_label")
all$sample <- word(all$column_label,c(1), sep = fixed("."))
n_occur <- data.frame(table(all$id))
n_occur[n_occur$Freq > 1,]

# check ids from flnc report is within the clustered.csv
cat("Processing cluster report:",input_cluster_report, "\n")
cluster_report <- read.csv(input_cluster_report)

# clusterd report contains a subset of id from flnc report
# not included ids are lowly clustered
cat("Same id in both cluster report and from flnc report of all the samples \n")
setdiff(all$id, cluster_report$read_id)
setdiff(cluster_report$read_id,all$id)

# check ids from flnc report same as the read_stats.txt from All merged data
cat("Processing cupcake collapsed read_stat file:",input_cupcake_readstat, "\n")
read_stats <- read.table(input_cupcake_readstat, header = TRUE)

# id from clustered report same as the read_stats
cat("Same id in both cluster report and from cupcake collapsed read_stat file generated in the merged data \n")
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

cat("Output file:",output_abuncance_file, "\n")
write.csv(final_output, output_abuncance_file, row.names = FALSE, quote = FALSE)
