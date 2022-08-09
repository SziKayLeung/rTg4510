# Szi Kay Leung 
# 11/01/2020: Check FL read counts across samples from demultiplexing one merged batch set 

##### Analysis #####
# 1. Number of FL reads in read_stat from cupcake collapse file per sample
# 2. Number of FL reads in abundance file created from the demultiplexing 
# 3. Number of FL reads from sqanti classification.txt (pre-filtering)
# 4. Number of FL reads from sqanti classification.txt (post-filtering)
# Using K17 sample: m54082_190405_063832, for reference

library("stringr")
'%ni%' <- Negate('%in%')

########## 1.read_stat from cupcake collapse file
# read in collapsed file from merged data
stats <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.read_stat.txt", header = TRUE)

# extract the run_id from the first part of the id
stats$run_id <- word(stats$id,c(1), sep = fixed("/"))

# number of FL reads 
K17_stats <- stats %>% filter(run_id == "m54082_190405_063832") %>% filter(stat == "unique") 

########## 2.abundance file
abundance <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt", header = TRUE, sep = ",") 

K17_abundance <- abundance[,c("K17_WT")] 

cupcake_filtered <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.filtered.abundance.txt", header = TRUE, as.is = T, sep = "\t") 

K17_cupcake_filtered <- abundance %>% filter(id %in% cupcake_filtered$pbid) %>% .[,c("K17_WT")] 
K17_cupcake_filtered_ig <- abundance %>% filter(id %ni% cupcake_filtered$pbid) %>% .[,c("K17_WT")] 

########## 3. sqanti classification.txt (pre-filtering)
class <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI/WholeIsoSeq.collapsed.filtered_classification.txt",header=T, as.is = T, sep = "\t")

K17_sqanti <- class[,c("isoform","FL.K17_WT")] 


########## 3. sqanti classification.txt (post-filtering)
class_filtered <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.txt",header=T, as.is = T, sep = "\t")
class_filtered_ig <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite_reasons.txt",header=T, as.is = T, sep = ",")


K17_sqanti_filtered <- class_filtered[,c("isoform","FL.K17_WT")]
K17_sqanti_filtered_ig <- class[class$isoform %in% class_filtered_ig$filtered_isoform,] %>% .[,c("isoform","FL.K17_WT")]

########## counts
cat("FL reads in read_stat file", nrow(K17_stats),"\n")
cat("FL reads in Cupcake collapsed abundance file", sum(K17_abundance),"\n")
cat("FL reads in Cupcake filtered abundance file", sum(K17_cupcake_filtered),"\n")
cat("FL reads in Sqanti file", sum(K17_sqanti$FL.K17_WT),"\n")
cat("FL reads in Sqanti filtered file", sum(K17_sqanti_filtered$FL.K17_WT),"\n")
cat("FL reads in Cupcake filtered removed file", sum(K17_cupcake_filtered_ig),"\n")
cat("FL reads in cupcake filtered remove file + cupcake filtered abundance file = read_stat file:",
    sum(K17_cupcake_filtered_ig) + sum(K17_cupcake_filtered) == nrow(K17_stats))
cat("FL reads in cupcake filtered remove file + sqanti + sqanti_filtered = read_stat file:",
    sum(K17_cupcake_filtered_ig) + sum(K17_sqanti_filtered$FL.K17_WT) + sum(K17_sqanti_filtered_ig$FL.K17_WT) == nrow(K17_stats))


