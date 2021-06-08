# 22/09/2020: Validation of human-specific MAPT sequence in TG
# Prequisite: Tg4510_hMAPT.sh to count the number of human-specific MAPT sequences in all 12 samples (whole transcriptome only)

# 14/01/2021: Method 2 to test the genes affected by human MAPT

# 3 Methods of determining human-specific MAPT sequences in WT vs TG
# Method1: human-specific MAPT counts in CCS reads
  #1. Input: hMAPT2.header = two columns: V1: transcript name, V2: number of FL counts: length \n
  #2. The genes thought to be affected by transgene
  #3. Input2: cluster_report.csv --> for the number of CCS reads
# Method2:


# Plots:
# p1: Scatter plot of ratio of ccs reads/total ccs reads against age
# p2: WT and TG of all the human genes that were detected from alignment of mouse Iso-Seq dataset to human genome (hg38) by TG and WT
# p3: MAPT from SQANTI2 annotation post hg38 alignment
# p4: WT and TG of all the human genes that were detected from alignment of mouse Iso-Seq dataset to human genome (hg38) by age
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Rmarkdown_Input.R")

################################################################## Method 1
## Determine the number of Cluster reads with human-specific MAPT sequence for each sample
#hMAPT.header for each file contains multiple HQ, FL-polished transcripts (different transcript names e.g "@transcriptX", "@transcriptY"), but which all contain the same hMAPT sequence. Reason that there are multiple transcripts is due to collapsed properly (redundancy) from Iso-Seq3. For this reason, the count of human-specific MAPT in each sample is calculated by the sum of FL counts for all these multiple transcripts.

########### Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
<<<<<<< HEAD:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
hMAPT_input <- list.files(path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/QC/Human_Mapt/Sequences/Whole/Individual_Clustered", pattern = "hMAPT2.header", full.names = T)
cluster_reads <- list.files("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/CLUSTER", pattern = "cluster_report.csv", full.names = T)
=======
hMAPT_input <- list.files(path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER/hMAPT", pattern = "hMAPT2.header", full.names = T)
ccs_reads <- list.files("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER", pattern = "cluster_report.csv", full.names = T)
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R

# only input hMAPT2.header files with data entry i.e. from TG mice
# Note: files of WT mice should be empty i.e file.size == 0 given no human-specific MAPT sequence
# first save file path into file_input_names, then read table and assign name to table based on filename
for (file in hMAPT_input){
  if (file.size(file) == 0) next
  put(file)
}
file_input_names <- magic_result_as_vector()
file_input <- lapply(file_input_names, read.table)
names(file_input) <- lapply(file_input_names, function(x) word(x,c(13),  sep = fixed ('/')))

# Tabulate the number of CCS reads for each @transcriptX/@transcriptY and then sum across all transcripts per sample (as above, same transcript even though named diferently due to redundancy of collapse)
# split V2 from hMAPT.header to extract the number of CCS reads
<<<<<<< HEAD:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
file_input <- lapply(file_input, function(x) x %>% mutate(cluster_reads = word(V2, c(1), sep = ";")) %>% mutate(cluster_reads = as.numeric(word(cluster_reads, c(2), sep = "="))))
=======
file_input <- lapply(file_input, function(x) x %>% mutate(CCS_reads = word(V2, c(1), sep = ";")) %>% mutate(CCS_reads = as.numeric(word(CCS_reads, c(2), sep = "="))))
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
# for each file (/sample), sum the number of CCS reads across all transcripts
for(i in 1:length(file_input)){
  sample <- names(file_input)[i]
  sum_cluster_reads <- sum(file_input[[i]]$cluster_reads)
  put(sample, sum_cluster_reads)
}
hMAPT_reads <- magic_result_as_dataframe() %>% mutate(sample = word(sample, c(1), sep = fixed(".")))


## Determine the number of CCS reads with human-specific MAPT sequence for each sample
########### Tally the total number of successfully generated CCS reads per sample for normalisation
# read in cluster_report.csv from all 12 samples, combine as one big dataframe and use "sample" identifier from file name
<<<<<<< HEAD:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
cluster_reads_input <- lapply(cluster_reads, function(x) read.table(x, header = T))
names(cluster_reads_input) <- lapply(cluster_reads, function(x) word(x,c(12),  sep = fixed ('/')))
cluster_reads_input_all <- do.call(rbind, cluster_reads_input)
cluster_reads_input_all <- setDT(cluster_reads_input_all, keep.rownames = TRUE)[] %>% mutate(sample = word(rn, c(1), sep = fixed(".")))

# tally the number of CCS reads per sample
total_cluster_reads <- cluster_reads_input_all %>% group_by(sample) %>% tally() %>% as.data.frame()
=======
ccs_reads_input <- lapply(ccs_reads, function(x) read.table(x, header = T))
names(ccs_reads_input) <- lapply(ccs_reads, function(x) word(x,c(12),  sep = fixed ('/')))
ccs_reads_input_all <- do.call(rbind, ccs_reads_input)
ccs_reads_input_all <- setDT(ccs_reads_input_all, keep.rownames = TRUE)[] %>% mutate(sample = word(rn, c(1), sep = fixed(".")))

# tally the number of CCS reads per sample
total_ccs_reads <- ccs_reads_input_all %>% group_by(sample) %>% tally() %>% as.data.frame()
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R


########### Plots
# Classifiers
twomos <- c("K18","O18","S18", "K17","M21","Q21")
WT <- c("K17","M21","Q21","K23","O23","S23")

# Prepare plot by merging counts of human-specific MAPT reads and total reads
<<<<<<< HEAD:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
plot_reads <- merge(hMAPT_reads, total_cluster_reads, by = "sample", all = T) %>% select(-i) %>%
=======
plot_reads <- merge(hMAPT_reads, total_ccs_reads, by = "sample", all = T) %>% select(-i) %>%
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
  # do not include other J20 samples
  filter(!sample %in% c("C21","E18","C20","B21")) %>%
  # classifiers of Age and Genotype
  mutate(Age = ifelse(grepl(paste(twomos,collapse="|"), sample),"2","8")) %>%
  mutate(Genotype = ifelse(grepl(paste(WT,collapse="|"), sample),"WT","TG")) %>%
  # note WT would not have human-specifici MAPT reads therefore replace NA with 0
<<<<<<< HEAD:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
  mutate(sum_cluster_reads = replace_na(sum_cluster_reads, 0)) %>%
=======
  mutate(sum_ccs_reads = replace_na(sum_ccs_reads, 0)) %>%
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:1_Transcriptome_Annotation/Isoseq3_QC/Tg4510_hMAPT.R
  # normalise with ratio
  mutate(normalised = sum_cluster_reads/n)

# p1: Scatter plot of ratio of ccs reads/total ccs reads against age
p1 <- ggplot(plot_reads, aes(x = Age, y = normalised, color = Genotype)) + geom_jitter(width = 0.09) + mytheme +
  scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) +
  labs(x = "Age (months)", y = "Ratio CCS reads/ total CCS reads") +
  mytheme

################################################################## Method 2
## The genes thought to be affected by transgene
# demultiplexed after merging all files
# cupcake collapse using lowered threshold (95%) for coverage and, applying TAMA filter after SQANTI filter
Sqanti_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER/"
Sqanti_input_file <- read.table(paste0(Sqanti_input_dir,"WholeIsoSeq.collapsed_sqantitamafiltered.classification.txt"), header = T)

# Genes affected by transgene insertion
#Vipr2 (vasoactive intestinal peptide receptor 2),
#Wdr60 (WD repeat-containing protein 60),
#Esyt2 (extended synaptotagmin-like protein 2),
#Ncapg2 (non-SMC condensin II complex, subunit G2),
#Ptprn2 (protein tyrosine phosphatase, receptor type, N polypeptide 2)

affected_genes <- c("Vipr2","Wdr60","Esyt2","Ncapg2","Ptprn2")
WT <- c("FL.K17_WT","FL.K23_WT","FL.M21_WT","FL.O23_WT","FL.Q21_WT","FL.S23_WT")
TG <- c("FL.K18_TG","FL.K24_TG","FL.L22_TG","FL.O18_TG","FL.Q20_TG","FL.S18_TG")
Sqanti_input_file %>% filter(associated_gene %in% affected_genes) %>% select(structural_category,associated_gene,FL.K17_WT:FL.S23_WT) %>%
  mutate(WT_reads = row)
  melt() %>% mutate(phenotype = word(.$variable, c(2), sep = fixed("_"))) %>%
  ggplot(., aes(x = structural_category, y = value, group = phenotype, colour = phenotype)) + geom_point() +
  facet_grid(~associated_gene)

################################################################## Method 3
## Determine the number of FL counts mapped to human MAPT after hg38 genome alignment and SQANTI2 annotation

# Read in classification files from SQANTI2 and bind to one big dataframe, with "sample" as identifier
Sqanti_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TO_HUMAN/SQANTI2_v7"
Sqanti_input_file <- list.files(path = Sqanti_input_dir, pattern = "filtered_lite_classification.txt", full.names = T)
Sqanti <- lapply(Sqanti_input_file, function(x) read.table(x, header = T))
names(Sqanti) <- lapply(Sqanti_input_file, function(x) word(x,c(12),  sep = fixed ('/')))
Sqanti_all <- do.call(rbind, Sqanti)
Sqanti_all  <- setDT(Sqanti_all, keep.rownames = TRUE)[] %>% mutate(Sample = word(rn, c(1), sep = fixed(".")))

# Further classify the input from SQANTI2 based on age and genotype, and remove J20 input
Sqanti_all <- Sqanti_all %>%
  mutate(Age = ifelse(grepl(paste(twomos,collapse="|"), Sample),"2","8")) %>%
  mutate(Genotype = ifelse(grepl(paste(WT,collapse="|"), Sample),"WT","TG")) %>%
  filter(!Sample %in% c("C21","E18","C20","B21"))

########### Plots
# p2: WT and TG of all the human genes that were detected from alignment of mouse Iso-Seq dataset to human genome (hg38) by TG and WT
p2 <- Sqanti_all %>%
  group_by(associated_gene, Sample, Age, Genotype) %>% tally(FL) %>%
  ggplot(., aes(x = reorder(associated_gene, n), y = n, color = Genotype)) + geom_boxplot() + geom_jitter() + coord_flip() + mytheme +
  labs(y = "FL Reads", x = "Gene (hg38)") + scale_color_manual(values = c(label_colour("TG"),label_colour("WT")))

# p3: MAPT from SQANTI2 annotation post hg38 alignment
# Note missing two samples?
p3 <- Sqanti_all %>%
  # tally the number of FL reads per detected human gene and filter for MAPT
  group_by(associated_gene, Sample, Age, Genotype) %>% tally(FL) %>% filter(associated_gene == "MAPT") %>%
  ggplot(., aes(x = Age, y = n, color = Genotype, label = Sample)) + geom_jitter(width = 0.085) + mytheme + geom_text() +
  scale_color_manual(values = c(label_colour("TG"))) +
  labs(x = "Age (months)", y = "FL Reads mappted to MAPT") + scale_y_continuous(breaks = seq(0, 12, by = 2)) + mytheme +
  theme(legend.position = "bottom")

# p4: WT and TG of all the human genes that were detected from alignment of mouse Iso-Seq dataset to human genome (hg38) by age
p4 <- Sqanti_all %>% group_by(associated_gene, Sample, Age) %>% tally(FL) %>%
  ggplot(., aes(x = reorder(associated_gene, n), y = n, color = Age)) + geom_boxplot() + geom_jitter() + coord_flip() + mytheme +
  labs(y = "FL Reads", x = "Gene (hg38)") + mytheme

################################################################## Output
pdf(paste0(output_plot_dir,"/Tg4510_hMAPT.pdf"), width = 11, height = 8.5)
p1
p2
p3
p4
dev.off()
