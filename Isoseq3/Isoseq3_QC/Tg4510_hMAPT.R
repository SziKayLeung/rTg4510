source("/gpfs/ts0/home/sl693/Scripts/IsoSeq3_Tg4510/Isoseq3/Isoseq3_QC/Rmarkdown_Input.R")

hMAPT_input <- list.files(path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER/hMAPT", pattern = "hMAPT2.header", full.names = T)
ccs_reads <- list.files("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER", pattern = "cluster_report.csv", full.names = T)

ccs_reads_input <- lapply(ccs_reads, function(x) read.table(x, header = T))
names(ccs_reads_input) <- lapply(ccs_reads, function(x) word(x,c(12),  sep = fixed ('/')))
ccs_reads_input_all <- do.call(rbind, ccs_reads_input)
ccs_reads_input_all <- setDT(ccs_reads_input_all, keep.rownames = TRUE)[] %>% mutate(sample = word(rn, c(1), sep = fixed(".")))

total_ccs_reads <- ccs_reads_input_all %>% group_by(sample) %>% tally() %>% as.data.frame()


# only input files with data entry i.e TG
for (file in hMAPT_input){
  if (file.size(file) == 0) next
  put(file)
}

file_input_names <- magic_result_as_vector() 
file_input <- lapply(file_input_names, read.table)
names(file_input) <- lapply(file_input_names, function(x) word(x,c(13),  sep = fixed ('/')))

file_input <- lapply(file_input, function(x) x %>% mutate(CCS_reads = word(V2, c(1), sep = ";")) %>% mutate(CCS_reads = as.numeric(word(CCS_reads, c(2), sep = "=")))) 
for(i in 1:length(file_input)){
  sample <- names(file_input)[i]
  sum_ccs_reads <- sum(file_input[[i]]$CCS_reads)
  put(sample, sum_ccs_reads)
}

twomos <- c("K18","O18","S18", "K17","M21","Q21")
WT <- c("K17","M21","Q21","K23","O23","S23")
hMAPT_reads <- magic_result_as_dataframe() %>% mutate(sample = word(sample, c(1), sep = fixed(".")))

plot_reads <- merge(hMAPT_reads, total_ccs_reads, by = "sample", all = T) %>% select(-i) %>% 
  filter(!sample %in% c("C21","E18","C20","B21")) %>%
  mutate(Age = ifelse(grepl(paste(twomos,collapse="|"), sample),"2","8")) %>%
  mutate(Genotype = ifelse(grepl(paste(WT,collapse="|"), sample),"WT","TG")) %>%
  mutate(sum_ccs_reads = replace_na(sum_ccs_reads, 0)) %>% 
  mutate(normalised = sum_ccs_reads/n)

ggplot(plot_reads, aes(x = Genotype, y = normalised, fill = Age)) + geom_bar(stat = "identity", position = position_dodge()) + mytheme + 
  scale_fill_manual(values = c(label_colour("TG"),label_colour("WT"))) + 
  labs(x = "Age (months)", y = "Ratio ccs reads/ total ccs reads")

ggplot(plot_reads, aes(x = Age, y = normalised, color = Genotype)) + geom_jitter(width = 0.09) + mytheme + 
  scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) + 
  labs(x = "Age (months)", y = "Ratio ccs reads/ total ccs reads")
  
################ MAPT after SQANTI2
Sqanti_input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TO_HUMAN/SQANTI2_v7"
Sqanti_input_file <- list.files(path = Sqanti_input_dir, pattern = "filtered_lite_classification.txt", full.names = T)
Sqanti <- lapply(Sqanti_input_file, function(x) read.table(x, header = T))
names(Sqanti) <- lapply(Sqanti_input_file, function(x) word(x,c(12),  sep = fixed ('/')))

Sqanti_all <- do.call(rbind, Sqanti)
Sqanti_all  <- setDT(Sqanti_all, keep.rownames = TRUE)[] %>% mutate(Sample = word(rn, c(1), sep = fixed(".")))

twomos <- c("K18","O18","S18", "K17","M21","Q21")
WT <- c("K17","M21","Q21","K23","O23","S23")

Sqanti_all <- Sqanti_all %>%  
  mutate(Age = ifelse(grepl(paste(twomos,collapse="|"), Sample),"2","8")) %>%
  mutate(Genotype = ifelse(grepl(paste(WT,collapse="|"), Sample),"WT","TG")) %>%
  filter(!Sample %in% c("C21","E18","C20","B21")) 


# WT vs TG 
Sqanti_all %>% group_by(associated_gene, Sample, Age, Genotype) %>% tally(FL) %>%
  ggplot(., aes(x = reorder(associated_gene, n), y = n, color = Genotype)) + geom_boxplot() + geom_jitter() + coord_flip() + mytheme + 
  labs(y = "FL Reads", x = "Gene (hg38)") + scale_color_manual(values = c(label_colour("TG"),label_colour("WT")))

# MAPT specifically 
Sqanti_all %>% group_by(associated_gene, Sample, Age, Genotype) %>% tally(FL) %>% filter(associated_gene == "Mapt") %>% 
  ggplot(., aes(x = Age, y = n, color = Genotype, label = Sample)) + geom_jitter(width = 0.085) + mytheme + geom_text() +
  scale_color_manual(values = c(label_colour("TG"))) + 
  labs(x = "Age (months)", y = "FL Reads mappted to MAPT") + scale_y_continuous(breaks = seq(0, 12, by = 2))

# Age
Sqanti_all %>% group_by(associated_gene, Sample) %>% tally(FL) %>% 
  mutate(Age = ifelse(grepl(paste(twomos,collapse="|"), Sample),"2","8")) %>%
  filter(!Sample %in% c("C21","E18","C20","B21")) %>%
  mutate(Genotype = ifelse(grepl(paste(WT,collapse="|"), Sample),"WT","TG")) %>%
  ggplot(., aes(x = reorder(associated_gene, n), y = n, color = Age)) + geom_boxplot() + geom_jitter() + coord_flip() + mytheme + 
  labs(y = "FL Reads", x = "Gene (hg38)")
