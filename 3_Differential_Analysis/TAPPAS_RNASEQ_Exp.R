library("stringr")
#rnaseq_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/ALLRNASEQ"
rnaseq_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/ALLRNASEQ"

abundance = lapply(list.files(path = rnaseq_dir, recursive = T, include.dirs = T, pattern = "abundance.tsv", full.names = T),
                   function(x) read.table(x,header = T))
samples = word(list.files(path = rnaseq_dir, recursive = T, include.dirs = T, pattern = "abundance.tsv"),c(1),sep=fixed("/"))

mod_abundance = lapply(abundance, function(x) x[,c("target_id","est_counts")])
names(mod_abundance) = samples
abundance %>% reduce(full_join, by = "V1")

all = dplyr::bind_rows(mod_abundance, .id = "sample")
counts = all %>% spread(sample, est_counts)

colnames(counts)[1] <- ""
#write.table(counts, paste0(rnaseq_dir,"/WholeMouseRNASeq_sqantitamafiltered.expression.txt"), quote = F, row.names = F,sep = "\t")


#### RNASeq to IsoSeq Targeted 
colnames(counts)[2:length(colnames(counts))] = paste0("FL.",colnames(counts)[2:length(colnames(counts))])
TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")
targeted.class.names.files = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt", as.is = T, sep = "\t", header = T)
TargetIsoforms <- targeted.class.names.files[toupper(targeted.class.names.files$associated_gene) %in% TargetGene,"isoform"]
counts = counts[counts[[1]] %in% TargetIsoforms,]
write.table(counts, paste0(rnaseq_dir,"/TargetedMouseRNASeq_sqantitamafiltered.expression.txt"), quote = F, row.names = F,sep = "\t")

