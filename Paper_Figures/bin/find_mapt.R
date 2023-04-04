# Aim: find the mapt transgene sequencing after grep mapt human and mouse sequence in clustered.fasta
# while using clustered.fasta, plot shows the occurence of transgene in original raw reads 
# since use the clustered_report.csv description detailing the number of raw reads clustered to each transcript
# Input:
# maptdir = str: directory path of files containing hmapt1_all_reads.csv, mmapt1_all_reads.csv, pre_cluster_read.csv
# phenotype = df: phenotype data of samples <Sample.ID, Age_in_months, Genotype>
# Output:
# p1: ratio of human and mouse transgene reads compared to all reads

find_mapt_isoseq <- function(maptdir, phenotype){
  
  # input 
  mapt_files <- list.files(path = maptdir, pattern = ".csv", full.names = T)
  mapt_files <- lapply(mapt_files, function(x) read.csv(x))
  names(mapt_files) <- list.files(path = maptdir, pattern = ".csv")
  mapt_files$pre_cluster_read.csv <- mapt_files$pre_cluster_read.csv %>% mutate(sample = word(file,c(1),sep=fixed(".")))
  
  species_specific <- function(all_reads, phenotype, cluster_reads, species){
    # all_reads format <matched_reads;sample>
    # matched_reads = transcript/X full_length_coverage=X;length=Y; want to extract X
    # X = number of raw CCS reads clustered to transcript
    all_reads <- all_reads %>% 
      mutate(readID = word(matched_reads, c(1), sep = ";")) %>% 
      mutate(matched_num_reads = as.numeric(word(readID, c(2), sep = "=")),
             readID = word(readID, c(1), sep = " "))
    
    # tally the number of matched reads per sample 
    tallied <- all_reads %>% group_by(sample) %>% tally(matched_num_reads)
    
    # assign phenotype data to samples
    tallied <- merge(tallied, phenotype, by = "sample", by.y = "Sample.ID", all = T)
    
    # tabulate with the total number of raw reads in each sample 
    tallied <- merge(tallied, cluster_reads, by = "sample", all = T)
    
    # some samples in phenotype have 0 matched reads, therefore no entry in all_reads file
    # therefore in merged dataset, replace NA in n (tally matched_num_reads) wih 0
    tallied <- tallied %>% mutate(n = ifelse(is.na(n), 0, n),
                                  normalised = n/num_reads, 
                                  mapt_specific = species)
    
    # calculate the mean across age and genotype
    mean_reads <- tallied %>% group_by(Age_in_months, Genotype) %>%  
      summarise_at(vars(normalised), funs(mean(., na.rm=TRUE))) %>% 
      mutate(mapt_specific = species) %>% as.data.frame()
    
    # remove samples with no phenotype
    tallied <- tallied %>% filter(!is.na(Genotype))
    
    output = list(all_reads, tallied, mean_reads)
    names(output) = c("all","tallied","mean")
    
    return(output)
  }
  
  humanMAPT <- species_specific(mapt_files[["hmapt1_all_reads.csv"]], phenotype, mapt_files[["pre_cluster_read.csv"]], "Human")
  mouseMAPT <- species_specific(mapt_files[["mmapt1_all_reads.csv"]], phenotype, mapt_files[["pre_cluster_read.csv"]], "Mouse")
  
  p1 <- bind_rows(humanMAPT$tallied,mouseMAPT$tallied) %>% 
    ggplot(., aes(x = Age_in_months, y = normalised, color = Genotype)) + geom_jitter(width = 0.09) + mytheme +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) +
    labs(x = "Age (months)", y = "Ratio of \nMAPT reads to total reads") +
    mytheme + #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))  +
    theme(legend.position = c("right"), axis.text.y = element_text(angle = 90)) + 
    stat_summary(data=bind_rows(humanMAPT$mean, mouseMAPT$mean), aes(x=Age_in_months, y=normalised, group=Genotype), 
                 fun.y="mean", geom="line", linetype = "dotted") + 
    facet_wrap(~mapt_specific) + theme(strip.background = element_blank())
  
  p1a <- p1 + theme(legend.position = "none")
  
  return(list(p1,p1a))
}


find_mapt_ont <- function(maptdir, phenotype){
  
  # input all_counts.csv 
  mapt_files <- list.files(path = maptdir, pattern = "all_counts.csv", full.names = T, recursive = TRUE)
  mapt_files <- lapply(mapt_files, function(x) read.csv(x))
  names(mapt_files) <- list.files(path = maptdir, pattern = "all_counts.csv", recursive = TRUE)
  
  # combine different datasets across human and mouse grepped reads
  all_counts <- data.frame(do.call(rbind, mapt_files)) %>% rownames_to_column(., var = "Dataset")
  
  # datawrnagle to generate barcodedSample for downstream merging with phenotype file
  # keep only output files grepped to hmapt1 and mmapt1
  all_counts <- all_counts %>% 
    mutate(Sample = word(sample,c(1), sep = fixed("_")),
           Batch =  word(Dataset,c(1), sep = fixed("/")),
           BarcodedSample = paste0(Batch,Sample),
           mapt_specific = word(word(Dataset,c(2),sep = fixed("/")),c(1),sep=fixed("_"))) %>% 
    filter(mapt_specific %in% c("hmapt1","mmapt1")) %>%
    mutate(mapt_specific = recode(mapt_specific, hmapt1 = 'Human', mmapt1 = 'Mouse'))
  
  # tally the number of matched reads (across human and mapt sequence) per sample 
  tallied <- all_counts %>% group_by(BarcodedSample, all_reads, mapt_specific) %>% tally(matched_reads)
  
  # assign phenotype data to samples
  tallied <- merge(tallied,phenotype,by = "BarcodedSample", all = T) %>%
    # note WT would not have human-specific MAPT reads therefore replace NA with 0
    mutate(matched_reads = replace_na(n, 0)) %>%
    # normalise with ratio
    mutate(normalised = n/all_reads)
  
  # remove samples with no phenotype
  tallied <- tallied %>% filter(!is.na(Phenotype))
  
  p1 <- ggplot(tallied, aes(x = Age, y = normalised, color = Phenotype)) + geom_jitter(width = 0.09) + mytheme +
    scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) +
    labs(x = "Age (months)", y = "Ratio of \nMAPT reads to total reads") +
    mytheme + #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))  +
    theme(legend.position = c("None"), axis.text.y = element_text(angle = 90)) + 
    stat_summary(data = tallied, aes(x=Age, y=normalised, group=Phenotype), 
                 fun="mean", geom="line", linetype = "dotted") + 
    facet_wrap(~mapt_specific) + theme(strip.background = element_blank()) 
  
  return(p1)
}
