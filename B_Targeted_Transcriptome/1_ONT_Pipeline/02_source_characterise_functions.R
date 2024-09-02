## ---------- Script -----------------
##
## Purpose: project-related functions for characterising ONT mouse targeted transcriptome datasets
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ----------Plot colours-----------------

label_colour <- function(genotype){
  if(genotype == "WT"){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "TG"){colour = wes_palette("Royal1")[2]}else{
      if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
        if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
          if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
            if(genotype == "targeted"){colour = wes_palette("Darjeeling1")[2]}else{
              if(genotype == "whole"){colour = wes_palette("Darjeeling1")[1]}else{
                if(genotype == "whole+targeted"){colour = wes_palette("Darjeeling2")[1]}else{
                  if(genotype %in% c("isoseq","IsoSeq")){colour = wes_palette("Darjeeling2")[2]}else{
                    if(genotype == "rnaseq"){colour = wes_palette("Darjeeling2")[1]}else{
                      if(genotype == "Plus"){colour = wes_palette("GrandBudapest1")[2]}else{
                        if(genotype == "Minus"){colour = wes_palette("GrandBudapest1")[3]}else{
                          if(genotype == "Minus"){colour = wes_palette("GrandBudapest1")[3]}else{
                            if(genotype == c("ONT")){colour = wes_palette("Darjeeling1")[4]}else{
                              if(genotype == c("Yes")){colour = alpha(wes_palette("Cavalcanti1")[4],0.8)}else{
                                if(genotype == c("No")){colour = alpha(wes_palette("Cavalcanti1")[5],0.5)}else{
                                  if(genotype == c("BothTech")){colour = alpha(wes_palette("Moonrise3")[2],0.5)}else{
                                  }}}}}}}}}}}}}}}}}
  return(colour)
}


label_colour_cate <- function(genotype){
  if(genotype == "FSM"){colour =alpha("#00BFC4",0.8)}else{
    if(genotype == "ISM"){colour = alpha("#00BFC4",0.3)}else{
      if(genotype == "NIC"){colour = alpha("#F8766D",0.8)}else{
        if(genotype == "NNC"){colour = alpha("#F8766D",0.3)}else{
          if(genotype == "Genomic"){colour = alpha("#808080",0.3)}
        }}}}
  return(colour)
}


## ---------- Dataset defined functions -----------------

# Aim: find the mapt transgene sequencing after grep mapt human and mouse sequence in clustered.fasta
# while using clustered.fasta, plot shows the occurence of transgene in original raw reads 
# since use the clustered_report.csv description detailing the number of raw reads clustered to each transcript
# Input:
  # maptdir = str: directory path of files containing hmapt1_all_reads.csv, mmapt1_all_reads.csv, pre_cluster_read.csv
  # phenotype = df: phenotype data of samples <Sample.ID, Age_in_months, Genotype>
# Output:
  # p1: ratio of human and mouse transgene reads compared to all reads

find_mapt <- function(maptdir, phenotype){
  
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
    labs(x = "Age (months)", y = "Ratio of species-specific \n MAPT transcripts/ total transcripts", title = "ONT") +
    mytheme + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))  +
    theme(legend.position = c("bottom"), axis.text.y = element_text(angle = 90)) + 
    stat_summary(data = tallied, aes(x=Age, y=normalised, group=Phenotype), 
                 fun="mean", geom="line", linetype = "dotted") + 
    facet_wrap(~mapt_specific) + theme(strip.background = element_blank())
  
  return(p1)
}
