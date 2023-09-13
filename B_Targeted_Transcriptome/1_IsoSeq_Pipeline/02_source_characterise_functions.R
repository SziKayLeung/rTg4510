## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
## 
##   
##
## 

## ---------- Packages -----------------

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(rjson)) # json files
suppressMessages(library(plyr)) # revalue
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(stringr)) 
suppressMessages(library(viridis)) 
suppressMessages(library(wesanderson)) 
suppressMessages(library(extrafont))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(tibble))
suppressMessages(library(VennDiagram))
suppressMessages(library(directlabels))
suppressMessages(library(cowplot))
suppressMessages(library(readxl))
suppressMessages(library(ggdendro))
suppressMessages(library(pheatmap))
suppressMessages(library(extrafont))
suppressMessages(loadfonts())


## ----------Plot colours-----------------

# plot label colour
label_colour <- function(genotype){
  if(genotype %in% c("Control","WT")){colour = wes_palette("Royal1")[1]}else{
    if(genotype %in% c("Case", "TG")){colour = wes_palette("Royal1")[2]}else{
      if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
        if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
          if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
            if(genotype == "targeted"){colour = wes_palette("Darjeeling1")[2]}else{
              if(genotype == "whole"){colour = wes_palette("Darjeeling1")[1]}else{
                if(genotype == "whole+targeted"){colour = wes_palette("Darjeeling2")[1]}else{
                  if(genotype == "isoseq"){colour = wes_palette("Darjeeling2")[2]}else{
                    if(genotype == "rnaseq"){colour = wes_palette("Darjeeling2")[1]}else{
                    }}}}}}}}}}
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
    
    output = list(all_reads, tallied, mean_reads)
    names(output) = c("all","tallied","mean")
    
    return(output)
  }
  
  humanMAPT <- species_specific(mapt_files[["hmapt1_all_reads.csv"]], phenotype, mapt_files[["pre_cluster_read.csv"]], "Human")
  mouseMAPT <- species_specific(mapt_files[["mmapt1_all_reads.csv"]], phenotype, mapt_files[["pre_cluster_read.csv"]], "Mouse")
  
  p1 <- bind_rows(humanMAPT$tallied,mouseMAPT$tallied) %>% 
    filter(Genotype != c("WT","TG")) %>%
    ggplot(., aes(x = Age_in_months, y = normalised, color = Genotype)) + geom_jitter(width = 0.09) + mytheme +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) +
    labs(x = "Age (months)", y = "Ratio of species-specific \n MAPT reads/ total reads") +
    mytheme + 
    theme(legend.position = c("bottom"), axis.text.y = element_text(angle = 90)) + 
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    stat_summary(data=bind_rows(humanMAPT$mean, mouseMAPT$mean), aes(x=Age_in_months, y=normalised, group=Genotype), 
                 fun.y="mean", geom="line", linetype = "dotted") + 
    facet_wrap(~mapt_specific) + theme(strip.background = element_blank())
  
  return(p1)
}
