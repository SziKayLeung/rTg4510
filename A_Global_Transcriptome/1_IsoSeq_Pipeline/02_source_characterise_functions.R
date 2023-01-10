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
suppressMessages(library(data.table))
suppressMessages(loadfonts())

## ----------Functions-----------------

# load all the functions
source_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise"
file.sources = list.files(path = source_dir, pattern="*.R", full = T)
sapply(file.sources,source,.GlobalEnv)

## ----------Plot colours-----------------

# plot label colour
label_colour <- function(genotype){
  if(genotype == "WT"){colour = wes_palette("Royal1")[2]}else{
    if(genotype == "WT_2mos"){colour = alpha(wes_palette("Royal1")[2],0.5)}else{
      if(genotype == "TG"){colour = wes_palette("Royal1")[1]}else{
        if(genotype == "TG_2mos"){colour = alpha(wes_palette("Royal1")[1],0.5)}else{
          if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
            if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
              if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
              }}}}}}}
  return(colour)
}

label_name <- function(variable){
  if(variable == "control"){name = "WT"}else{
    if(variable == "case"){name = "TG"}
  }
  return(name)
}



### Human MAPT #################################################################
## Determine the number of Cluster reads with human-specific MAPT sequence for each sample
#hMAPT.header for each file contains multiple HQ, FL-polished transcripts (different transcript names e.g "@transcriptX", "@transcriptY"), but which all contain the same hMAPT sequence. Reason that there are multiple transcripts is due to collapsed properly (redundancy) from Iso-Seq3. For this reason, the count of human-specific MAPT in each sample is calculated by the sum of FL counts for all these multiple transcripts.

########### Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
find_mapt <- function(){
  
  twomos <- c("K18","O18","S18", "K17","M21","Q21")
  WT <- c("K17","M21","Q21","K23","O23","S23")
  
  species_specific_mapt <- function(MAPT_input){
    # only input MAPT2.header files with data entry i.e. from TG mice
    # Note: files of WT mice for hMAPT should be empty i.e file.size == 0 given no human-specific MAPT sequence
    # first save file path into file_input_names, then read table and assign name to table based on filename
    count = 1
    file_input_names <- c()
    for (file in MAPT_input){
      if (file.size(file) == 0) next
      file_input_names[[count]] <- file
      count = count + 1
    }
    file_input <- lapply(file_input_names, read.table)
    names(file_input) <- lapply(file_input_names, function(x) word(word(x,c(12),  sep = fixed ('/')),c(1), sep = fixed(".")))
    
    # Tabulate the number of CCS reads for each @transcriptX/@transcriptY and then sum across all transcripts per sample (as above, same transcript even though named diferently due to redundancy of collapse)
    # split V2 from MAPT.header to extract the number of CCS reads
    file_input <- lapply(file_input, function(x) x %>% mutate(cluster_reads = word(V2, c(1), sep = ";")) %>% mutate(cluster_reads = as.numeric(word(cluster_reads, c(2), sep = "="))))
    # for each file (/sample), sum the number of CCS reads across all transcripts
    MAPT_reads <- data.frame()
    count = 1
    for(i in 1:length(file_input)){
      sample <- names(file_input)[i]
      sum_cluster_reads <- sum(file_input[[i]]$cluster_reads)
      MAPT_reads[count,1] <- sample
      MAPT_reads[count,2] <- sum_cluster_reads
      count = count + 1
    }
    colnames(MAPT_reads) <- c("sample","sum_cluster_reads")
    return(MAPT_reads)
  }
  
  hMAPT_reads <- species_specific_mapt(hMAPT_input)
  mMAPT_reads <- species_specific_mapt(mMAPT_input)

  # read in cluster_report.csv from all 12 samples, combine as one big dataframe and use "sample" identifier from file name
  cluster_reads_input <- lapply(cluster_reads, function(x) read.table(x, header = T))
  names(cluster_reads_input) <- lapply(cluster_reads, function(x) word(x,c(12),  sep = fixed ('/')))
  cluster_reads_input_all <- do.call(rbind, cluster_reads_input)
  cluster_reads_input_all <- setDT(cluster_reads_input_all, keep.rownames = TRUE)[] %>% mutate(sample = word(rn, c(1), sep = fixed(".")))
  
  # tally the number of CCS reads per sample
  total_cluster_reads <- cluster_reads_input_all %>% group_by(sample) %>% tally() %>% as.data.frame()
  
  # Prepare plot by merging counts of human-specific MAPT reads and total reads
  plot_species_MAPT <- function(MAPT_reads,species){
    plot_reads <- merge(MAPT_reads, total_cluster_reads, by = "sample", all = T) %>%
      # do not include other J20 samples
      filter(!sample %in% c("C21","E18","C20","B21")) %>%
      # classifiers of Age and Genotype
      mutate(Age = ifelse(grepl(paste(twomos,collapse="|"), sample),"2","8")) %>%
      mutate(Genotype = ifelse(grepl(paste(WT,collapse="|"), sample),"WT","TG")) %>%
      # note WT would not have human-specifici MAPT reads therefore replace NA with 0
      mutate(sum_cluster_reads = replace_na(sum_cluster_reads, 0)) %>%
      # normalise with ratio
      mutate(normalised = sum_cluster_reads/n, mapt_specific = species)
    
    mean_reads = plot_reads %>% group_by(Age,Genotype) %>%  summarise_at(vars(normalised), funs(mean(., na.rm=TRUE))) %>% 
      mutate(mapt_specific = species)
    
    return(list(plot_reads, as.data.frame(mean_reads)))
  }
  
   
  final_humanMAPT <- plot_species_MAPT(hMAPT_reads,"Human")
  final_mouseMAPT <- plot_species_MAPT(mMAPT_reads,"Mouse")
  
  p1 <- bind_rows(final_humanMAPT[[1]],final_mouseMAPT[[1]]) %>% 
   ggplot(., aes(x = Age, y = normalised, color = Genotype)) + geom_jitter(width = 0.09) + mytheme +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) +
    labs(x = "Age (months)", y = "Ratio of species-specific \n MAPT transcripts/ total transcripts") +
    mytheme + scale_y_continuous(lim = c(0,0.006), labels = function(x) format(x, scientific = TRUE))  +
    theme(legend.position = c("bottom"), axis.text.y = element_text(angle = 90)) + 
    stat_summary(data=bind_rows(final_humanMAPT[[2]],final_mouseMAPT[[2]]), aes(x=Age, y=normalised, group=Genotype), 
                 fun.y="mean", geom="line", linetype = "dotted") + 
    facet_wrap(~mapt_specific) + theme(strip.background = element_blank())
  
  
  # aligment identity and length of hMAPT reads from all clustered
  p2 <- map %>% 
    filter(name_of_read %in% hMAPT_input_all$V1) %>%
    mutate(identity =alignment_identity * 100) %>% 
    mutate(length = alignment_length_perc * 100) %>%
    ggplot(., aes(x = length, y = identity)) +
    geom_point() + mytheme + labs(x = "Alignment Length (%)", y = "Alignment Identity (%)") +
    geom_hline(yintercept = 85, linetype = 2) +
    geom_vline(xintercept = 95, linetype = 2) +
    annotate("text", x = 11, y = 82, label = "85% Threshold") +
    annotate("text", x = 92, y = 50, label = "95% Threshold", angle = 90) +
    annotate("rect", xmin = 0, xmax = 75, ymin = 95, ymax = 101, fill = wes_palette("Zissou1")[3], alpha = 0.2) +
    annotate("rect", xmin = 25, xmax = 75, ymin = 25, ymax = 70, fill = wes_palette("Zissou1")[1], alpha = 0.2) +
    annotate("rect", xmin = 95, xmax = 100, ymin = 85, ymax = 100, fill = "palegreen", alpha = 0.2) 
  
  #View(map %>% filter(name_of_read == c("transcript/6355")))
  
  return(list(p1,p2))
}

