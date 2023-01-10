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

find_mapt <- function(){
  # read in cluster_report.csv from all 24 samples, combine as one big dataframe and use "sample" identifier from file name
  cluster_reads_input <- lapply(cluster_reads, function(x) read.table(x, header = T))
  names(cluster_reads_input) <- lapply(cluster_reads, function(x) word(x,c(12),  sep = fixed ('/')))
  cluster_reads_input_all <- do.call(rbind, cluster_reads_input)
  cluster_reads_input_all <- setDT(cluster_reads_input_all, keep.rownames = TRUE)[] %>% mutate(sample = word(rn, c(1), sep = fixed(".")))
  # tally the number of CCS reads per sample
  total_cluster_reads <- cluster_reads_input_all %>% group_by(sample) %>% tally() %>% as.data.frame()
  
  species_specific_mapt <- function(MAPT_input, species){
    MAPT_reads <- merge(MAPT_input,targetedpheno, by = "Sample")
    
    plot_reads <- merge(MAPT_reads, total_cluster_reads, by.x = "Sample", by.y = "sample", all.x = T) %>%
      # normalise with ratio
      mutate(normalised = Counts/n, mapt_specific = species)
    
    mean_reads = plot_reads %>% group_by(Age,Phenotype) %>%  summarise_at(vars(normalised), funs(mean(., na.rm=TRUE))) %>% 
      mutate(mapt_specific = species)
    
    return(list(plot_reads, as.data.frame(mean_reads)))
  }
  
  final_humanMAPT <- species_specific_mapt(hMAPT_input,"Human")
  final_mouseMAPT <- species_specific_mapt(mMAPT_input,"Mouse")
  
  p <- ggplot(final_humanMAPT[[1]], aes(x = Age, y = normalised, color = Phenotype)) + geom_jitter(width = 0.09) + mytheme +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) +
    labs(x = "Age (months)", y = "Ratio of species-specific \n MAPT transcripts/ total transcripts") +
    mytheme + scale_y_continuous(lim = c(0,0.006), labels = function(x) format(x, scientific = TRUE))  +
    theme(legend.position = c("bottom"), axis.text.y = element_text(angle = 90)) + 
    stat_summary(data=final_humanMAPT[[2]], aes(x=Age, y=normalised, group=Phenotype), 
                 fun.y="mean", geom="line", linetype = "dotted") 
  
  return(p)
}