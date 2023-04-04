## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: sources functions for generating downstream plots for SFARI dataset 
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
suppressMessages(library(ggrepel))
suppressMessages(library(forcats))
suppressMessages(library(extrafont))
suppressMessages(loadfonts())


## ----------Functions-----------------

# load all the functions
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/draw_venn.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/draw_density.R"))
source(paste0(LOGEN_ROOT, "compare_datasets/base_comparison.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_usage.R"))
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)

## ----------Theme-----------------

label_colour <- function(genotype){
  if(genotype %in% c("WT","Control","CONTROL")){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "WT_2mos"){colour = alpha(wes_palette("Royal1")[2],0.5)}else{
      if(genotype %in% c("TG","Case","CASE")){colour = wes_palette("Royal1")[2]}else{
        if(genotype == "TG_2mos"){colour = alpha(wes_palette("Royal1")[1],0.5)}else{
          if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
            if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
              if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
              }}}}}}}
  return(colour)
}

label_group <- function(genotype){
  if(genotype %in% c("TG","Case","CASE")){group = "TG"}else{
    if(genotype %in% c("WT","Control","CONTROL")){group = "WT"}}
  return(group)
}

vennONTvsIso <- function(classf){
  ont <- classf %>% filter(Dataset %in% c("ONT","Both")) %>% mutate(Dataset = "ONT") 
  isoseq <- classf  %>% filter(Dataset %in% c("Both","Iso-Seq")) %>% mutate(Dataset = "Iso-Seq") 
  
  p <- twovenndiagrams(ont$isoform,isoseq $isoform, "ONT", "Iso-Seq")
  return(p)
}

pSensitivity <- function(classf){
  
  dat <- classf %>% arrange(-nreads) %>% mutate(cumreads = cumsum(nreads), relative = prop.table(nreads), cumrel = cumsum(relative)) 
  p <- ggplot(dat, aes(x = isoform, y = cumrel, label = paste0(associated_gene,", ",associated_transcript))) + 
    geom_point() + 
    aes(x = fct_reorder(isoform, cumrel)) + 
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "XX") +
    geom_label_repel(data          = subset(dat, cumrel < 0.48),
                     size          = 4,
                     box.padding   = 0.5,
                     point.padding = 0.5,
                     force         = 100,
                     segment.size  = 0.2,
                     segment.color = "grey50",
                     direction     = "x") +
    labs(x = "Isoform", y = "Cumulative read proportion") + mytheme
  
  return(p)
  
  #densityfill <- function(x){
  #  if(x <= 10){return("<10")
  #  }else if (10 < x & x <= 100){return("10-100")
  #  }else if (100 < x & x < 200){return("100-200")
  #  }else {return(">200")}
  #}
  #class.files$targ_all$ndensity <- factor(unlist(lapply(class.files$targ_all$nreads, function(x) densityfill(x))),
  #                                        levels = c("<10","10-100","100-200",">200"))
  #ggplot(class.files$targ_all, aes(x = associated_gene, fill = as.factor(ndensity))) + 
  #  geom_bar() + labs(x = "Target genes", y = "Number of isoforms") + 
  #  scale_fill_discrete(name = "Number of total reads") + mytheme +
  #  theme(legend.position="bottom")
  
}
