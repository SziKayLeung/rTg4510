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
suppressMessages(library(ggtranscript))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer))


## ----------Functions-----------------

# load all the functions
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/draw_venn.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/draw_density.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/plot_basic_stats.R"))
source(paste0(LOGEN_ROOT, "compare_datasets/base_comparison.R"))
source(paste0(LOGEN_ROOT, "compare_datasets/whole_vs_targeted.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_usage.R"))
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)

## ----------Theme-----------------

label_colour <- function(genotype){
  if(genotype %in% c("WT","Control","CONTROL")){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "WT_2mos"){colour = alpha(wes_palette("Royal1")[2],0.5)}else{
      if(genotype %in% c("TG","Case","CASE")){colour =wes_palette("IsleofDogs1")[4]}else{
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


reportStats <- function(res,stats,isoList){
  for(i in isoList){
    print(res %>% select(!c("lfcSE","stat")) %>% filter(isoform == i))  
    cat("Wald statistic: ")
    cat(stats[rownames(stats) == i,"WaldStatistic_group_CASE_vs_CONTROL"],"\n")
    cat("******************\n")
  }
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


#' @title draw_heatmap_gene_level
#' @description This function uses hierarchal clustering to plot a heatmap of the gene expression across samples
#' @param diff_genes: tappAS output of list of differentially expressed genes 
#' @param genelevel_exp: tappAS output of normalised gene expression values
#' @param type: "glob_isoseq" used to determine annnotation of colours
#' @param type: "yes" or "no" whether to filter differentially expressed genes 

draw_heatmap_gene_level <- function(diff_genes,genelevel_exp, type, diff){
  
  dat <- genelevel_exp[,c("associated_gene", "Exp","time","variable")]
  
  if(diff == "yes"){
    dat <- dat %>% filter(associated_gene %in% rownames(diff_genes))
  }
  
  dat <- aggregate(Exp ~ time + variable + associated_gene, data = dat, mean) %>% 
    mutate(Exp = log2(Exp)) %>% 
    select(variable, associated_gene, Exp) %>% 
    spread(., variable, Exp) %>% tibble::column_to_rownames(var = "associated_gene")  
  
  # remove isoforms that have been removed by tappAS due to very low count 
  # replace infinity value from log value with 0 
  # rotate the dataframe for visualisation ease
  dat <- dat[,colSums(is.na(dat))<nrow(dat)]
  dat[dat == "-Inf"] <- 0
  
  # set the order for the column (Age, Genotype)
  coldata = genelevel_exp[,c("sample", "time", "group")] %>% 
    distinct(.keep_all = TRUE) %>% 
    column_to_rownames(var = "sample") %>% 
    mutate(time = as.factor(time)) %>% 
    mutate(group = ifelse(group == "CONTROL","WT","TG"))
  colnames(coldata) = c("Age (months)","Genotype")
  
  if(type == "glob_isoseq"){
    annotation_colors = list(
      Genotype = c(WT=wes_palette("Royal1")[1], TG=wes_palette("Royal1")[2]),
      `Age (months)` = c("2"="white", "8"="black"))
  }else{
    annotation_colors = list(
      Genotype = c(WT=wes_palette("Royal1")[1], TG=wes_palette("Royal1")[2]),
      `Age (months)` = c("2"="white","4"="#CFCFCF","6"="#777777","8"="black"))
  }
  
  if(diff == "yes"){
    p = pheatmap(dat, annotation_col=coldata, annotation_legend = TRUE,annotation_names_col = FALSE,
                 show_colnames = FALSE,show_rownames = TRUE, color = viridis(10),annotation_colors = annotation_colors,
                 fontsize_col = 20)
    
  }else{
    p = pheatmap(dat, annotation_col=coldata, annotation_legend = TRUE,annotation_names_col = FALSE,
                 show_colnames = FALSE,show_rownames = FALSE, color = viridis(10),annotation_colors = annotation_colors,
                 fontsize_col = 20)
  }
  
  return(p)
}

# Draw the heatmap for the Isoform expression of the gene 

draw_heatmap_gene <- function(gene, cf, normCounts, type){
  
  # Subset the normalised expression count to gene, and datawrangle for plot
  dat = normCounts[normCounts$associated_gene == gene,c("isoform","normalised_counts","sample")] %>%
    mutate(log2normalised = log2(normalised_counts)) %>% 
    select(isoform, log2normalised, sample) %>%
    spread(., isoform, log2normalised) %>% tibble::column_to_rownames(var = "sample") 
  
  # remove isoforms that have been removed by tappAS due to very low count 
  # replace infinity value from log value with 0 
  # rotate the dataframe for visualisation ease
  dat <- dat[,colSums(is.na(dat))<nrow(dat)]
  dat[dat == "-Inf"] <- 0
  dat.t <- t(dat)
  
  # set the order for the column (Age, Genotype)
  coldata = normCounts %>% 
    dplyr::select(sample, time, group) %>% distinct(.keep_all = TRUE) %>% column_to_rownames(var = "sample") %>% 
    mutate(time = as.factor(time))
  colnames(coldata) = c("Age (months)","Genotype")
  
  
  # set the order for the row (isoform structural category)
  rowdata = cf[cf$isoform %in% colnames(dat),c("isoform","structural_category")] %>% dplyr::select(-isoform)
  colnames(rowdata) = c("Category")
  
  # set annotation colours
  if(type == "targeted"){
    annotation_colors = list(
      Genotype = c(CONTROL=wes_palette("Royal1")[1], CASE=wes_palette("Royal1")[2]),
      `Age (months)` = c("2"="white","4"="#CFCFCF","6"="#777777","8"="black"), 
      Category = c(FSM = alpha("#00BFC4",0.8),ISM = alpha("#00BFC4",0.3),NIC = alpha("#F8766D",0.8),NNC = alpha("#F8766D",0.3)))
  }else{
    annotation_colors = list(
      Genotype = c(WT=wes_palette("Royal1")[1], TG=wes_palette("Royal1")[2]),
      `Age (months)` = c("2"="white", "8"="black"), 
      Category = c(FSM = alpha("#00BFC4",0.8),ISM = alpha("#00BFC4",0.3),NIC = alpha("#F8766D",0.8),NNC = alpha("#F8766D",0.3)))
  }
  
  
  # draw heatmap
  if(nrow(dat.t) > 1){
    p <- pheatmap(dat.t, 
             annotation_col = coldata, 
             annotation_row = rowdata, 
             annotation_legend = FALSE,
             show_colnames = FALSE, 
             show_rownames = FALSE, 
             color = viridis(10),
             annotation_colors = annotation_colors,
             fontsize_col = 20,
             labels_row = FALSE, 
             labels_col = FALSE,
             legend = FALSE,
             annotation_names_row = FALSE)
  }else{
    p = ggplot()
  }
  
  return(p)
}
