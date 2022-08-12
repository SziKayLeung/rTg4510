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

## ----------Functions-----------------

# load all the functions
source_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential"
file.sources = list.files(path = source_dir, pattern="*.R", full = T)
sapply(file.sources,source,.GlobalEnv)

## ----------Plot colours-----------------

# plot label colour
label_colour <- function(genotype){
  if(genotype == "WT"){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "TG"){colour = wes_palette("Royal1")[2]}else{
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

# apply Genotype category
genotype_classification <- function(sample){
  classified_sample <- if(sample %in% c("O18","K18","S18","L22","Q20","K24")){"TG"
  } else if (x %in% c("Q21","K17","M21","O23","S23","K23")){ "WT"
  } else {"J20"
  }
  
  classified_sample <- unlist(classified_sample)
  return(classified_sample)
}


# number_of_reads
# Aim: huge wrapper function to read in multiple files from CCS, LIMA and REFINE directory, and output 5 plots
number_of_reads <- function(){
  
  ## Input CCS, LIMA, REFINE files
  # CCS, LIMA summary
  Reads <- input_isoseq_files("AllMouseTargeted_CCS_output.csv","AllMouseTargeted_LIMA_summary.csv")
  
  # classify Genotype in Reads
  Reads$Genotype <- lapply(Reads$sample, function(x)
    if(x %in% c("K18","K20","K24","L22","O18","O22","T20","Q20","S18","Q18","L18","T18")){"TG"
    } else if (x %in% c("K19","K23","K21","K17","S19","M21","O23","P19","Q21","S23","Q17","Q23")){ "WT"
    } else {"Batch"
    }
  )
  Reads$Genotype <- unlist(Reads$Genotype)
  
  # classify Genotype in CCS_values_mod for downstream plotting
  CCS_values_mod$Genotype <- lapply(CCS_values_mod$sample, function(x)
    if(x %in% c("K18","K20","K24","L22","O18","O22","T20","Q20","S18","Q18","L18","T18")){"TG"
    } else if (x %in% c("K19","K23","K21","K17","S19","M21","O23","P19","Q21","S23","Q17","Q23")){ "WT"
    } else {"Batch"
    }
  )
  CCS_values_mod$Genotype <- unlist(CCS_values_mod$Genotype)
  CCS_values_mod <<- CCS_values_mod
  
  return(Reads)
  
}


QC_yield_plot <- function(){
  #Reads <- number_of_reads()
  #write.csv(Reads,paste0(OUTPUT_DIR,"/Tg4510_IsoSeqTargetedReadsStats.csv"))
  Reads <- read.csv(paste0(OUTPUT_DIR,"/Tg4510_IsoSeqTargetedReadsStats.csv"))
  Reads$Description <- factor(Reads$Description, levels = c("Polymerase Reads","CCS Reads","FL Reads","FLNC Reads","Poly-A FLNC Reads","Transcripts"))
  
  Reads_plot <- Reads %>% mutate(Batch = word(Reads$variable, c(3), sep = fixed("_"))) %>% 
    full_join(., targetedpheno, by = c("Sample.ID" = "Sample")) %>%
    unite("Batch", Batch.x,Batch.y, na.rm = TRUE) %>%
    mutate(Batch = recode(Batch, "3b" = "3", "3a" = "3 (partial run)")) %>%
    mutate(Genotype = Phenotype)
  
  p1 <- Reads_plot %>% filter(Description != "Transcripts") %>% filter(Batch != "3 (partial run)") %>%
    ggplot(., aes(x = Description, y = value, colour = Batch, group = Batch)) +
    geom_line() + geom_point(size = 3) +  mytheme + theme(legend.position = c(0.8,0.8)) + labs(x = "", y = "Number of Reads (Thousands)") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_colour_discrete(name = "",labels = c("Batch 1 (n = 6)","Batch 2 (n = 9) ","Batch 3 (n = 9)"))
  
  p2 <- Reads_plot %>% filter(Description == "Transcripts") %>% 
    filter(!is.na(Genotype)) %>%
    ggplot(., aes(x = Genotype, y = value, colour = Genotype)) +
    geom_boxplot() + geom_point(size = 3) +  mytheme + 
    labs(x = "Genotype", y = "Number of FL Transcripts (Thousands)") +
    scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) + theme(legend.position = "none") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3,accuracy = 1))
  
  p3 <- Reads_plot[Reads_plot$Description == "Poly-A FLNC Reads" ,] %>% 
    filter(!is.na(Genotype)) %>%
    ggplot(., aes(x = Batch, y = value)) + geom_boxplot() + geom_point(aes(colour = Genotype),size = 3) + mytheme +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    labs(x = "Batch", y = "Number of Poly-A FLNC Reads (Thousands)") +
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) 
  
  cat("Number of CCS Reads in Batch 1, 2 and 3:", 
      sum(Reads[Reads$Description == "CCS Reads" & Reads$sample != "Targeted_Seq_3b_ccs_report" ,"value"]))
  cat("Number of CCS Reads in Batch 1, 2 and 3:\n")
  Reads[Reads$Description == "CCS Reads" & Reads$sample != "Targeted_Seq_3b_ccs_report" ,]
  
  cat("Sum number (Thousand) of Poly-A FLNC Reads across all Batches:", 
      sum(Reads_plot[Reads_plot$Description == "Poly-A FLNC Reads","value"])/1000,"\n")
  cat("Mean number (Thousand) of Poly-A FLNC Reads across all Batches:", 
      mean(Reads_plot[Reads_plot$Description == "Poly-A FLNC Reads","value"])/1000,"\n")
  cat("Min number (Thousand) of Poly-A FLNC Reads across all Batches:", 
      min(Reads_plot[Reads_plot$Description == "Poly-A FLNC Reads","value"])/1000,"\n")
  cat("Max number (Thousand) of Poly-A FLNC Reads across all Batches:", 
      max(Reads_plot[Reads_plot$Description == "Poly-A FLNC Reads","value"])/1000,"\n")
  for(i in 1:3){cat("Sum number (Thousand) of Poly-A FLNC Reads in Batch",i,":", 
                    sum(Reads_plot[Reads_plot$Description == "Poly-A FLNC Reads" & Reads_plot$Batch == i,"value"])/1000,"\n")}
  
  # not normally distributed therefore wilcoxon rank sum test
  #with(Reads_plot %>% filter(Description == "Transcripts"), shapiro.test(value[Phenotype == "WT"]))
  #with(Reads_plot %>% filter(Description == "Transcripts"), shapiro.test(value[Phenotype == "TG"]))
  #var.test(value ~ Phenotype,Reads_plot %>% filter(Description == "Transcripts")) #cannot assume variance
  wilcox.test(value ~ Phenotype,Reads_plot %>% filter(Description == "Transcripts")) 
  
  # correlation of FL transcripts and RIN
  transcript_RIN <- merge(Reads_plot %>% filter(Description == "Transcripts"),
                          tg4510_samples, by = "Sample.ID", all.x = T)
  #shapiro.test(transcript_RIN$value) # spearman's rank
  #shapiro.test(transcript_RIN$RIN) 
  cor.test(transcript_RIN$value,transcript_RIN$RIN, method = "spearman", exact = FALSE)
  
  return(list(p1,p2,p3))
}


on_target_plot <- function(){
  Probes <- merge(ldply(Probes_files , function(x) nrow(x)),
                  ldply(Probes_files , function(x) length(which(x$num_base_overlap != "0"))),by = ".id") %>%
    `colnames<-`(c("file", "Total_mapped_reads", "reads_probe_hit")) %>% 
    mutate(perc = reads_probe_hit/Total_mapped_reads * 100) %>%
    mutate(sample = word(.$file, c(1), sep = fixed("."))) %>% 
    full_join(., targetedpheno, by = c("sample" = "Sample"))
  
  p1<- ggplot(Probes, aes(x = as.factor(Batch), y = perc, fill = as.factor(Phenotype))) + geom_boxplot() +
    geom_point(aes(colour = as.factor(Phenotype)), position = position_jitterdodge(), size = 3) + 
    mytheme + labs(y = "On-Target Rate (%)", x = "Batch") + 
    scale_fill_manual(values = c(alpha(label_colour("TG"),0.4),alpha(label_colour("WT"),0.4)), name = "Genotype") + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), guide="none") 
  
  for(i in 1:3){cat("Mean on target rate in Batch",i,":", 
                    mean(Probes[Probes$Batch == i,"perc"]),"\n")}
  return(p1)
}


whole_vs_targeted_plots <- function(){
  cuff_tmap_exact = cuff_tmap[cuff_tmap$class_code == "=",]
  whole.class.files = whole.class.files %>% mutate(Matching = ifelse(isoform %in% cuff_tmap_exact$ref_id,"Both","Whole"))
  subsettargeted.class.files = subsettargeted.class.files %>% mutate(Matching = ifelse(isoform %in% cuff_tmap_exact$qry_id,"Both","Targeted")) 
  
  cols = c("isoform","associated_gene","Matching","structural_category")
  p1 = rbind(whole.class.files[whole.class.files$Matching != "Both",cols],subsettargeted.class.files[,cols]) %>% 
    filter(associated_gene %in% TargetGene) %>%
    group_by(associated_gene, Matching) %>% tally() %>%
    ggplot(., aes(x = reorder(associated_gene, -n), fill = Matching, y = n)) + geom_bar(stat = "identity") +
    mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "Number of Isoforms") +
    scale_fill_manual(name = "", values = c(label_colour("whole"),label_colour("whole+targeted"),label_colour("targeted")))
  
  p2 = subsettargeted.class.files %>% filter(associated_gene %in% TargetGene) %>%
    group_by(structural_category, Matching) %>% tally() %>% 
    ggplot(., aes(x = structural_category, y = n, fill = Matching)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "Structural Category", y = "Number of Isoforms \n Targeted Transcriptome") +
    scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("targeted"))) + 
    theme(legend.position = "top")
  
  print(rbind(whole.class.files[whole.class.files$Matching != "Both",cols],subsettargeted.class.files[,cols]) %>% 
          filter(associated_gene %in% TargetGene) %>%
          group_by(Matching) %>% tally())
  
  plot_exp <- function(dat, dataset){
    if(dataset == "Targeted"){colour <-  scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("targeted")))
    }else{colour <-  scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("whole")))}
    
    dat$ad_f = factor(dat$ADGene, levels=c("Target Genes","Not Target Genes"))
    p = ggplot(dat, aes(x = Matching, y = log10(FL), fill = Matching)) + geom_boxplot() + facet_grid(~ad_f) + mytheme + 
      labs(x = "", y = paste0("FL Read Counts \n",dataset," Transcriptome (Log10)")) + colour + theme(legend.position = "none")
    
    return(p)
  }
  
  p3 = plot_exp(whole.class.files,"Whole")
  p4 = plot_exp(subsettargeted.class.files,"Targeted")
  
  return(list(p1,p2,p3,p4))
}

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