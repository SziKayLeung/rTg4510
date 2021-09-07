# Szi Kay Leung
# Functions script for Thesis Chapter on Whole Transcriptome IsoSeq

# plot theme
loadfonts()
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=18,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 18,family="CM Roman"),
                 axis.text.x= element_text(size=16,family="CM Roman"),
                 axis.text.y= element_text(size=16,family="CM Roman"))

legend_theme <- theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.position="right",
                      legend.justification=c(0,1),
                      legend.margin=unit(1,"cm"),
                      legend.box="vertical",
                      legend.box.just = "left",
                      legend.key.size=unit(1,"lines"),
                      legend.text.align=0,
                      legend.background=element_blank())





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

# To scale axis into 1000s
ks <- function(x){ format(x/1000, big.mark=",")} 
perc_lab <- function(x){ format(x* 100, big.mark=",")} 

# apply Genotype category
genotype_classification <- function(sample){
  classified_sample <- if(sample %in% c("O18","K18","S18","L22","Q20","K24")){"TG"
  } else if (x %in% c("Q21","K17","M21","O23","S23","K23")){ "WT"
  } else {"J20"
  }
  
  classified_sample <- unlist(classified_sample)
  return(classified_sample)
}

density_plot <- function(dat,x.var,y.var, x_lab, y_lab,title){
  
  
  print(paste0(title))
  print(paste0("Correlation between", x.var, "and", y.var))
  
  print(cor.test(dat[[x.var]],dat[[y.var]]))
  cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  
  corr.value <- cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  p.value <- cor.test(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")$p.value 
  
  
  # corr.value <- cor(FSM_TPM$ISOSEQ_TPM_Normalised,FSM_TPM$RNASeq_TPM) # normalised ISOSEQ FL counts to length
  corr <- grobTree(textGrob(paste("r = ", round(corr.value, 2)), 
                            x = 0.05, y = 0.80, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 14, fontface = "italic",family="CM Roman")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  print(paste0("corr.value", corr.value))
  print(paste0("p.value", p.value))
  
  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1, name = "Density") +
    theme_bw() +
    labs(x = x_lab, y = y_lab, title = paste(title,"")) + 
    geom_smooth(method=lm, colour = "black") + 
    mytheme + 
    theme(legend.position = "none")
  
  return(p)
}



# number_of_reads
# Aim: huge wrapper function to read in multiple files from CCS, LIMA and REFINE directory, and output 5 plots
number_of_reads <- function(){
  
  ## Input CCS, LIMA, REFINE files
  # CCS, LIMA summary
  CCS <- read.csv(paste0(CCS_input_dir, "/WholeIsoSeqAll_CCS_output.csv"), header = T)
  LIMA <- read.csv(paste0(LiMA_input_dir, "/WholeIsoSeqAll_LIMA_summary.csv"), header = T)
  # REFINE summary input
  REFINE_json_list <- list.files(paste0(REFINE_input_dir), pattern = "flnc.filter_summary.json", full.names = T)
  REFINE_list <- lapply(REFINE_json_list , function(x) as.data.frame(fromJSON(file = x)))
  names(REFINE_list) <- list.files(paste0(REFINE_input_dir), pattern = "flnc.filter_summary.json")
  for(i in 1:length(names(REFINE_list))){
    REFINE_list[[i]]$Sample <- names(REFINE_list)[[i]]
  }
  # CLUSTER summary input
  CLUSTER_list_names <- list.files(paste0(CLUSTER_input_dir), pattern = ".cluster_report.csv", full.names = T)
  CLUSTER_list <- lapply(CLUSTER_list_names, function(x) read.csv(x))
  names(CLUSTER_list) <- list.files(paste0(CLUSTER_input_dir), pattern = ".cluster_report.csv")
  for(i in 1:length(names(CLUSTER_list))){
    CLUSTER_list[[i]]$Sample <- names(CLUSTER_list)[[i]]
  }
  
  # Cluster merged samples
  CLUSTER_merge <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq/MERGED_CLUSTER/WholeIsoSeq.clustered.cluster_report.csv")
  
  ## total number of transcripts (hq and lq)
  #CLUSTER_merge %>% group_by(cluster_id) %>% tally %>% nrow()
  ## number of transcripts per sample 
  transcript_per_sample <- CLUSTER_merge %>% mutate(sample = word(read_id, c(1), sep = "/")) %>% group_by(sample) %>% tally()
  transcript_per_sample <- merge(transcript_per_sample,sample_run, by.x = "sample", by.y = "run_id")
  res <- t.test(n ~ phenotype , data = transcript_per_sample, var.equal = TRUE)
  res
  
  ####################### CCS
  # Extract only values from mix of values and percentage
  CCS_values <- cbind(as.character(CCS[,1]), apply(CCS[,-1], 2, function(x) word(x, c(1), sep = fixed("("))))
  colnames(CCS_values)[1] <- "Description"
  CCS_values_mod <- as.data.frame(CCS_values) %>% melt(., id = "Description") %>%
    mutate(sample = word(variable, c(1), sep = fixed("_")))
  CCS_values <<- CCS_values
  
  
  ####################### LIMA
  # Extract only values from the mix of values and percentage from LIMA output
  LIMA_values <- data.frame(LIMA[,1],lapply(LIMA[,2:ncol(LIMA)], function(x) as.numeric(word(x, c(1),  sep = fixed ('(')) )))
  colnames(LIMA_values)[1] <- "Description"
  LIMA_values_mod <- as.data.frame(LIMA_values) %>% melt(., id = "Description") %>%
    mutate(sample = word(variable, c(1), sep = fixed("_")))
  
  
  # Refine into one dataframe and data-wrangled for easy merging
  REFINE <- do.call("rbind", REFINE_list) %>%
    rownames_to_column(., var = "variable") %>%
    mutate(sample = word(variable, c(1), sep = fixed("."))) %>%
    as.data.frame() %>%
    gather(., Description, value, 2:4, factor_key=TRUE) %>%
    .[,c("Description","variable","value","sample")] %>%
    mutate(value = as.numeric(value))
  
  
  CLUSTER <- do.call("rbind",CLUSTER_list) %>% 
    rownames_to_column(., var = "variable") %>%
    mutate(sample = word(variable, c(1), sep = fixed("."))) %>%
    as.data.frame() %>% 
    group_by(cluster_id,sample) %>%
    tally() %>% 
    group_by(sample) %>% tally() %>% mutate(Description = "Clustered_transcripts") %>% mutate(variable = "clustered.csv") %>% 
    .[,c(3,4,2,1)] 
  colnames(CLUSTER)[3] <- "value"
  
  Reads <-   
    CCS_values_mod[CCS_values_mod$Description == "ZMWs input               ",] %>%
    bind_rows(CCS_values_mod[CCS_values_mod$Description == "ZMWs pass filters        ",],)  %>%
    mutate(value = as.numeric(value)) %>%
    bind_rows(REFINE) %>%
    bind_rows(CLUSTER) %>%
    mutate(Description = as.character(Description))
  
  Reads$Description <- revalue(Reads$Description, c("ZMWs input               "="Polymerase Reads", "ZMWs pass filters        "="CCS Reads",
                                                    "num_reads_fl"="FL Reads", "num_reads_flnc"="FLNC reads",
                                                    "num_reads_flnc_polya" = "Poly-A FLNC reads",
                                                    "Clustered_transcripts" = "Transcripts"))
  levels(Reads$Description) <- c("Polymerase Reads","CCS Reads","FL Reads","FLNC reads","Poly-A FLNC reads","Transcripts")
  
  ### To calculate proportions to generate plots
  # total failed ccs reads
  failed_CCS_reads <- as.data.frame(CCS_values) %>% melt(., id = "Description") %>% filter(Description %in% c("ZMWs filtered       (C)  "))
  failed_LIMA_reads <- as.data.frame(LIMA_values) %>% melt(., id = "Description") %>% filter(Description %in% c("ZMWs below any threshold  (C) "))
  
  
  # classify Genotype in Reads
  Reads$Genotype <- lapply(Reads$sample, function(x)
    if(x %in% c("O18","K18","S18","L22","Q20","K24")){"TG"
    } else if (x %in% c("Q21","K17","M21","O23","S23","K23")){ "WT"
    } else {"J20"
    }
  )
  Reads$Genotype <- unlist(Reads$Genotype)
  
  # classify Genotype in CCS_values_mod for downstream plotting
  CCS_values_mod$Genotype <- lapply(CCS_values_mod$sample, function(x)
    if(x %in% c("O18","K18","S18","L22","Q20","K24")){"TG"
    } else if (x %in% c("Q21","K17","M21","O23","S23","K23")){ "WT"
    } else {"J20"
    }
  )
  CCS_values_mod$Genotype <- unlist(CCS_values_mod$Genotype)
  CCS_values_mod <<- CCS_values_mod
  
  return(Reads)
  
}

QC_yield_plot <- function(){
  # Reads <- number_of_reads()
  # write.csv(Reads,paste0(output_helpfig_dir,"/Tg4510_IsoSeqReadsStats.csv"))
  sequenced <- sequenced %>% filter(Genotype != "NA")
  Reads <- read.csv(paste0(output_helpfig_dir,"/Tg4510_IsoSeqReadsStats.csv"))
  
  p1 <- sequenced %>% ggplot(., aes(x = Genotype, y = Total.Bases..GB., fill = Genotype)) + geom_boxplot() +
    geom_jitter(shape=17, position=position_jitter(0)) +
    theme_bw() +
    labs(x = "Genotype", y = "Total Bases (Gb)") +
    mytheme + theme(legend.position = "none") +
    scale_fill_manual(values = c(label_colour("TG"),label_colour("WT")))
  
  p2 <- sequenced %>% ggplot(., aes(x = RIN, y = Total.Bases..GB., color = Genotype)) + geom_point(size = 2, shape = 17) +
    scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) +
    mytheme + labs(x = "RIN", y = "Total Bases (Gb)") 
  
  p3 <-
    Reads %>%
    filter(Genotype != "J20") %>% 
    filter(Description != "Transcripts") %>% 
    ggplot(., aes(x = reorder(Description, -value), y = value, group = sample, colour = Genotype)) +
    geom_line() + geom_point() +  mytheme + theme(legend.position = "top") + labs(x = "", y = "Number of Reads (Thousands)") +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p4 <- Reads %>%
    filter(Genotype != "J20") %>% 
    filter(Description == "Transcripts") %>% 
    ggplot(., aes(x = Genotype, y = value, colour = Genotype)) +
    geom_boxplot() + geom_point() +  mytheme + 
    labs(x = "", y = "Number of FL Transcripts (Thousands)", title = "\n\n") +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) + theme(legend.position = "none") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3,accuracy = 1), limits = c(30000,35000))
  
  ### Correlations 
  # Difference between WT and TG yield
  #with(sequenced, shapiro.test(Total.Bases..GB.[Genotype == "WT"])) # p = 0.4
  #with(sequenced, shapiro.test(Total.Bases..GB.[Genotype == "TG"])) # p = 0.7
  #res.ftest <- var.test(Total.Bases..GB.~ Genotype , data = sequenced)
  #res.ftest # p = 0.29
  res <- t.test(RIN ~ Genotype , data = sequenced, var.equal = TRUE)
  res
  
  res <- t.test(Total.Bases..GB.~ Genotype , data = sequenced, var.equal = TRUE)
  res
  
  return(list(p1,p2,p3,p4))
}

Mapping_stats_plots <- function(){
  Merged_mapping <- read.table(paste0(mapping_input_dir,"/WholeIsoSeq_reads_with_alignment_statistics.txt"))
  
  # duplicated multi-mapping
  Merged_Mapping_dup <- Merged_mapping %>% group_by(V1) %>% tally()
  nrow(Merged_Mapping_dup[Merged_Mapping_dup$n > 1,])
  
  # Threshold 
  Merged_mapping %>% filter(V8 > 0.85 & V6 > 0.95) %>% nrow()
  
  # plot of Alignable length and identity
  pMergedMap <- Merged_mapping %>%
    mutate(identity = V8 * 100) %>%
    mutate(length = V6 * 100) %>%
    ggplot(., aes(x = length, y = identity)) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1, name = "Density") +
    mytheme + labs(x = "Alignment Length (%)", y = "Alignment Identity (%)", title = "\n\n") +
    theme(legend.position = "none")
  
  return(pMergedMap)
}

lengths_plots <- function(){
  
  parse_lengths <- function(lengths_input_dir,suffix_name, type){
    
    # read in list of files ending with the input suffix name
    filenames <- as.list(list.files(path = lengths_input_dir, pattern = paste0(suffix_name,"$"), full.names = TRUE))
    files <- lapply(filenames, read.table)
    names(files) <- list.files(path = lengths_input_dir, pattern = paste0(suffix_name,"$"))
    
    # differentiate files prior to merging by creating column with the name of the file
    for(i in 1:length(files)){files[[i]]$File <- names(files)[[i]]}
    
    # bind all files
    all <- bind_rows(files)
    
    # create column for sample (note the separator varies between files for CCS and collapsed)
    if(type == "CCS"){
      all <- all %>% mutate(Sample = word(File, c(1), sep = fixed("."))) %>% mutate(Type = "CCS")
    } else if(type == "clustered"){
      all <- all %>% mutate(Sample = word(File, c(1), sep = fixed("."))) %>% mutate(Type = "Clustered")
    } else if(type == "collapsed"){
      all <- all %>% mutate(Sample = word(File, c(1), sep = fixed("."))) %>% mutate(Type = "Collapsed")
    } else {
      print("Type argument either CCS or clustered or collapsed")
    }
    
    # create column for genotype depending on the sample column
    all$Genotype <- ifelse(all$Sample %in% c("O18","K18","S18","L22","Q20","K24"), "TG","WT")
    
    return(all)
  }
  
  violin_plot <- function(input_dat, type){
    
    y_var <- if (type == "CCS"){ "CCS Read Length (kb)"
    } else if (type == "clustered"){ "Poly-A FLNC Read Length (kb)"
    } else if (type == "collapsed"){ "Collapsed Reads (Transcripts) (kb)"
    } else { print("Type argument either CCS or clustered or collapsed")
    }
    
    p <- input_dat %>%
      filter(!Sample %in% c("All_Merged","WT_Merged","TG_Merged")) %>%
      ggplot(., aes(x=reorder(Sample,V2, median), y=V2, fill = Genotype)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      labs(y = y_var, x = "") +
      theme_bw() + mytheme +
      theme(legend.title = element_blank(), legend.position = c(0.9, 0.9)) +
      scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
      scale_fill_manual(values = c(label_colour("WT"),label_colour("TG")))
    
    return(p)
    
  }
  
  #all_ccs <- parse_lengths(paste0(CCS_input_dir,"/Lengths"),"ccs.fasta.seqlengths.txt","CCS")
  #all_clustered <- parse_lengths(paste0(CLUSTER_input_dir,"/Lengths"),"clustered.hq.fasta.seqlengths.txt","clustered")
  #write.csv(all_ccs,paste0(output_helpfig_dir,"/Tg4510_CCSLengthsReadsStats.csv"))
  #write.csv(all_clustered,paste0(output_helpfig_dir,"/Tg4510_CLusteredLengthsReadsStats.csv"))
  
  all_ccs <- read.csv(paste0(output_helpfig_dir,"/Tg4510_CCSLengthsReadsStats.csv"))
  all_clustered <- read.csv(paste0(output_helpfig_dir,"/Tg4510_CLusteredLengthsReadsStats.csv"))
  p1 <- violin_plot(all_ccs, "CCS")
  p2 <- violin_plot(all_clustered, "clustered")
  
  # merged
  p3 <- bind_rows(all_ccs, all_clustered) %>%
    ggplot(., aes(x = Type, y = V2, fill = Genotype)) + geom_violin() +
    geom_violin() +
    labs(y = "Read Length (kb)", x = "") +
    theme_bw() + mytheme +
    theme(legend.title = element_blank(), legend.position = c(0.9, 0.9)) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    scale_fill_manual(values = c(label_colour("WT"),label_colour("TG"))) +
    scale_x_discrete(breaks=c("CCS","Clustered","Collapsed"),
                     labels=c("CCS Reads", "Poly-A FLNC Reads", "Collapsed Reads (Transcripts)"))
  
  return(list(p1,p2,p3))
}

iso_length <- function(class){
  class <- class %>% filter(subcategory != "mono-exon")
  p <- ggplot(class, aes(x = length)) + geom_histogram(bins = 15, fill="gray", col="black") + 
    labs(x = "Transcript Length (kb)", y = "Number of Isoforms (Thousand)") + mytheme +
    scale_x_continuous(labels = ks) + 
    scale_y_continuous(labels = ks) 
  
  two <- class[which(class$length >= 2000 & class$length <= 4000),] %>% nrow()
  print(paste0("Number of isoforms 2-4kb:", two, "(",round(two/nrow(class),2) *100,"%)"))
  
  return(p)
}

rarefaction_distribution <- function(){
  
  # Merge input files pertainng to genes and isoforms 
  all_rarefaction_levels <- bind_rows(all_rarefaction_genes, all_rarefaction_isoforms)
  
  ## Plots 
  # p1 <- isoform and gene level 
  # p2 <- isoform per category for individual datasets

  p1 <- all_rarefaction_levels %>% 
    ggplot(., aes(x = size, y = mean, linetype = type)) + 
    geom_line(size = 1.5) + 
    labs(x ="Number of Subsampled Reads (Thousand)", y = "Number of Genes/Isoforms (Thousand)") + 
    theme_bw() + mytheme + 
    scale_y_continuous(labels = ks) + scale_x_continuous(labels = ks) + 
    #scale_color_manual(values=c(label_colour("mouse")), name = "") + 
    scale_linetype_manual(values=c("dotted","solid")) +
    theme(legend.position = c(0.8, 0.6), legend.spacing.y = unit(-0.1, "cm"),legend.title = element_blank())
  
   
  p2 <- ggplot(all_rarefaction_isoforms_category, aes(x=size, y=mean, color=category)) + geom_line(aes(linetype = type), size = 1.5) + 
    labs(x = "Number of Subsampled Reads (Thousand)", 
         y = "Number of Isoforms (Thousand)", title = "") +
    mytheme + 
    scale_y_continuous(labels = ks, limits = c(0, 20000)) + 
    scale_x_continuous(labels = ks, expand = c(0.15, 0)) + 
    theme(legend.position = "none") +
    scale_linetype_manual(values=c("solid", "longdash","dotted")) +
    geom_dl(aes(label = category),  method = list(dl.trans(x = x + 0.5), "last.bumpup", cex = 1.3, hjust = .1))
    
  
  return(list(p1,p2))
}


no_of_isoforms_sample <- function(class){
  # number of samples with detected expression of isoform
  # prerequisite: demultiplexing with numver of counts per sample
  # across each row (i.e isoform, count the number of occurences where reads are != 0)
  dat <- class %>% dplyr::select(starts_with("FL.")) %>% 
    mutate(median_FL = apply(.,1, function(x) median(x)), num_samples = apply(.,1, function(x) length(x[which(x != "0")]))) 
  table(dat$num_samples)
  table(dat$num_samples)/sum(table(dat$num_samples))
  p1 <- ggplot(dat, aes(x = as.factor(num_samples))) + geom_bar(aes(y = (..count..)/sum(..count..))) + 
    scale_y_continuous(labels = perc_lab)  + mytheme + labs(x = "Number of Samples", y = "Isoforms (%)")
  
  p2 <- ggplot(dat, aes(x = as.factor(num_samples), y = log(median_FL))) + geom_boxplot() + 
    mytheme + labs(x = "Number of Samples", y = "Median FL Read Count(Log10)")
  
  return(list(p1,p2))
}

novel_annotated <- function(){
  subset_feature <- function(col_name_feature, ylabel, category){
  
  # subset classification files (filtered) by annotated genes, and novel transcripts 
  annotated.class.files <- class.files[!grepl("NOVELGENE",class.files$associated_gene),]
  annotated.class.files.novel.transcripts <- annotated.class.files[grepl("novel",annotated.class.files$associated_transcript),] %>% 
                                                      filter(subcategory != "mono-exon")
  annotated.class.files.annotated.transcripts <- annotated.class.files[!grepl("novel",annotated.class.files$associated_transcript),] %>% 
                                                          filter(subcategory != "mono-exon")
  
  subset_transcripts <- function(dataset, type){
    if(dataset == "Novel"){
      dat <- annotated.class.files.novel.transcripts %>% filter(structural_category == type)
    } else {
      dat <- annotated.class.files.annotated.transcripts %>% filter(structural_category == type)
    }
    return(dat)
  }
  
  #annotated.class.files.NIC.transcripts <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% filter(structural_category == "NIC"))
  #annotated.class.files.NNC.transcripts <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% filter(structural_category == "NNC"))
  # annotated.class.files.Fusion.transcripts <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% filter(structural_category == "Fusion"))
  
  
  # Subset length from classification file (bind_rows as input is list of dataframes)
  extract_feature <- function(dat, dat_name){
    dat1 <- dat %>% bind_rows() %>% 
      mutate(Transcripts = dat_name) %>% .[,c(col_name_feature, "Transcripts")]
  }
  
  y.var <- rlang::sym(quo_name(enquo(col_name_feature)))
  
  if(category == "all_novel") { 
    set1 <- extract_feature(annotated.class.files.annotated.transcripts, "Known")
    set2 <- extract_feature(annotated.class.files.novel.transcripts, "Novel") 
    merge <- bind_rows(list(set1, set2))
    p <- merge %>%  
      ggplot(., aes(x = Transcripts, y = !! y.var, fill = Transcripts)) + 
      geom_boxplot() + theme_bw() + 
      mytheme + labs(x = "", y= ylabel) + 
      scale_fill_manual(values = c(label_colour("known"),label_colour("novel"))) +
      theme(legend.position = "none") + 
      guides(fill=guide_legend(nrow=3,byrow=TRUE))
    
  } else { 
    set1 <- extract_feature(subset_transcripts("Known", "FSM"), "FSM")
    set2 <- extract_feature(subset_transcripts("Known", "ISM"), "ISM")
    set3 <- extract_feature(subset_transcripts("Novel", "NIC"), "NIC")
    set4 <- extract_feature(subset_transcripts("Novel", "NNC"), "NNC")
    set5 <- extract_feature(subset_transcripts("Novel", "Fusion"), "Fusion")
    merge <- bind_rows(list(set1, set2, set3, set4, set5)) 
    merge$Transcripts <- factor(merge$Transcripts, levels = c("FSM", "ISM", "NIC", "NNC","Fusion"))
    p <- merge %>%  
      ggplot(., aes(x = Transcripts, y = !! y.var, fill = Transcripts)) + 
      geom_boxplot() + theme_bw() + 
      scale_fill_manual(values = c(alpha(label_colour("known"),0.8),alpha(label_colour("known"),0.3),
                                   alpha(label_colour("novel"),0.8),alpha(label_colour("novel"),0.5),alpha(label_colour("novel"),0.3))) +
      mytheme + labs(x = "", y= ylabel) + 
      theme(legend.position = "none") + 
      guides(fill=guide_legend(nrow=3,byrow=TRUE))
  }
  return(p)
  }
  
  expression <- subset_feature("Log_ISOSEQ_TPM", "Iso-Seq Expression (Log10 TPM)", "all_novel")
  split_expression <- subset_feature("Log_ISOSEQ_TPM", "Iso-Seq Expression (Log10 TPM)", "all_split")
  
  length <- subset_feature("length", "Transcript Length (kB)", "all_novel") + scale_y_continuous(labels = ks)
  split_length <- subset_feature("length", "Transcript Length (kB)", "all_split") + scale_y_continuous(labels = ks)
  
  exon <- subset_feature("exons", "Number of Exons", "all_novel")
  split_exon <- subset_feature("exons", "Number of Exons", "all_split")
  
  output <- list(expression,split_expression,length,split_length,exon,split_exon)
  return(output)
}

# summary_info <read.clasification.file> 
# Aim: summarise the number of isoforms, minimum number of exons and max etc per gene 
summary_info <- function(dat){
  
  total_fl <- sum(dat$FL, na.rm=T)
  
  info <- list(
    # Number of isoforms 
    dat %>% count(associated_gene, name = "Num_of_Isoforms"),
    # Min Exons 
    dat %>% group_by(associated_gene) %>% summarise(Min_exon = min(exons)),
    # Max Exons 
    dat %>% group_by(associated_gene) %>% summarise(Max_exon = max(exons)),
    # Min Length 
    dat %>% group_by(associated_gene) %>% summarise(Min_length = min(length)),
    # Max Length
    dat %>% group_by(associated_gene) %>% summarise(Max_length = max(length)), 
    # Num of FL reads 
    dat %>% group_by(associated_gene) %>% summarise(Total_Reads = sum(FL))
  )
  
  final <- Reduce(function(...) merge(..., by='associated_gene', all.x=TRUE), info) %>%
    filter() 
  
  # TPM by gene (sum of FL reads)
  final <- final %>% mutate("FL_TPM" = round(Total_Reads*(10**6)/total_fl)) %>%
    mutate("Log10_FL_TPM" = log10(FL_TPM))
  
  return(final)
}


# exon_length_isoform_correlation <read_classification_file>
# Aim: Correlate the number of exons/length with number of isoforms (only multiexonic only)
# Mulitple plots and .txt file for the output of correlation of stats
exon_length_isoform_correlation <- function(){
  
  # remove mono-exonic transcripts 
  dat1 <- class.files %>% filter(subcategory != "mono-exon")
  
  # summary_info = function to summarise the number of isoforms, etc
  # note TPM calculated with the removal of monoexonic transcripts 
  dat2 <- data.frame(summary_info(dat1))
  
  # Gene expression cut off threshold at Log10_FL_TPM at >2.5 
  dat3 <- dat2 %>% filter(Log10_FL_TPM > 2.5)
  
  
  # Plots 
  # P1: transcript length vs number of exons for all transcripts 
  # P2: transcrpt length vs number of exons (using only representative transcript per gene)
  # P3: same plot as P2 but filtered at expression level 
  # P4: transcript length vs number of isoforms (using only representative length) --> SUPP
  # P5: same plot as P4 but filtered at expression level --> SUPP
  # P6: exon vs number of isoforms (exon defined by representative transcript) --> SUPP
  # P7: same plot as P7 but filtered at expression level --> SUPP
  p1 <- density_plot(dat1, "length","exons","Transcript length (kb)","Number of exons", 
                     "Transcript length\nvs Number of Exons (all transcripts)")  + labs(title="\n\n\n")
  
  p2 <- density_plot(dat2, "Max_length","Max_exon", "Transcript length (kb)", "Number of exons",
                    "Transcript length \nvs Number of Exons (representative transcript per gene)") + labs(title="") +
    scale_x_continuous(labels = ks)
  
  p3 <- density_plot(dat3, "Max_length","Max_exon", "Transcript length (bases)", "Number of exons",
                    "Transcript Length\nvs Number of Exons (representative transcript per gene):\nFiltered at 25TPM")  + labs(title="\n\n\n")
  
  # Transcript length (max representative)\nvs Number of Isoforms
  p4 <- density_plot(dat2, "Max_length","Num_of_Isoforms", 
                     "Transcript Length (kb)", "Number of Isoforms", "Transcript length (max representative) vs Number of Isoforms") + 
    scale_x_continuous(labels = ks)  + labs(title="")
  
  # Transcript length (max representative)\nvs Number of Isoforms:\nFiltered for high gene expression
  p5 <- density_plot(dat3, "Max_length","Num_of_Isoforms", 
                     "Gene Length (kb)", "Number of Isoforms",
                     " Transcript length (max representative) vs Number of Isoforms:Filtered for high gene expression") + 
    scale_x_continuous(labels = ks) + labs(title="\n\n\n")
  
  # Number of Exons (max representative)\nvs Number of Isoforms
  p6 <- density_plot(dat2, "Max_exon","Num_of_Isoforms", "Number of Exons", "Number of Isoforms","Number of Exons (max representative) vs Number of Isoforms") + labs(title="")
  
  # Number of Exons (max representative)\nvs Number of Isoforms:\nFiltered for high gene expression
  p7 <- density_plot(dat3, "Max_exon","Num_of_Isoforms", "Number of Exons", "Number of Isoforms","Number of Exons (max representative) vs Number of Isoforms:Filtered for high gene expression") + labs(title="\n\n\n")

  
  # linear regression 
  # summary(lm(Num_of_Isoforms ~ Max_length + Max_exon, dat2))
  return(list(p1,p2,p3,p4,p5,p6,p7))
}

### ERCC #################################################################
run_ERCC <- function(class){
  
  ERCC_conc <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/ERCC/ERCC_calculations.csv", header = T)[-1,]
  
  cat("Total unique ERCCs:", length(unique(class$chrom)), paste0("(", round(length(unique(class$chrom))/92 * 100,2), "%)"))
  
  # redundant 
  redundant <- class[,c("chrom")] %>% table() %>% reshape2::melt(.)
  colnames(redundant) <- c("ERCC","num_isoforms")
  
  p1 <- redundant %>% 
    group_by(num_isoforms) %>% tally() %>% ggplot(., aes(x = as.factor(num_isoforms), y = n)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "ERCC", y = "Number of Isoforms")
  
  # isoform vs concentration
  isoform_conc <- merge(redundant, ERCC_conc, by.x = "ERCC", by.y = "ERCC_ID", all = TRUE) %>%
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>%
    replace_na(list(num_isoforms = 0)) %>% 
    mutate(num_isoforms = as.factor(num_isoforms))
  
  p2 <- ggplot(isoform_conc, aes(x = num_isoforms, y = log2_amount_of_ERCC, colour = num_isoforms)) + geom_jitter(width = 0.2) + theme(legend.position = "none") + 
    mytheme + labs(x = "Number of Isoforms", y = "Amount of ERCC (Log2)") + theme(legend.position = "none")
  
  # correlation 
  ERCC_corr <- merge(class,ERCC_conc, by.x = "chrom", by.y = "ERCC_ID") %>% 
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>%
    mutate(log2_FL_reads = log2(FL)) 
  
  p3 <- density_plot(ERCC_corr,"log2_amount_of_ERCC","log2_FL_reads", "Amount of ERCC (Log2)", "Number of FL Reads (Log2)","") 
  p3
  
  return(list(p1,p2,p3))
}

rnaseq_isoseq_transcriptome <- function(cuffrefmap_input,cufftmap_input){
  #### READ Files from Gff compare 
  # note only include isoforms from Iso-Seq that are partially or fully matched to RNA-Seq, not all isoforms from Iso-Seq dataset 
  cuffrefmap <- read.table(cuffrefmap_input, header = T)
  # no duplicates of isoform - QC checked
  # cuffrefmap[cuffrefmap$class_code == "=",] %>% group_by(ref_id) %>% tally() %>% filter(n > 1)
  
  cufftmap <- read.table(cufftmap_input, header = T)
  # Note not all isoforms are detected by RNA-Seq therefore not listed in cufftmap
  # novel_gene <- class.files[grepl("NOVE", class.files$associated_gene),"isoform"]
  # for(i in novel_gene){print(cufftmap[cufftmap $ref_id == i,])}
  
  # classification of rnaseq reads to isoseq 
  match = wes_palette("Darjeeling1")[2]
  pot_match = wes_palette("Darjeeling1")[4]
  unknown = wes_palette("Darjeeling1")[1]
  frag = wes_palette("Royal1")[1]
  class_cuff <- cufftmap %>% group_by(class_code) %>% tally() %>% mutate(perc = n/sum(n) * 100) 
  for(n in 1:nrow(class_cuff)){
    class_cuff$Classification[n] = if(class_cuff$class_code[n] == "="){"match"}else if(class_cuff$class_code[n] == "u"){"unknown"}else if(
      class_cuff$class_code[n] == "j"){"pot"} else {"frag"}}
  class_cuff$Classification <- factor(class_cuff$Classification,levels = c("pot","unknown","match","frag"),)
  
  p0 <- ggplot(class_cuff, aes(x = reorder(class_code, -perc), y = perc, fill = Classification)) + geom_bar(stat = "identity") +
    scale_fill_manual(values = c(pot_match, unknown, match,frag), 
                      labels = c("Close Match to Iso-Seq Isoforms","No alignment to Iso-Seq Isoforms",
                                 "Complete Match","Likely fragments of Iso-Seq Isoforms")) + mytheme + 
    labs(x = "GffCuffcompare's Transcript Classification Codes", y = "RNA-Seq Isoforms (%)") + 
    theme(legend.position = c(0.8,0.8))
  # cufftmap %>% filter(class_code %in% c("u")) %>% ggplot(.,aes(x = num_exons)) + geom_bar(aes(y = (..count..)/sum(..count..)))
  
  ## 
  # Replace id of those rnaseq isoforms that fully match with isoseq isoforms with pbid for venn diagram
  rnaseq <- c(as.character(cufftmap[cufftmap$class_code != "=","qry_id"]),as.character(cuffrefmap[cuffrefmap$class_code == "=","ref_id"])) 
  p1 <- venn.diagram(
    x = list(rnaseq, class.files$isoform), category.names = c("RNA-Seq","Iso-Seq"), filename = NULL, output=TRUE,
    lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC"), main = "\n",
    cex = 1,fontface = "bold",fontfamily = "ArialMT",
    cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-27, 27),  cat.dist = c(0.055, 0.055),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n",
    print.mode = "raw"
  )
  
  
  ## Number of Isoforms per dataset 
  num_iso <- list(class.files %>% mutate(Sample = "Iso-Seq") %>%  .[,c("isoform", "associated_gene", "novelGene","FSM_class","gene_exp","Sample")],
                  rnaseq.class.files %>%  .[,c("isoform","associated_gene", "novelGene","FSM_class","gene_exp","Sample")])
  isoPerGene <- lapply(num_iso, function(x) SQANTI_gene_preparation(x)) %>% bind_rows()
  
  # Total Number of Genes per Type 
  Total_Num <- isoPerGene %>% group_by(Sample) %>% count(Sample)
  # Number of isoform cateogories per Sample
  p2 <- isoPerGene %>% group_by(Sample) %>% count(nIsoCat) %>% full_join(Total_Num,., by = "Sample") %>% mutate(Perc = n.y/n.x * 100) %>% 
    ggplot(., aes(x=nIsoCat, fill=Sample)) +
    geom_bar(stat="identity", aes(y= Perc, group = as.factor(Sample)), color="black", size=0.3, width=0.7, 
             position="dodge") + 
    labs(x ="Number of Isoforms", y = "Genes (%)", fill = "", title = "\n") +
    #scale_fill_manual(values=c(label_colour(dataset1), label_colour(dataset2)), 
    #                 labels=c(paste(dataset1,"Cortex"), paste(dataset2,"Cortex")))  + 
    mytheme + 
    theme(legend.position = c(0.75,0.95))
  
  # Length
  Lengthcomp <- bind_rows(class.files %>% mutate(Sample = "Iso-Seq") %>% .[,c("length","Sample")], rnaseq.class.files %>% .[,c("length","Sample")])
  p3 <-   ggplot(Lengthcomp, aes(x = Sample, y = log10(length/1000))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Isoform Length (Log10 kb)", title = "\n")
  res = t.test(length ~ Sample, data = Lengthcomp, var.equal = TRUE)
  res 
  #res$p.value
  
  # Exons
  Exoncomp <- bind_rows(class.files %>% mutate(Sample = "Iso-Seq") %>% .[,c("exons","Sample")], rnaseq.class.files %>% .[,c("exons","Sample")]) 
  p4 <- ggplot(Exoncomp, aes(x = Sample, y = log10(exons))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Number of Exons (Log10)", title = "\n")
  res = t.test(exons ~ Sample, data = Exoncomp, var.equal = TRUE)
  res 
  #res$p.value
  
  
  # cage peak
  cage <- data.frame(dataset = c("Iso-Seq","RNA-Seq"),
                     within_50 = c(nrow(class.files %>% filter(abs(dist_to_cage_peak) <= 50)),
                                   nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50))),
                     # include NAs
                     without_50 = c(nrow(class.files) - nrow(class.files %>% filter(abs(dist_to_cage_peak) <= 50)),
                                    nrow(rnaseq.class.files) - nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50)))) %>% 
    remove_rownames %>% column_to_rownames(var="dataset")
  cage
  res <- fisher.test(cage)
  res 
  #res$p.value
  
  
  p5 <- data.frame(dataset = c("Iso-Seq","RNA-Seq"),
                   within_50 = c(nrow(class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(class.files) * 100,
                                 nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(rnaseq.class.files) * 100),
                   # include NAs
                   without_50 = c(100 - nrow(class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(class.files) * 100,
                                  100 - nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(rnaseq.class.files) * 100)) %>%
    melt() %>%
    mutate(variable = factor(variable, levels = c("without_50","within_50"))) %>%
    ggplot(., aes(x = dataset, y = value, fill = variable)) + geom_bar(stat = "identity") + mytheme +
    labs(y ="Isoforms within 50bp CAGE (%)", x = "", fill = "", title = "\n") + theme(legend.position = "right") +
    scale_fill_discrete(labels = c("No","Yes"))
  
  # structural_category
  p6 <- bind_rows(class.files %>% mutate(Sample = "Iso-Seq") %>% .[,c("structural_category","Sample")] %>% 
                    group_by(structural_category,Sample) %>% tally() %>% group_by(Sample) %>% mutate(perc = n/sum(n)*100),
                  rnaseq.class.files %>% .[,c("structural_category","Sample")] %>% group_by(structural_category,Sample) %>% tally() %>% group_by(Sample) %>% mutate(perc = n/sum(n)*100)) %>% 
    ggplot(., aes(x = Sample, y = perc, fill = structural_category)) + geom_bar(stat = "identity") + mytheme +
    labs(x = "", y = "Isoforms (%)", fill = "Structural Category", title = "\n") + theme(legend.position = "bottom")
  
  output <- plot_grid(grobTree(p1),p2,p3,p4,p5,p6, labels = "auto", label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
  
  return(output)
}


### Alternative Splicing #################################################################
AS_genes_events <- function(){
  cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  splicing_events <- read.csv(paste0(input_table_dir,"/AS_IR/ALL_SUPPA2_Genes_Output_Updated2.csv"))
  dataset_tally_events <- splicing_events %>% group_by(Sample) %>% tally(n)  
  
  anno_novel_AS <- read.csv(paste0(input_table_dir,"/AS_IR/Mouse_AS_NovelAnnoIsoforms.csv")) %>% gather(., Type, n, known:novel, factor_key=TRUE)
  anno_novel_AS_tally <- anno_novel_AS %>% group_by(Type) %>% tally(n)
  
  # Number of splicing events - all (inc Novel Genes)
  splicing_events_tally <- splicing_events %>% group_by(Event, Sample) %>% tally(n) %>% left_join(dataset_tally_events, by = "Sample") %>% mutate(perc = n.x/n.y * 100) %>% filter(Sample == "Mouse") %>% mutate(Type = "All") %>% mutate(Group = "1") %>% select(!Sample)
  
  print(splicing_events_tally)
  
  # Number of splicing events - annotated vs novel isoforms of annotated genes
  anno_novel_AS_tally_present <- anno_novel_AS  %>% left_join(anno_novel_AS_tally, by = "Type") %>% mutate(perc = n.x/n.y * 100) %>% mutate(Group = "2")
  
  # Number of splicing events - annotated genes
  knowngenes_AS <- read.csv(paste0(input_table_dir,"/AS_IR/Mouse_AS_KnownGenes.csv")) %>% mutate(Group = "1", Type = "Annotated Genes")
  
  p1 <- bind_rows(knowngenes_AS[,c("Event","Type","perc","Group")],anno_novel_AS_tally_present[,c("Event","Type","perc","Group")]) %>%
    ggplot(., aes(x = Type, y = perc, fill = reorder(Event, -perc))) + geom_bar(stat = "identity") + 
    theme_bw() + labs(y = "Splicing Events (%)", x = "") + mytheme + 
    theme(legend.position = "bottom",legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1)) +
    facet_grid(~ Group,scales='free') +
    scale_x_discrete(labels = c("known" = "Known Isoforms", "novel" =  "Novel Isoforms")) +
    theme(strip.background = element_blank(),strip.text.x = element_blank())
  
  dataset_tally_gene <- splicing_events %>% group_by(Sample) %>% count(associated_gene) %>% tally()
  splicing_events_number <- splicing_events %>% group_by(Sample, associated_gene) %>% tally() %>% group_by(n, Sample) %>% 
  tally() %>% left_join(dataset_tally_gene, by = "Sample") %>%
  mutate(perc = nn/n.y * 100)  %>% `colnames<-`(c("Number_of_Splicing_Events", "Sample", "Genes", "Total_Genes","perc")) %>% as.data.frame(.) 
  splicing_events_number <<- splicing_events_number
  p2 <- splicing_events_number %>%
    filter(Sample %in% c("Mouse")) %>%
    ggplot(., aes(x = Number_of_Splicing_Events, y = perc)) + 
    geom_bar(stat = "identity", position = position_dodge()) + 
    theme_bw() + labs(y = "AS Genes (%)", x = "Number of Splicing Events") + mytheme + 
    theme(legend.position = c(0.85,0.85), legend.title = element_blank()) + 
    scale_x_continuous(breaks = 1:7) 
  
  return(list(p1,p2))
}
### IR and NMD#################################################################
IR_NMD_run <- function(df){
  
    library(RColorBrewer)
    myCol <- brewer.pal(3, "Set2")
    
    annotated_genes <- df[!grepl("NOVEL", df$associated_gene), ]
    IR_NMD <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(subcategory == "intron_retention" & predicted_NMD == "TRUE")
    IR_Not_NMD <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(subcategory == "intron_retention" & predicted_NMD == "FALSE")
    NMD <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(predicted_NMD == "TRUE")
    NMD_Not_IR <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(predicted_NMD == "TRUE") %>% filter(subcategory != "intron_retention")
    IR_coding <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(subcategory == "intron_retention" & coding == "coding")
    
    #plot[[i]] <- venn_diagram_plot_twocircles(unique(IR_coding$associated_gene), unique(NMD$associated_gene), "Genes with Intron Retention", "Genes with NMD")
    p <- venn.diagram(
      x = list(unique(IR_coding$associated_gene), unique(NMD$associated_gene), unique(IR_NMD$associated_gene)),
      category.names = c("IR", "NMD", "IR-NMD"),
      filename = NULL,output=TRUE, print.mode = "raw",fill = myCol, cex = 2,fontface = "bold",fontfamily= "CM Roman",
      # Set names
      cat.cex = 2,cat.fontfamily = "CM Roman")
    
  return(p)
}

NMD_vs_NonNMD <- function(){

    # subset transcripts whether from annotated or novel genes
    annotated_genes <- class.files[!grepl("NOVEL", class.files$associated_gene),]
    novel_genes <- class.files[grepl("NOVEL", class.files$associated_gene),]
    
    # Annotated Genes, NMD transcripts
    NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE")
    non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE")
    
    IR_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE" & subcategory == "intron_retention")
    Non_IR_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE" & subcategory != "intron_retention")
    IR_Non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE" & subcategory == "intron_retention")
    Non_IR_Non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE" & subcategory != "intron_retention")
    
    
    # Novel Genes, NMD transcripts 
    NMD_transcripts_novel_genes <- novel_genes %>% filter(predicted_NMD == "TRUE")

    
    p1 <- rbind(data.frame("Num" = IR_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM, "Type" = "IR-NMD"),
                data.frame("Num" = IR_Non_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM, "Type" = "IR"),
                data.frame("Num" = Non_IR_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM,"Type" = "NMD"),
                data.frame("Num" = Non_IR_Non_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM, "Type" = "Non IR-NMD")) %>%
      ggplot(., aes(y = Num, x = Type)) + geom_boxplot() +
      theme_bw() +
      labs(x = "", y = "Iso-Seq Expression (Log10 TPM)") + mytheme +
      theme(legend.position = c(.90, 0.90), legend.title = element_blank()) 

  return(p1)
}


### lncRNA #################################################################
lncRNA <- function(){
  
  # Internal Function
  # statstical test between lncRNA and non-lncRNA
  lncRNA_tests <- function(all_dat, lnc_dat, title){

    if(title == "Transcript Expression (Annotated Genes only)"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(y = Num, x = Sample)) + geom_boxplot() +
        theme_bw() +
        labs(x = "", y = "Iso-Seq Expression (Log10 TPM)") +
        mytheme + scale_fill_manual(values = wes_palette("Darjeeling2")) +
        theme(legend.position = c(.15, 0.85), legend.title = element_blank()) + theme(legend.direction="horizontal")
    }else if(title == "Number of Exons (Annotated Genes only)"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.2) +
        theme_bw() +
        labs(x = "Number of Exons", y = "Density") +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.15, 0.85), legend.title = element_blank()) + theme(legend.direction="horizontal")
    }else if(title == "ORF Length (Annotated Genes only)"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.2) +
        theme_bw() +
        labs(x = "ORF Length (kb)", y = "Density") +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.15, 0.85), legend.title = element_blank()) + theme(legend.direction="horizontal")
    }else if(title == "Number of Isoforms per Gene"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(y = Num, x = Sample)) + geom_boxplot() +
        theme_bw() +
        labs(x = "", y = "Number of Isoforms per Gene (Log10)") +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.25, 0.85), legend.title = element_blank()) + theme(legend.direction="horizontal")
    }
    else{
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.2) +
        theme_bw() +
        labs(x = "Transcript Length (kb)", y = "Density") +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.45, 0.85), legend.title = element_blank()) + theme(legend.direction="horizontal")
    }
  }
  
 
    # input file = annotated genes only 
    input_file <- class.files[!grepl("NOVEL", class.files$associated_gene),]
    # lncfiles = classifiction files aligned to lncRNA reference annotation
    input_lnc_file <- lnc.files
    # transcrips annotated genes = known transcripts of annotated genes 
    transcripts_annotated_genes <- input_file[!grepl("NOVEL", input_file$associated_gene),]
    # genes_annotated_genes = tally of the known transcripts of annotated genes
    genes_annotated_genes <- transcripts_annotated_genes %>% group_by(associated_gene) %>% tally()
    # coding_transcripts_annotated_genes = coding, known transcripts of annotated genes 
    coding_transcripts_annotated_genes <- input_file[!grepl("NOVEL", input_file$associated_gene) & input_file$coding == "coding",]
    coding_genes_annotated_genes <- coding_transcripts_annotated_genes %>% group_by(associated_gene) %>% tally()
    
    # annotated lncRNA transcripts and genes from lnc.files 
    lncRNA_transcripts <- input_lnc_file[!grepl("novel",input_lnc_file$associated_gene),]
    lncrna_genes <- unique(input_lnc_file[!grepl("novel",input_lnc_file$associated_gene),"associated_gene"])
    
    # subset of lncRNA transcripts from original classification file using genes extracted from lnc.files 
    class_non_lncRNA_transcripts <- input_file[!input_file$associated_gene %in% lncrna_genes,]
    class_lncRNA_transcripts <- input_file[input_file$associated_gene %in% lncrna_genes,]
    
    # number of isoforms (n) of lncRNA genes vs nonlncRNA genes for each gene
    isoform_diversity <- bind_rows(class_non_lncRNA_transcripts %>% group_by(associated_gene) %>% tally() %>% mutate(gene_type = "non_lncRNA"),
                                   class_lncRNA_transcripts %>% group_by(associated_gene) %>% tally() %>% mutate(gene_type = "lncRNA")) 
    isoform_diversity_lncRNA <- isoform_diversity[isoform_diversity$gene_type == "lncRNA",]
    isoform_diversity_nonlncRNA <- isoform_diversity[isoform_diversity$gene_type == "non_lncRNA",]
    
    
    # Ensure non-lncRNA dataset first
    # p1: Transcripts Length (Annotated Genes only)
    # p2: Number of Exons (Annotated Genes only)
    # p3: Non-Coding Transcripts Length (Annotated Genes only)
    # p4: Coding Transcripts Length (Annotated Genes only)
    # p5: Transcript Expression (Annotated Genes only)
    # p6: ORF Length (Annotated Genes only)
    # p7: Number of Isoforms per Gene
    
    
    # comparing transcript length of lncRNA vs non-lncRNA (irrespective of coding or not)
    p1 <- lncRNA_tests(class_non_lncRNA_transcripts$length, class_lncRNA_transcripts$length, 
                       "Transcripts Length (Annotated Genes only)") + 
      scale_x_continuous(labels = ks) + theme(legend.position="none")
    
    p2 <- lncRNA_tests(class_non_lncRNA_transcripts$exons, class_lncRNA_transcripts$exons, 
                       "Number of Exons (Annotated Genes only)") + theme(legend.position="none")
    
    p3 <- lncRNA_tests(class_non_lncRNA_transcripts[class_non_lncRNA_transcripts$coding == "non_coding","length"],
                       class_lncRNA_transcripts[class_lncRNA_transcripts$coding == "non_coding","length"], 
                       "Non-Coding Transcripts Length (Annotated Genes only)") +
      scale_x_continuous(labels = ks) 
    
    p4 <- lncRNA_tests(class_non_lncRNA_transcripts[class_non_lncRNA_transcripts$coding == "coding","length"],
                       class_lncRNA_transcripts[class_lncRNA_transcripts$coding == "coding","length"], 
                       "Coding Transcripts Length (Annotated Genes only)") +
      scale_x_continuous(labels = ks) + theme(legend.position="none")
    
    
    p5 <- lncRNA_tests(class_non_lncRNA_transcripts$Log_ISOSEQ_TPM, class_lncRNA_transcripts$Log_ISOSEQ_TPM, 
                       "Transcript Expression (Annotated Genes only)") + theme(legend.position="none")
    
    p6 <- lncRNA_tests(class_non_lncRNA_transcripts$ORF_length, class_lncRNA_transcripts$ORF_length, 
                       "ORF Length (Annotated Genes only)") +   scale_x_continuous(labels = ks) + theme(legend.position="none")
    
    # Log number for plotting, makes no difference to statstical test
    p7 <-lncRNA_tests(log10(isoform_diversity_nonlncRNA$n), log10(isoform_diversity_lncRNA$n), 
                      "Number of Isoforms per Gene") + theme(legend.position="none")
    
    
    plots <- list(p1,p2,p3,p4,p5,p6,p7)
  
  return(plots)
}

### Human MAPT #################################################################
## Determine the number of Cluster reads with human-specific MAPT sequence for each sample
#hMAPT.header for each file contains multiple HQ, FL-polished transcripts (different transcript names e.g "@transcriptX", "@transcriptY"), but which all contain the same hMAPT sequence. Reason that there are multiple transcripts is due to collapsed properly (redundancy) from Iso-Seq3. For this reason, the count of human-specific MAPT in each sample is calculated by the sum of FL counts for all these multiple transcripts.

########### Read in hMAPT.header from TG mice (counts of human-specific MAPT sequences)
find_mapt <- function(){
  
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


### RNASeq vs IsoSeq #################################################################

rnaseq_isoseq_transcriptome <- function(cuffrefmap_input,cufftmap_input){
  #### READ Files from Gff compare 
  # note only include isoforms from Iso-Seq that are partially or fully matched to RNA-Seq, not all isoforms from Iso-Seq dataset 
  cuffrefmap <- read.table(cuffrefmap_input, header = T)
  # no duplicates of isoform - QC checked
  # cuffrefmap[cuffrefmap$class_code == "=",] %>% group_by(ref_id) %>% tally() %>% filter(n > 1)
  
  cufftmap <- read.table(cufftmap_input, header = T)
  # Note not all isoforms are detected by RNA-Seq therefore not listed in cufftmap
  # novel_gene <- class.files[grepl("NOVE", class.files$associated_gene),"isoform"]
  # for(i in novel_gene){print(cufftmap[cufftmap $ref_id == i,])}
  
  # classification of rnaseq reads to isoseq 
  # cufftmap %>% group_by(class_code) %>% tally() %>% mutate(perc = n/sum(n) * 100) %>% ggplot(., aes(x = reorder(class_code, -perc), y = perc)) + geom_bar(stat = "identity")
  # cufftmap %>% filter(class_code %in% c("u")) %>% ggplot(.,aes(x = num_exons)) + geom_bar(aes(y = (..count..)/sum(..count..)))
  
  ## 
  # Replace id of those rnaseq isoforms that fully match with isoseq isoforms with pbid for venn diagram
  rnaseq <- c(as.character(cufftmap[cufftmap$class_code != "=","qry_id"]),as.character(cuffrefmap[cuffrefmap$class_code == "=","ref_id"])) 
  p1 <- venn.diagram(
    x = list(rnaseq, class.files$isoform), category.names = c("RNA-Seq","Iso-Seq"), filename = NULL, output=TRUE,
    lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC"), main = "\n",
    cex = 1,fontface = "bold",fontfamily = "ArialMT",
    cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-27, 27),  cat.dist = c(0.055, 0.055),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n",
    print.mode = "raw"
  )
  
  
  ## Number of Isoforms per dataset 
  num_iso <- list(class.files %>% mutate(Sample = "Iso-Seq") %>%  .[,c("isoform", "associated_gene", "novelGene","FSM_class","gene_exp","Sample")],
                  rnaseq.class.files %>%  .[,c("isoform","associated_gene", "novelGene","FSM_class","gene_exp","Sample")])
  isoPerGene <- lapply(num_iso, function(x) SQANTI_gene_preparation(x)) %>% bind_rows()
  
  # Total Number of Genes per Type 
  Total_Num <- isoPerGene %>% group_by(Sample) %>% count(Sample)
  # Number of isoform cateogories per Sample
  p2 <- isoPerGene %>% group_by(Sample) %>% count(nIsoCat) %>% full_join(Total_Num,., by = "Sample") %>% mutate(Perc = n.y/n.x * 100) %>% 
    ggplot(., aes(x=nIsoCat, fill=Sample)) +
    geom_bar(stat="identity", aes(y= Perc, group = as.factor(Sample)), color="black", size=0.3, width=0.7, 
             position="dodge") + 
    labs(x ="Number of Isoforms", y = "Genes (%)", fill = "", title = "\n") +
    #scale_fill_manual(values=c(label_colour(dataset1), label_colour(dataset2)), 
    #                 labels=c(paste(dataset1,"Cortex"), paste(dataset2,"Cortex")))  + 
    mytheme + 
    theme(legend.position = c(0.75,0.95))
  
  # Length
  p3 <- bind_rows(class.files %>% mutate(Sample = "Iso-Seq") %>% .[,c("length","Sample")], rnaseq.class.files %>% .[,c("length","Sample")]) %>% 
    ggplot(., aes(x = Sample, y = log10(length/1000))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Isoform Length (Log10 kb)", title = "\n")
  
  # Exons
  p4 <- bind_rows(class.files %>% mutate(Sample = "Iso-Seq") %>% .[,c("exons","Sample")], rnaseq.class.files %>% .[,c("exons","Sample")]) %>% 
    ggplot(., aes(x = Sample, y = log10(exons))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Number of Exons (Log10)", title = "\n")
  
  # cage peak
  p5 <- data.frame(dataset = c("Iso-Seq","RNA-Seq"),
                   within_50 = c(nrow(class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(class.files) * 100,
                                 nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(rnaseq.class.files) * 100),
                   # include NAs
                   without_50 = c(100 - nrow(class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(class.files) * 100,
                                  100 - nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(rnaseq.class.files) * 100)) %>%
    melt() %>%
    mutate(variable = factor(variable, levels = c("without_50","within_50"))) %>%
    ggplot(., aes(x = dataset, y = value, fill = variable)) + geom_bar(stat = "identity") + mytheme +
    labs(y ="Isoforms within 50bp CAGE (%)", x = "", fill = "", title = "\n") + theme(legend.position = "right") +
    scale_fill_discrete(labels = c("No","Yes"))
  
  # structural_category
  p6 <- bind_rows(class.files %>% mutate(Sample = "Iso-Seq") %>% .[,c("structural_category","Sample")] %>% 
                    group_by(structural_category,Sample) %>% tally() %>% group_by(Sample) %>% mutate(perc = n/sum(n)*100),
                  rnaseq.class.files %>% .[,c("structural_category","Sample")] %>% group_by(structural_category,Sample) %>% tally() %>% group_by(Sample) %>% mutate(perc = n/sum(n)*100)) %>% 
    ggplot(., aes(x = Sample, y = perc, fill = structural_category)) + geom_bar(stat = "identity") + mytheme +
    labs(x = "", y = "Isoforms (%)", fill = "Structural Category", title = "\n") + theme(legend.position = "bottom")
  
  output <- plot_grid(grobTree(p1),p2,p3,p4,p5,p6, labels = "auto", label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
  
  return(output)
}

# rnaseq_isoseq_counts
# correlate the number of reads from RNASeq alignment to Iso-Seq and RNA-Seq defined transcriptome using only the 23,761 isoforms considered matching from gffcompare (cuffrefmap file from gffcompare output)
rnaseq_isoseq_counts <- function(class_file){
  
  #### READ Files from Gff compare 
  # note only include isoforms from Iso-Seq that are partially or fully matched to RNA-Seq, not all isoforms from Iso-Seq dataset 
  cuffrefmap <- read.table(cuffrefmap_input, header = T)
  # no duplicates of isoform - QC checked
  # cuffrefmap[cuffrefmap$class_code == "=",] %>% group_by(ref_id) %>% tally() %>% filter(n > 1)
  
  cufftmap <- read.table(cufftmap_input, header = T)
  
  # only consider the isoforms that are "matching" between Iso-Seq and RNA-Seq defined transcriptome (output from gff compare)
  matching <- cuffrefmap[cuffrefmap$class_code == "=",] %>% mutate(qry_id = word(qry_id_list,c(2),sep = fixed("|")))
  
  # number of isoforms matching = 23,761
  # nrow(cuffrefmap[cuffrefmap$class_code == "=",]) 
  
  # ref_id = pacbio id from isoseq-defined transcriptome gtf
  matching_pbisoform <- cuffrefmap[cuffrefmap$class_code == "=","ref_id"]
  
  # qry_id_list = id from rnaseq-defined transcriptome gtf
  matching_rnaseqisoform <- unique(word(cuffrefmap[cuffrefmap$class_code == "=","qry_id_list"],c(2),sep = fixed("|")))
  
  # QC no repeated refid (PBid from isoforms)
  #n_occur <- data.frame(table(cuffrefmap[cuffrefmap$class_code == "=","ref_id"]))
  #n_occur[n_occur$Freq > 1,]
  
  # repeated qry_id due to imperfect match from cuffrefmap i.e. two PB isoforms refer to the same RNA-Seq isoform if 5' and 3' end different but internal junction the same
  #n_occur <- data.frame(table(word(cuffrefmap[cuffrefmap$class_code == "=","qry_id_list"],c(2),sep = fixed("|"))))
  #n_occur[n_occur$Freq > 1,]
  
  ### FeatureCounts
  # subset Iso-Seq Defined and RNA-Seq defined transcriptome with corresponding matching isoform id for counts  
  #IsoSeq_FL <- class_file[class_file$isoform %in% matching_pbisoform,c("isoform","FL","ISOSEQ_TPM")]
  #IsoSeq_Defmatching <- IsoSeq_Def %>% filter(Geneid %in% matching_pbisoform) %>% select(Geneid, RNA2IsoSeq.sorted.bam) %>% `colnames<-`(c("PBID", "RNA2Isoseq_Counts"))
  #RNASeq_Defmatching <- RNASeq_Def %>% filter(Geneid %in% matching_rnaseqisoform) %>% select(Geneid, RNA2RNASeq.sorted.bam) %>% `colnames<-`(c("RNASeqID", "RNA2RNAseq_Counts"))
  
  ### Kallisto
  IsoSeq_FL <- class_file[class_file$isoform %in% matching_pbisoform,c("isoform","FL","ISOSEQ_TPM")]
  IsoSeq_Defmatching <- IsoSeq_Def %>% filter(target_id %in% matching_pbisoform) %>% select(target_id, est_counts,tpm) %>% `colnames<-`(c("PBID", "RNA2Isoseq_Counts","RNA2Isoseq_TPM"))
  RNASeq_Defmatching <- RNASeq_Def %>% filter(target_id %in% matching_rnaseqisoform) %>% select(target_id, est_counts,tpm) %>% `colnames<-`(c("RNASeqID", "RNA2RNAseq_Counts","RNA2RNAseq_TPM"))
  
  
  # Merge all the counts with the isoforms that are considered matching 
  final <- merge(matching,IsoSeq_Defmatching ,by.x = "ref_id", by.y = "PBID", all = T)
  final <- merge(final,RNASeq_Defmatching ,by.x = "qry_id", by.y = "RNASeqID", all = T)
  final <- merge(final,IsoSeq_FL ,by.x = "ref_id", by.y = "isoform", all = T)
  
  # Correlation tests 
  ## Featurecounts
  #cor.test(final$RNA2Isoseq_Counts,final$RNA2RNAseq_Counts)
  #cor.test(final$RNA2Isoseq_Counts,final$FL)
  #cor.test(final$RNA2RNAseq_Counts,final$FL)
  
  cor.test(final$RNA2Isoseq_TPM,final$RNA2RNAseq_TPM, method = "pearson")
  cor.test(final$RNA2Isoseq_TPM,final$ISOSEQ_TPM, method = "pearson")
  cor.test(final$RNA2RNAseq_TPM,final$ISOSEQ_TPM,method = "pearson")
  
  
  p1 <- density_plot(final,"RNA2Isoseq_TPM","ISOSEQ_TPM", "RNA-Seq Transcript Expression from \n Iso-Seq-defined transcriptome (TPM)", "Iso-Seq Transcript Expression \n from Iso-Seq transcriptome (TPM)","")
  p2 <- density_plot(final,"RNA2Isoseq_TPM","RNA2RNAseq_TPM", "RNA-Seq Expression from \n Iso-Seq-defined transcriptome (TPM)", "RNA-Seq expression from RNA-Seq transcriptome (TPM)","")
  p3 <- density_plot(final,"RNA2RNAseq_TPM","ISOSEQ_TPM", "RNA-Seq expression from RNA-Seq transcriptome (TPM)", "Iso-Seq Transcript Expression \n from Iso-Seq transcriptome (TPM)","")
  
  return(list(p1,p2,p3))
}