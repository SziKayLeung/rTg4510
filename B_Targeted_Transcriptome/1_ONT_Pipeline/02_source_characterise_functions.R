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

find_mapt_ont <- function(type){
  if(type == "human"){
    ONT_humanMapt_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON_Human/MAPT"
  }else{
    ONT_humanMapt_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/All/HumanMapt/"
  }
  
  
  cluster_reads_B2 = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/Basecalled_trimmed_reads/Batch2_Demultiplex/Batch2_count.txt") %>% mutate(Sample = word(V2,c(1), sep = fixed("_"))) %>% group_by(Sample) %>% tally(V1) %>% mutate(BarcodedSample = paste0("Batch2",Sample))
  
  cluster_reads_B3 = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/Basecalled_trimmed_reads/Batch3_Demultiplex/Batch3_count.txt") %>% mutate(Sample = word(V2,c(1), sep = fixed("_"))) %>% group_by(Sample) %>% tally(V1) %>% mutate(BarcodedSample = paste0("Batch3",Sample))
  
  cluster_reads_input_all = rbind(cluster_reads_B2,cluster_reads_B3)
  total_cluster_reads = merge(cluster_reads_input_all,ONTBarcodedPhenotype,by = "BarcodedSample") %>% select(sample,n)
  
  ONT_humanMapt = list(paste0(ONT_humanMapt_dir,"/ONT_Batch2/hMAPT2.count.txt"),
                       paste0(ONT_humanMapt_dir,"/ONT_Batch2/mMAPT1.count.txt"),
                       paste0(ONT_humanMapt_dir,"/ONT_Batch3/hMAPT2.count.txt"),
                       paste0(ONT_humanMapt_dir,"/ONT_Batch3/mMAPT1.count.txt"))
  ONT_humanMapt = lapply(ONT_humanMapt, function(x) read.table(x))
  names(ONT_humanMapt) = c("Human_Batch2","Mouse_Batch2","Human_Batch3","Mouse_Batch3")
  ONT_humanMapt = data.frame(do.call(rbind, ONT_humanMapt)) %>% rownames_to_column(., var = "Dataset") %>% 
    mutate(mapt_specific = word(Dataset,c(1),sep = fixed("_")), 
           Batch= word(word(Dataset,c(2), sep = fixed("_")),c(1), sep = fixed(".")),
           Sample = word(V2,c(1), sep = fixed(".")),
           BarcodedSample = paste0(Batch,Sample)) 
  
  AllMapt = merge(ONT_humanMapt,ONTBarcodedPhenotype,by = "BarcodedSample") %>% dplyr::rename("sum_cluster_reads" = "V1") %>% 
    select(sample, sum_cluster_reads, Age, Genotype,mapt_specific) %>%
    full_join(., total_cluster_reads, by = "sample") %>%
    filter(!sample %in% c("C21","E18","C20","B21","AllMouseTargeted")) %>%
    # note WT would not have human-specifici MAPT reads therefore replace NA with 0
    mutate(sum_cluster_reads = replace_na(sum_cluster_reads, 0)) %>%
    # normalise with ratio
    mutate(normalised = sum_cluster_reads/n)
  
  p1 = ggplot(AllMapt, aes(x = Age, y = normalised, color = Genotype)) + geom_jitter(width = 0.09) + mytheme +
    scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) +
    labs(x = "Age (months)", y = "Ratio of species-specific \n MAPT transcripts/ total transcripts", title = "ONT") +
    #mytheme + scale_y_continuous(lim = c(0,0.06), labels = function(x) format(x, scientific = TRUE))  +
    theme(legend.position = c("bottom"), axis.text.y = element_text(angle = 90)) + 
    stat_summary(data=AllMapt, aes(x=Age, y=normalised, group=Genotype), 
                 fun="mean", geom="line", linetype = "dotted") + 
    facet_wrap(~mapt_specific) + theme(strip.background = element_blank())
  
  return(p1)
}

find_mapt_isoseq <- function(){
  humanmapt_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/All/HumanMapt/IsoSeq"
  hMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "hMAPT2.header", full.names = T)
  mMAPT_input <- list.files(path = humanmapt_input_dir, pattern = "mMAPT1.header", full.names = T)
  
  CLUSTER_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/IsoSeq/CLUSTER"
  cluster_reads <- list.files(CLUSTER_input_dir, pattern = "cluster_report.csv", full.names = T)
  
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
    names(file_input) <- lapply(file_input_names, function(x) word(word(x,c(13),  sep = fixed ('/')),c(1), sep = fixed(".")))
    
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
  plot_species_MAPT <- function(MAPT_reads,species,type){
    
    if(type == "IsoSeq"){
      MAPT_reads = merge(MAPT_reads, IsoPhenotype, by = "sample",all = T) %>% `colnames<-`(c("sample", "sum_cluster_reads", "Age","Genotype"))
      MAPT_reads = merge(MAPT_reads, total_cluster_reads, by = "sample", all = T) 
    }
    
    plot_reads <- MAPT_reads %>%
      #  do not include other J20 samples
      filter(!sample %in% c("C21","E18","C20","B21","AllMouseTargeted")) %>%
      # note WT would not have human-specifici MAPT reads therefore replace NA with 0
      mutate(sum_cluster_reads = replace_na(sum_cluster_reads, 0)) %>%
      # normalise with ratio
      mutate(normalised = sum_cluster_reads/n, mapt_specific = species)
    
    #plot_reads_old <- merge(MAPT_reads, total_cluster_reads, by = "sample", all = T) %>%
    #  # do not include other J20 samples
    #  filter(!sample %in% c("C21","E18","C20","B21","AllMouseTargeted")) %>%
    #  # classifiers of Age and Genotype
    #  mutate(Age = ifelse(grepl(paste(twomos,collapse="|"), sample),"2","8")) %>%
    #  mutate(Genotype = ifelse(grepl(paste(WT,collapse="|"), sample),"WT","TG")) %>%
    #  # note WT would not have human-specifici MAPT reads therefore replace NA with 0
    #  mutate(sum_cluster_reads = replace_na(sum_cluster_reads, 0)) %>%
    #  # normalise with ratio
    #  mutate(normalised = sum_cluster_reads/n, mapt_specific = species)
    
    mean_reads = plot_reads %>% group_by(Age,Genotype) %>%  summarise_at(vars(normalised), funs(mean(., na.rm=TRUE))) %>% 
      mutate(mapt_specific = species)
    
    return(list(plot_reads, as.data.frame(mean_reads)))
  }
  
  final_humanMAPT <- plot_species_MAPT(hMAPT_reads,"Human","IsoSeq")
  final_mouseMAPT <- plot_species_MAPT(mMAPT_reads,"Mouse","IsoSeq")
  AllMapt = rbind(final_humanMAPT[[1]],final_mouseMAPT[[1]])
  
  p1 <- ggplot(AllMapt, aes(x = Age, y = normalised, color = Genotype)) + geom_jitter(width = 0.09) + mytheme +
    scale_color_manual(values = c(label_colour("TG"),label_colour("WT")), labels = c("TG","WT")) +
    labs(x = "Age (months)", y = "Ratio of species-specific \n MAPT transcripts/ total transcripts", title = "Iso-Seq") +
    mytheme + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))  +
    theme(legend.position = c("bottom"), axis.text.y = element_text(angle = 90)) + 
    stat_summary(data=bind_rows(AllMapt), aes(x=Age, y=normalised, group=Genotype), 
                 fun="mean", geom="line", linetype = "dotted") + 
    facet_wrap(~mapt_specific) + theme(strip.background = element_blank())
  
  return(p1)
}