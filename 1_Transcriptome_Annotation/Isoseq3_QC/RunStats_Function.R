# 18/01/2020: Functions to for QC of targeted and whole transcriptome runs


# number_of_reads
# Aim: huge wrapper function to read in multiple files from CCS, LIMA and REFINE directory, and output 5 plots
number_of_reads <- function(){
  
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
  
  ####################### CLUSTER
  CLUSTER <- ldply(CLUSTER_list, function(x) nrow(x)) %>% 
      mutate(Description = "num_clustered_reads") %>% 
      mutate(sample = word(.id, c(1), sep = fixed("."))) %>% 
      select(Description, .id, V1, sample) %>%
      `colnames<-`(c("Description", "variable", "value","sample"))
      
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
                                                    "num_clustered_reads" = "Clustered reads(Transcripts)"))
  
  Reads$Description <- factor(Reads$Description , levels = c("Polymerase Reads","CCS Reads","FL Reads","FLNC reads","Poly-A FLNC reads",
                                                             "Clustered reads(Transcripts)"))
  
  ### To calculate proportions to generate plots
  # total failed ccs reads
  failed_CCS_reads <- as.data.frame(CCS_values) %>% melt(., id = "Description") %>% filter(Description %in% c("ZMWs filtered       (C)  "))
  failed_LIMA_reads <- as.data.frame(LIMA_values) %>% melt(., id = "Description") %>% filter(Description %in% c("ZMWs below any threshold  (C) "))
  
  return(Reads)
  
}

# plot of Alignable length and identity
plot_mapping_quality <- function(dat){
  p <- Merged_mapping %>% 
    mutate(identity = V8 * 100) %>% 
    mutate(length = V6 * 100) %>%
    ggplot(., aes(x = length, y = identity)) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1, name = "Density") +
    mytheme + labs(x = "Alignment Length (%)", y = "Alignment Identity (%)", title = "\n\n") +
    geom_hline(yintercept = 95, linetype = 2) +
    geom_hline(yintercept = 85, linetype = 2) +
    geom_vline(xintercept = 95, linetype = 3) + 
    theme(legend.position = "none") +
    annotate("text", x = 12, y = 97, label = "95% Threshold (Reduced)") + 
    annotate("text", x = 10, y = 87, label = "85% Threshold") + 
    annotate("text", x = 94, y = 50, label = "95% Threshold", angle = 90) +
    annotate("rect", xmin = 95, xmax = 100, ymin = 95, ymax = 100, fill = "palegreen", alpha = 0.2)
  
  return(p)
}
