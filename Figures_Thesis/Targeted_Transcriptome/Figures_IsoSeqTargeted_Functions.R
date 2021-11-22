# Szi Kay Leung
# Functions script for Thesis Chapter on Whole Transcriptome IsoSeq

# plot theme
loadfonts()
library("ggplot2")
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
  CCS <- read.csv(paste0(CCS_input_dir, "/Batched_CCS_output.csv"), header = T)
  LIMA <- read.csv(paste0(LiMA_input_dir, "/Batched_LIMA_summary.csv"), header = T)
  # REFINE summary input
  REFINE_json_list <- list.files(paste0(REFINE_input_dir), pattern = "flnc.filter_summary.json", full.names = T)  %>% .[!grepl("Targeted_Seq",.)]
  REFINE_list <- lapply(REFINE_json_list , function(x) as.data.frame(fromJSON(file = x)))
  names(REFINE_list) <- list.files(paste0(REFINE_input_dir), pattern = "flnc.filter_summary.json")  %>% .[!grepl("Targeted_Seq",.)]
  for(i in 1:length(names(REFINE_list))){
    REFINE_list[[i]]$Sample <- names(REFINE_list)[[i]]
  }
  # CLUSTER summary input
  CLUSTER_list_names <- list.files(paste0(CLUSTER_input_dir), pattern = ".cluster_report.csv", full.names = T) %>% .[!grepl("Targeted_Seq",.)]
  CLUSTER_list <- lapply(CLUSTER_list_names, function(x) read.csv(x))
  names(CLUSTER_list) <- list.files(paste0(CLUSTER_input_dir), pattern = ".cluster_report.csv") %>% .[!grepl("Targeted_Seq",.)]
  for(i in 1:length(names(CLUSTER_list))){
    CLUSTER_list[[i]]$Sample <- names(CLUSTER_list)[[i]]
  }
  
  ####################### CCS
  # Extract only values from mix of values and percentage
  CCS_values <- cbind(as.character(CCS[,1]), apply(CCS[,-1], 2, function(x) word(x, c(1), sep = fixed("("))))
  colnames(CCS_values)[1] <- "Description"
  CCS_values_mod <- as.data.frame(CCS_values) %>% melt(., id = "Description") %>%
    mutate(sample = word(variable, c(1), sep = fixed(".")))
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
    mutate(value = as.numeric(as.character(value))) %>%
    bind_rows(REFINE) %>%
    bind_rows(CLUSTER) %>%
    mutate(Description = as.character(Description))
  
  Reads$Description <- revalue(Reads$Description, c("ZMWs input               " ="Polymerase Reads", "ZMWs pass filters        "="CCS Reads",
                                                    "num_reads_fl"="FL Reads", "num_reads_flnc"="FLNC Reads",
                                                    "num_reads_flnc_polya" = "Poly-A FLNC Reads",
                                                    "Clustered_transcripts" = "Transcripts"))
  Reads$Description <- factor(Reads$Description, levels = c("Polymerase Reads","CCS Reads","FL Reads","FLNC Reads","Poly-A FLNC Reads","Transcripts"))

  ### To calculate proportions to generate plots
  # total failed ccs reads
  failed_CCS_reads <- as.data.frame(CCS_values) %>% melt(., id = "Description") %>% filter(Description %in% c("ZMWs filtered       (C)  "))
  failed_LIMA_reads <- as.data.frame(LIMA_values) %>% melt(., id = "Description") %>% filter(Description %in% c("ZMWs below any threshold  (C) "))
  
  
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
  #write.csv(Reads,paste0(output_helpfig_dir,"/Tg4510_IsoSeqTargetedReadsStats.csv"))
  Reads <- read.csv(paste0(output_helpfig_dir,"/Tg4510_IsoSeqTargetedReadsStats.csv"))
  Reads$Description <- factor(Reads$Description, levels = c("Polymerase Reads","CCS Reads","FL Reads","FLNC Reads","Poly-A FLNC Reads","Transcripts"))
  
  Reads_plot <- Reads %>% mutate(Batch = word(Reads$variable, c(3), sep = fixed("_"))) %>% 
    full_join(., targetedpheno, by = c("sample" = "Sample")) %>%
    unite("Batch", Batch.x,Batch.y, na.rm = TRUE) %>%
    mutate(Batch = recode(Batch, "3b" = "3", "3a" = "3 (partial run)"))
  
  p1 <- Reads_plot %>% filter(Description != "Transcripts") %>% filter(Batch != "3 (partial run)") %>%
    ggplot(., aes(x = Description, y = value, colour = Batch, group = Batch)) +
    geom_line() + geom_point() +  mytheme + theme(legend.position = c(0.8,0.8)) + labs(x = "", y = "Number of Reads (Thousands)") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_colour_discrete(name = "",labels = c("Batch 1 (n = 6)","Batch 2 (n = 9) ","Batch 3 (n = 9)"))
  
  p2 <- Reads_plot %>% filter(Description == "Transcripts") %>% 
    ggplot(., aes(x = Genotype, y = value, colour = Genotype)) +
    geom_boxplot() + geom_point() +  mytheme + 
    labs(x = "Genotype", y = "Number of FL Transcripts (Thousands)") +
    scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) + theme(legend.position = "none") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3,accuracy = 1))
  
  p3 <- Reads_plot[Reads_plot$Description == "Poly-A FLNC Reads" ,] %>% 
    ggplot(., aes(x = Batch, y = value)) + geom_boxplot() + geom_point(aes(colour = Genotype)) + mytheme +
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
  transcript_RIN <- merge(Reads_plot %>% filter(Description == "Transcripts"),tg4510_samples, by.x = "sample",by.y = "Sample.ID", all.x = T)
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
    geom_point(aes(colour = as.factor(Phenotype)), position = position_jitterdodge()) + 
    mytheme + labs(y = "On-Target Rate (%)", x = "Batch") + 
    scale_fill_manual(values = c(alpha(label_colour("TG"),0.4),alpha(label_colour("WT"),0.4)), name = "Genotype") + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), guide=FALSE) 
  
  for(i in 1:3){cat("Mean on target rate in Batch",i,":", 
                    mean(Probes[Probes$Batch == i,"perc"]),"\n")}
  return(p1)
}

level2filter <- function(){
 
  # Filter by targeted genes
  # targeted_transcripts = transcript/x of targeted genes (differented by having overlap with probes)
  # targeted_read_id = mx/x/ccs of targeted genes (origianal reads)
  # targeted_pbid = pbid of targeted genes (from readstat file after filtering for targeted_read_id)
  
  targeted <- Merged_probes_files %>% filter(Merged_probes_files$num_base_overlap != "0") %>%
    select(read_id,Gene) %>% dplyr::rename(.,transcript_id = read_id) %>%
    left_join(.,Merged_PostIso$cluster, by = c("transcript_id"= "cluster_id")) %>% 
    left_join(.,Merged_PostIso$cupcake_readstat, by = c("read_id"="id") )
  
  # Example of unmapped transcript but probed
  #Merged_probes_files %>% filter(read_id == "transcript/24")
  #Merged_PostIso$cluster %>% filter(cluster_id == "transcript/24")
  #Merged_PostIso$cupcake_readstat %>% filter(id == "m54082_200731_163617/50594276/ccs")
  
  tgcluster <- targeted %>% group_by(transcript_id,Gene) %>% tally() %>% # multiple same transcript_id from clustering ccs reads 
    group_by(Gene) %>% tally() %>% mutate(type = "Clustered")
  
  # PBids = NA for non-targeted genes
  tgcupcake <- Merged_PostIso$cupcake %>% left_join(., unique(targeted[,c("pbid","Gene")])) %>% group_by(Gene) %>% tally() %>% filter(!is.na(Gene)) %>% mutate(type = "Cupcake")
  tgcupcakefilter <- Merged_PostIso$cupcake_filter %>% left_join(., unique(targeted[,c("pbid","Gene")], by = "pbid")) %>% group_by(Gene) %>% tally() %>% filter(!is.na(Gene))  %>% mutate(type = "Cupcake Filter")
  
  # if use probes, also misannotation of other genes nearby "Gm49227"
  #tgsqfilter <- Merged_PostIso$sqanti_filter %>% left_join(., unique(targeted[,c("pbid","Gene")]), by = c("isoform" = "pbid")) %>% group_by(Gene) %>% tally() %>% filter(!is.na(Gene)) 
  tgsqfilter <- Merged_PostIso$sqanti_filter %>% filter(toupper(associated_gene) %in% targeted$Gene) %>% group_by(associated_gene) %>% tally() %>% mutate(type = "Sqanti filter") %>% mutate(Gene = toupper(associated_gene))
  
  pNumfil <- bind_rows(tgcluster, tgcupcake, tgcupcakefilter, tgsqfilter) %>%
    ggplot(., aes(x = Gene, y = n, fill = type)) + geom_bar(stat = "identity") +
    mytheme + labs(x = "", y = "Number of Transcripts (Thousand)") + theme(legend.position = c(0.8,0.8)) + 
    scale_fill_discrete(name = "", labels = c("Iso-Seq Cluster","Cupcake Collapse","Cupcake Filter","SQANTI filter")) + 
    scale_y_continuous(labels = ks) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  return(pNumfil)
}

# final_num_iso <sqanti_dataframe> <type == structural_category/associated_transcript
final_num_iso <- function(sq_df,type){
  dat = sq_df %>% filter(toupper(associated_gene) %in% TargetGene) %>% 
    select(isoform, structural_category, associated_gene, associated_transcript, subcategory) %>% 
    mutate(structural_category = factor(structural_category, levels = c("NNC","NIC","ISM","FSM"))) 
  
  if(type == "structural_category"){
    print("Structural Category")
    ### number of novel vs Known isoforms 
    dat1 <- dat %>% 
      mutate(transcript_type = ifelse(.$associated_transcript == "novel","Novel","Known")) %>% 
      mutate(transcript_type = factor(transcript_type, levels = c("Novel","Known"))) %>% 
      group_by(associated_gene,transcript_type) %>% tally() %>% spread(., transcript_type, n) %>%
      mutate(total = rowSums(.[,c("Novel","Known")])) 
    
    cat("Number of total isoforms across 20 genes:", sum(dat1$total))

    ### number of novel vs Known isoforms by subcategory
    dat2 <- dat %>% group_by(associated_gene, structural_category) %>% tally()
  }else if(type == "associated_transcript"){
    dat_at = dat %>%  mutate(Transcript_ID = ifelse(associated_transcript != "novel",associated_transcript,paste0(isoform,"_",structural_category))) %>% 
      group_by(associated_gene,Transcript_ID) %>% tally()
    
    dat1 <- dat_at %>% 
      mutate(transcript_type = ifelse(grepl("PB",Transcript_ID),"Novel","Known")) %>% 
      mutate(transcript_type = factor(transcript_type, levels = c("Novel","Known"))) %>% 
      group_by(associated_gene,transcript_type) %>% tally() %>% 
      spread(.,transcript_type,n) %>%
      mutate(total = rowSums(.[,c("Novel","Known")])) 
    
    cat("Number of total isoforms across 20 genes:", sum(dat1$total))
    
    for(i in 1:nrow(dat_at)){
      if(grepl("ENSMUST",dat_at$Transcript_ID[[i]])){
        cate = sq_df[sq_df$associated_transcript == dat_at$Transcript_ID[[i]],"structural_category"]
        dat_at$structural_category[[i]] = ifelse('FSM' %in% cate,"FSM","ISM")
      }else{
        dat_at$structural_category[[i]] = word(dat_at$Transcript_ID[[i]],c(2),sep = fixed("_"))  
      }}
      
      dat2 = dat_at %>% group_by(associated_gene, structural_category) %>% tally() 
  }else{
      print("Type Required")
    }

  
  # table for plot 
  dat_tab <- dat1 %>% mutate(Novel = paste0(Novel, " (", signif(Novel/total * 100,2),"%)"),
                 Known = paste0(Known, " (", signif(Known/total * 100,2),"%)")) %>% .[,c(1,4,2,3)] %>%
    `colnames<-`(c("Target Gene", "Total Isoforms","Novel Isoforms", "Known Isoforms"))
  
  p1 <- dat1 %>% reshape2::melt() %>%
    filter(variable != "total") %>%
    ggplot(., aes(x = reorder(associated_gene, -value), fill = variable, y = value)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "", y = "Number of Isoforms") + 
    scale_fill_discrete(name = "Isoform Classification",guide = guide_legend(reverse=TRUE)) + 
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    annotation_custom(tableGrob(dat_tab,rows=NULL, theme = ttheme_minimal(base_size = 12)),
                    xmin=25,xmax=3,ymin=65,ymax=90)
  
  ### number of novel vs Known isoforms by subcategory
  dat2_tab <- dat2 %>% spread(., structural_category, n) %>% 
    replace(., is.na(.), 0) %>%
    mutate(NIC = paste0(NIC, ifelse(NIC == "0","", paste0(" (", signif(NIC/rowSums(.[2:5]) * 100,2),"%)"))),
           FSM = paste0(FSM, ifelse(FSM == "0","", paste0(" (", signif(FSM/rowSums(.[2:5]) * 100,2),"%)"))),
           NNC = paste0(NNC, ifelse(NNC == "0","", paste0(" (", signif(NNC/rowSums(.[2:5]) * 100,2),"%)"))),
           ISM = paste0(ISM, ifelse(ISM == "0","", paste0(" (", signif(ISM/rowSums(.[2:5]) * 100,2),"%)"))))
  
  p2 <- dat2 %>% 
    mutate(structural_category = factor(structural_category, levels = c("NNC","NIC","ISM","FSM"))) %>%
    ggplot(., aes(x = reorder(associated_gene,-n), y = n, fill = structural_category)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "", y = "Number of Isoforms") +
    theme(legend.position = c(0.9,0.9)) + 
    scale_fill_manual(name = "Isoform Classification", labels = c("NNC","NIC","ISM","FSM"),
                      values = c(alpha("#F8766D",0.3),alpha("#F8766D",0.8),alpha("#00BFC4",0.3),alpha("#00BFC4",0.8)),
                      guide = guide_legend(reverse=TRUE)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "bottom")  +
    annotation_custom(tableGrob(dat2_tab,rows=NULL, theme = ttheme_minimal(base_size = 12)),
                      xmin=25,xmax=3,ymin=65,ymax=90)
  
  return(list(p1,p2))
}

whole_vs_targeted_plots <- function(){
  # from TAMA merge transcript file, filter transcripts associated with target gene 
  ADtrans <- TAMA_transfile %>% filter(toupper(gene_name) %in% TargetGene)
  # from TAMA merge transcript file, filter transcripts not associated with target gene 
  nonADtrans <- TAMA_transfile %>% filter(!toupper(gene_name) %in% TargetGene)
  
  #### plots ###
  #p1 = number of isoforms per target gene identified in whole transcriptome, targeted transcriptome, and in both
  #p2 = FL read counts in targeted transcriptome of isoforms annotated to AD genes (target) that are detected uniquely in targeted transcriptome and also in whole
  #p3 = FL read counts in whole transcriptome of isoforms not annotated to AD genes (off target) that are detected uniquely in whole and also in targeted
  
  ##### p1
  # calculate the total number of merged transcripts for plot in reorder 
  ADtrans_ttgene <- TAMA_transfile %>% group_by(gene_name) %>% tally()
  
  # Tally the merged transcripts by where it is from and also by gene name 
  # source from tama merge transcript file would either be "Targeted" or "Targeted,Whole" or "Whole"
  p1 <- ADtrans %>% group_by(sources,gene_name) %>% tally() %>%
    mutate(sources = recode(sources,"Targeted,Whole" = "Both")) %>%
    mutate(sources = factor(sources, levels = c("Whole","Both","Targeted"))) %>%
    left_join(.,ADtrans_ttgene, by = "gene_name") %>%
    ggplot(.,aes(x = reorder(gene_name,-n.y), y = n.x, fill = sources)) + geom_bar(stat = "identity") + mytheme + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "Number of Isoforms") +
    scale_fill_manual(name = "", values = c(label_colour("whole"),label_colour("whole+targeted"),label_colour("targeted")))
  
  ##### p2, p3
  # iso_filter
  # filter the transcripts in whole or targeted transcriptome to determine the difference in expression 
  # if type == ADGenes, then subset the transcripts that are in tama merge file associated with ADGenes and are either in targeted, or in targeted and equivalent whole, then find the targeted PacBio Id of that transcript, and find the FL; 
  # if type == NonADGenes, then subset the transcripts that are in tama merge file not associated with ADGenes and are either in whole, or in targeted and equivalent whole, then find the Whole PacBio Id of that transcript, and find the FL; 
  iso_filter <- function(input_dat, input_class,type){
    if(type == "ADGenes"){
      fil_group <- c("Targeted","Targeted,Whole")
      name <- "Targeted"
      colour <- scale_fill_manual(name = "", values = c(label_colour("targeted"),label_colour("whole+targeted")))
    }else if(type == "NonADGenes"){
      fil_group <- c("Whole","Targeted,Whole")
      name <- "Whole"
      colour <-  scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("whole")))
    }else{
      print("1")
    }
    
    dat <- input_dat %>% filter(sources %in% fil_group) %>% 
      # separate column by the different sources if detected in both whole and targeted
      # the order from TAMA is not consistent i.e one row would be wholeXXX,targetedXXX and another would be targetedXXX,wholeXXX
      # therefore use ifelse to grep the correct PBID reference to whole and targeted transcriptome
      separate(all_source_trans, c("match1", "match2"), ",") %>% 
      mutate(Whole = word(ifelse(word(match1, c(1), sep = fixed("_")) == "Whole",
                                 word(match1, c(2), sep = fixed("_")),
                                 word(match2, c(2), sep = fixed("_"))))) %>% 
      mutate(Targeted = word(ifelse(word(match1, c(1), sep = fixed("_")) == "Targeted",
                                    word(match1, c(2), sep = fixed("_")),
                                    word(match2, c(2), sep = fixed("_"))))) 
    
    # merge the pacbio ID from the whole or targeted transcriptome to the classification file to extract the FL reads
    dat2 <- merge(dat[,c("sources",name)],input_class[,c("isoform","FL")], by.x = name, by.y = "isoform", all.x = T) %>% 
      mutate(sources = recode(sources, "Targeted,Whole" = "Whole & Targeted", name = paste0(name," Only")))
    
    p <- ggplot(dat2, aes(x = sources, y = log10(FL), fill = sources)) + geom_boxplot() + mytheme + 
      labs(x = "Transcriptome Approach", y = paste0("FL Read Counts \n",name," Transcriptome (Log10)")) +
      colour + theme(legend.position = "none")
    
    return(p)
  }
  
  p2 <- iso_filter(ADtrans,targeted.class.files,"ADGenes")
  p3 <- iso_filter(nonADtrans,whole.class.files,"NonADGenes")
  
  return(list(p1,p2,p3))
}

# sqanti_filter_reason
# Aim: plot the number of isoforms that were removed using SQANTI filter under default settings (RNASeq coverage of 24 samples)
# Pre-requisite: targeted sqanti filtered_lite_reasons.txt, targeted sqanti pre-filtered classification files
sqanti_filter_reason <- function(){
  IP_col <- wes_palette("GrandBudapest1")[[1]]
  LC_col <- wes_palette("GrandBudapest1")[[2]]
  RT_col <- wes_palette("GrandBudapest1")[[3]]
  # plots
  # p1 = number of isoforms discarded by targeted sqanti filter (note this includes off target genes)
  # p2 = number of isoforms discarded by targeted sqanti filter, only AD-genes (target genes)
  # p3 = FL read count of isoforms that were discarded by sqanti filter by reason 
  
  p1 <- targeted_sqfil_reason %>% mutate(group = "Discarded Isoforms") %>%
    ggplot(., aes(x = group, fill = reason)) + geom_bar(stat = "count") + 
    mytheme + labs(x = "", y = "Number of Isoforms") +
    scale_fill_manual(values = c(IP_col,LC_col,RT_col), labels = c("Intrapriming","Low RNA-Seq Coverage","RT Switching"),
                      name = "Filtered Reason") + 
    ylim(0,15000) + theme(legend.position = "left")
  
  # Targeted_Pbid = the first part of PBID of target genes to subset the isoforms that are annotated to AD genes in the reason files
  Target_Pbid <- targeted.class.files %>% filter(toupper(associated_gene) %in% TargetGene) %>% mutate(PBId = paste0("PB.",word(.$isoform, c(2), sep = fixed(".")))) %>% select(PBId, associated_gene)
  
  # List of AD-associated isoforms that were removed in SQANTI filter 
  filtered_Target <- targeted_sqfil_reason %>% filter(PBId %in% unique(Target_Pbid[[1]])) %>% left_join(.,Target_Pbid) 
  
  p2 <- filtered_Target %>% group_by(associated_gene, reason) %>% tally() %>% 
    ggplot(., aes(x = reorder(associated_gene,-n), fill = reason, y = n)) + geom_bar(stat = "identity") + mytheme + 
    labs(y = "Number of Isoforms discarded \n by SQANTI filter", x = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual(values = c(IP_col,LC_col,RT_col), labels = c("Intrapriming","Low RNA-Seq Coverage","RT Switching"),
                      name = "Filtered Reason") + theme(legend.position = "none")
  
  p3 <- targeted.preclass.files %>% filter(toupper(associated_gene) %in% TargetGene) %>% 
    left_join(., targeted_sqfil_reason, by = c("isoform" = "filtered_isoform")) %>% 
    filter(reason != "NA") %>%
    mutate(reason = recode(reason, "LowCoverage/Non-Canonical" = "Low RNA-Seq Coverage", "IntraPriming" = "Intrapriming",
                           "RTSwitching" = "RT Switching")) %>%
    ggplot(., aes(x = reason, y = log10(FL), fill = associated_gene)) + geom_boxplot() + mytheme + 
    labs(x = "", y = "FL Read Count of discarded isoforms (Log10)") + 
    theme(legend.position = "bottom") + scale_fill_discrete(name = "Target Genes")
  
  return(list(p1,p2,p3))
}

# summary of the number of isoforms per target gene and the list of known isoforms
summary_isoform <- function(input_sqanti_file,gene){
  # subset gene from input_sqanti_file
  dat = input_sqanti_file[toupper(input_sqanti_file$associated_gene) == gene,]
  
  # list of known isoforms 
  known_isoforms = paste0(unique(dat[dat$associated_transcript != "novel","associated_transcript"]),sep = ",", collapse = '')
  
  # new data frame 
  df <- data.frame(
    Target_Gene = dat$associated_gene[1],
    Total_Isoforms_Num = length(dat$isoform),
    Known_Isoforms_Num = nrow(dat[dat$associated_transcript != "novel",]), 
    NIC = nrow(dat[dat$structural_category == "NIC",]),
    NNC = nrow(dat[dat$structural_category == "NNC",]), 
    Known_Isoforms = str_sub(known_isoforms,1,nchar(known_isoforms)-1) # remove the final comma
    
  )
  return(df)
}

### Human MAPT #################################################################

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

#### Validation of SQANTI filter 
# Aim: Work off the preSqanti filter data for targeted transcriptome, but first filter out for partial degraded transcripts using TAMA filtering script 
# applied TAMA filtering of partial transcripts to the prefiltered SQANTI filter data 
sqantifilter_valdidation <- function(){
  
  # Read in the isoforms that were retained after applying tama filtering to sqanti prefiltered data 
  tamaretained = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI2/TAMAFILTERSQANTI/tama_retained_pbid.txt")
  
  # list of isoforms that were filtered by SQANTI as low coverage, annotated to Target Genes
  lowcov_isoforms = targeted_sqfil_reason[targeted_sqfil_reason$reason == "LowCoverage/Non-Canonical","filtered_isoform"] # refer to all isoforms inc not annotated to target genes
  unfil_isoforms = targeted.preclass.files[toupper(targeted.preclass.files$associated_gene) %in% TargetGene,"isoform"] # all isoforms annotated to target genes including unfiltered isoforms
  AD_lowcov_isoforms = intersect(lowcov_isoforms,unfil_isoforms)
  
  # list of isoforms that were not filtered by SQANTI, annotated to Target Genes 
  AD_highcov_isoforms = targeted.class.files[toupper(targeted.class.files$associated_gene) %in% TargetGene,"isoform"]
  
  # All isoforms (low cov and high cov)
  all_AD_isoforms = bind_rows(data.frame(isoform = AD_highcov_isoforms) %>% mutate(type = "rnaseq_supported"),
                              data.frame(isoform = AD_lowcov_isoforms) %>% mutate(type = "not_rnaseq_supported"))
  
  # no common isoforms that were removed by SQANTI but kept by TAMA filter 
  #cat("Isoforms that were filtered by SQANTI but kept in by TAMA")
  intersect(tamaretained$V1, lowcov_isoforms)
  
  # merge the final sqanti classification file after filtering for RT switching and intrapriming 
  # note many isoforms were removed for RT switching, intrapriming and partial degradation
  final = merge(targeted.preclass.files[,c("isoform","associated_gene","associated_transcript")], all_AD_isoforms, all.y = T)
  
  p <- final %>% group_by(associated_gene,type) %>% tally %>%
    filter(associated_gene != "NA") %>%
    ggplot(.,aes(x = reorder(associated_gene,-n), y = n, fill = type)) + geom_bar(stat = "identity") + 
    labs(x = "", y = "Number of Novel Isoforms") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    mytheme + theme(legend.position = c(0.8,0.9)) + scale_fill_discrete(labels = c("No","Yes"), name = "Supported by RNA-Seq")
  
  # list of isoforms that were filtered for RT switching/intrapriming 
  IR_isoforms = targeted_sqfil_reason[targeted_sqfil_reason$reason == "IntraPriming","filtered_isoform"] # refer to all isoforms inc not annotated to target genes
  AD_IR_isoforms = intersect(IR_isoforms,unfil_isoforms)
  RTS_isoforms = targeted_sqfil_reason[targeted_sqfil_reason$reason == "RTSwitching","filtered_isoform"] # refer to all isoforms inc not annotated to target genes
  AD_RTS_isoforms = intersect(RTS_isoforms,unfil_isoforms)
  
  # list of partially degraded isoforms associated to Target Genes i.e not Intrapriming isoforms, not RT Switching isoforms, not retained isoforms
  AD_deg_isoforms = setdiff(unfil_isoforms,c(as.character(tamaretained$V1),as.character(AD_IR_isoforms),as.character(RTS_isoforms), as.character(final$isoform)))
  
  AD_artefacts = bind_rows(data.frame(isoform = AD_IR_isoforms) %>% mutate(`Filtered Reason` = "IntraPriming"),
                           data.frame(isoform = AD_RTS_isoforms) %>% mutate(`Filtered Reason` = "RT Switching"),
                           data.frame(isoform = AD_deg_isoforms) %>% mutate(`Filtered Reason` = "Partial Degradation"))
  
  # merge the final sqanti classification file after filtering for RT switching and intrapriming 
  AD_artefacts_des = merge(targeted.preclass.files[,c("isoform","associated_gene","associated_transcript")], AD_artefacts, all.y = T)
  
  p2 <- AD_artefacts_des %>% group_by(associated_gene,`Filtered Reason`) %>% tally %>%
    filter(associated_gene != "NA") %>%
    ggplot(.,aes(x = reorder(associated_gene,-n), y = n, fill = `Filtered Reason`)) + geom_bar(stat = "identity") + 
    labs(x = "", y = "Number of Isoforms") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + mytheme + theme(legend.position = c(0.8,0.9)) 
  
  cat("Total number of isoforms pre-sqanti filtering:",nrow(targeted.preclass.files),"\n")
  cat("Number of isoforms removed due to intrapriming:",length(AD_IR_isoforms),"\n")
  cat("Number of isoforms removed due to reverse transcription switching:",length(AD_RTS_isoforms),"\n")
  cat("Number of isoforms removed due to partial degradation:",length(AD_deg_isoforms),"\n")
  cat("Total number of isoforms after sqanti and tama filtering (not accounting for rnaseq support):",nrow(final),"\n")
  cat("Total number of isoforms after sqanti and tama filtering (rnaseq support):",nrow(final %>% filter(type == "rnaseq_supported")),"\n")
  cat("Total number of isoforms after sqanti and tama filtering (not rnaseq support):",nrow(final %>% filter(type == "not_rnaseq_supported")),"\n")
  cat("Total number of novel isoforms:", nrow(final %>% filter(associated_transcript == "novel")))
  
  # checking that the numbers make sense 
  #for(gene in TargetGene){
  #  dat3 = AD_artefacts_des %>% group_by(associated_gene,`Filtered Reason`) %>% tally 
  #  arten = sum(dat3[toupper(dat3$associated_gene) == gene,"n"])
  #  totaln = targeted.preclass.files[toupper(targeted.preclass.files$associated_gene) == gene,] %>% nrow()
  #  finaln = final[toupper(final$associated_gene) == gene,] %>% nrow()
  #  cat("Processing",gene,"***********\n")
  #  if(totaln == arten + finaln){print("TRUE")}
  #}
  
  return(list(p,p2))
}


### TappAS Differential Analysis #################################################################
input_tappasfiles <- function(tappas_input_dir){
  # read in files generated from TAPPAS
  tappasfiles <- list.files(path = tappas_input_dir, pattern = ".tsv", full.names = T)
  tappasfiles <- lapply(tappasfiles, function(x) read.table(x, sep = "\t", header = T))
  names(tappasfiles) <- list.files(path = tappas_input_dir, pattern = ".tsv")
  
  # InputExpression Matrix from TAPPAS of which transcripts are filtered during normalisation 
  # match the filtered transcripts with the associated gene uing the isoform id from classification file
  tappasfiles$tappAS_Transcripts_InputExpressionMatrix.tsv <- 
    merge(targeted.class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], 
          tappasfiles$tappAS_Transcripts_InputExpressionMatrix.tsv, by.x = "isoform", by.y = "Id")
  
  return(tappasfiles)
}

tappas_removediso <- function(filteredtappasfile){
  # number of transcripts that are filtered for statistical purposes
  dat = filteredtappasfile
  
  # tally of the number of transcripts filtered per gene 
  dat2 = dat %>% group_by(associated_gene, structural_category, Filtered) %>% tally()
  
  # plot the number of transcripts by gene only 
  p1 <- ggplot(dat2, aes(x = associated_gene, y = n, fill = Filtered)) + 
    geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms") + mytheme + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual(labels = c("Retained","Removed due to low coverage"), values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[5])) + 
    theme(legend.position = c(0.85,0.8))
  
  p2 <- dat2 %>% filter(Filtered != "NO") %>% 
    ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + 
    geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms Removed") + 
    mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = c(0.8,0.8))
  
  return(list(p1,p2))
  
}

# tappas_resultsanno
# prerequisite for plots
# aim1: annotate the tappas normalised expression file with the isoforms from the class file and join by phenotype
# aim2: deduce the gene expression by the sum of the expression of the filtered isoforms 
# output: norm_transcounts (normalised transcript counts) & GeneExp
tappas_resultsanno <- function(classification_file, tappas_normalised_expmatrix,phenotype){
  
  # Annotate the normalised expression matrix from tappas with the associated gene and sample phenotypes for plots 
  Norm_transcounts = 
    # annotate the tappas output normalised expression matrix of transcripts by the associated gene name and transcript name
    merge(classification_file [,c("isoform","associated_gene","associated_transcript","structural_category")], 
                            tappas_normalised_expmatrix, by.x = "isoform", by.y = 0) %>% reshape2::melt() %>% 
    # annotate the samples to the phenotype 
    left_join(., phenotype, by = c("variable" = "sample")) %>% 
    # change factor levels for plots
    mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("WT", "TG")),
           structural_category=recode(structural_category, `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM",`novel_in_catalog`="NIC"),
           Isoform = paste0(isoform,"_",structural_category))
  
  # Deduce gene expression from the sum of normalised transcript counts 
  GeneExp = Norm_transcounts %>% group_by(associated_gene,variable) %>% dplyr::summarise(Exp = sum(value)) %>%
    left_join(., phenotype, by = c("variable" = "sample"))
  
  output <- list(Norm_transcounts,GeneExp)
  names(output) <- c("Norm_transcounts","GeneExp")
  return(output)
}


# plot gene expression or by isoformID 
plot_mergedexp <- function(InputGene,IsoformID,GeneExp,Norm_transcounts){
  if (InputGene != "NA"){
    df <- GeneExp %>% filter(associated_gene == InputGene) 
    df$group <- factor(df$group, levels = c("CONTROL","CASE"))
    plot_title <- InputGene
  }else if(IsoformID != "NA"){
    df <- Norm_transcounts  %>% filter(isoform == IsoformID) %>% left_join(., phenotype, by = c("variable" = "sample"))
    colnames(df)[6] <- "Exp"
    df$group <- factor(df$group, levels = c("CONTROL","CASE"))
    plot_title <- paste0(df$associated_gene,": ",df$associated_transcript)
  }else{
    print("2 variables required")
  }
  
  p <- ggplot(df, aes(x = time, y = Exp, colour = group)) + geom_point() + 
    stat_summary(data=df, aes(x=time, y=Exp, group=group), fun="mean", geom="line", linetype = "dotted") + 
    labs(y = "Normalised Gene Expression", x = "Age (months)", title = paste0(plot_title,"\n\n")) + mytheme +
    # colours in the opposite direction of the levels 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) + 
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"))
  
  # subtitles
  if(IsoformID != "NA"){
    p <- p + labs(title = plot_title, subtitle = paste0(df$isoform,"\n\n"), y = "Normalised Isoform Expression") + theme(plot.subtitle = element_text(hjust = 0.5, size = 12,face = "italic"))
  }
  
  return(p)
}

# group_plots 
group_plots <- function(genegroup,plottype){
  myplots <- list()
  myplots2 <- list()
  
  if(length(genegroup) > 6){
    for(i in 1:6){myplots[[i]] <- plottype[[genegroup[[i]]]]}
    output1 = plot_grid(plotlist=myplots,labels = paste(letters[1:6]),label_size = 30, label_fontfamily = "CM Roman", nrow = 3,ncol = 2, scale = 0.9)
    count = 1
    for(i in 7:length(genegroup)){myplots2[[count]] <- plottype[[genegroup[[i]]]]; count = count + 1}
    output2 = plot_grid(plotlist=myplots2,labels = paste(letters[6:length(genegroup)]),label_size = 30, label_fontfamily = "CM Roman", nrow = 3,ncol = 2, scale = 0.9)
    output = list(output1,output2)
  }else{
    for(i in 1:length(genegroup)){myplots[[i]] <- plottype[[genegroup[[i]]]]}
    output = plot_grid(plotlist=myplots,labels = paste(letters[1:length(genegroup)]),label_size = 30, label_fontfamily = "CM Roman", nrow = 3,ncol = 2, scale = 0.9)
  }
  return(output)
}

group_plots_rnavsiso <- function(gene, plottype1, plottype2){
  myplots <- list()
  myplots[[1]] <- plottype1[[gene]] 
  myplots[[2]] <- plottype2[[gene]]
  output = plot_grid(plotlist=myplots,labels = c("a","b"),label_size = 30, label_fontfamily = "CM Roman",scale = 0.9)
  return(output)
}

# Differentially expressed genes 
# Input: tappassig with the sheet names "WholeIso_Geneexp" and "WholeRNA_Geneexp"
# Plots: 
# P1: venn diagram of genes that are differentially expressed between RNA+RNA(Isabel),Iso+RNA,Iso+Iso
tappas_genesig <- function(){
  
  # plot of the significant genes from Targeted Transcriptome using either Iso-Seq or RNA-Seq as expression input 
  p1 <- bind_rows(tappassig$TargetedIso_Genexp %>% mutate(type = "TargetedIso"),
                  tappassig$TargetedRNA_Genexp %>% mutate(type = "TargetedRNA")) %>% 
    ggplot(., aes(x = `...1`, y = `R-squared`, fill = type)) + geom_bar(position = position_dodge(preserve = "single"), stat = "identity") + mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme(legend.position = "bottom") + labs(x = "", y =  expression(paste(R^2))) + 
    geom_hline(yintercept=0.5,linetype="dashed") + scale_fill_manual(name = "Transcript Expression Input", labels = c("Iso-Seq","RNA-Seq"), values = c(label_colour("isoseq"),label_colour("rnaseq")))
  
  return(p1)
}

plot_transexp <- function(InputGene,Norm_transcounts,type,name){
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  df$time <- as.factor(df$time)
  
  #p <-
  #  ggplot(df, aes(x = reorder(isoform,-value), y = value, colour = group, shape = time)) + #geom_boxplot() + 
  #  geom_jitter(size = 3, position = position_jitterdodge()) +
  #   stat_summary(data=df, aes(x=group, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
  #  mytheme + labs(x = "", y = "Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
  #  theme(strip.background = element_blank(), 
  #        plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
  #        panel.spacing = unit(2, "lines"), 
  #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #  legend_theme +
  #  facet_grid(~structural_category,  scales = "free", space = "free") +
  #  scale_y_continuous(trans='log10') +
  #  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  
  p <- df %>% mutate(grouping = paste0(group,"_",time)) %>% group_by(grouping, isoform) %>% dplyr::summarise(structural_category,time,group, mean_exp = mean(value), .groups = 'drop') %>% mutate(group = factor(group, levels = c("TG","WT"))) %>% 
    ggplot(., aes(x = reorder(isoform,-mean_exp), y = mean_exp)) + geom_point(aes(colour = group, shape = time),size = 3) + mytheme + 
    labs(x = "", y = "Mean Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
    theme(strip.background = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    legend_theme +
    facet_grid(~structural_category,  scales = "free", space = "free") +
    scale_y_continuous(trans='log10') +
    guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  
  
  if(type == "isoseq"){
    p <- p + scale_shape_manual(name = "Age", values=c(1, 16), label = c("2 months", "8 months")) + 
      scale_colour_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT"))) 
  } else {
    p <- p + scale_shape_manual(name = "Age (months)", values=c(1, 16, 2, 17)) + 
      scale_colour_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT"))) 
  }
  
  return(p)
}

plot_transexp_overtime <- function(InputGene,Norm_transcounts){
  
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  
  # linear regression for significant changes over time
  #for(isoform in unique(df$isoform)){
  #  cat("Processing",isoform,"\n")
  #  df1 <- df[df$isoform == isoform & df$group == "WT",]
  #  df1_WTmean <- df1 %>% group_by(time) %>% summarise(mean_exp = mean(value))
  #  df1_TG <- df[df$isoform == isoform & df$group == "TG",]
  #  df2 <- merge(df1_TG, df1_WTmean, by = "time") %>% mutate(diff = abs(value - mean_exp))
  #print(summary(lm(diff~0 + time,df2)))
  #}
  
  #df_WT <-  df[df$group == "WT",]
  #df_WTmean <- df %>% group_by(time, isoform) %>% summarise(mean_exp = mean(value),  .groups = 'drop')
  #df_TG <- df[df$group == "TG",]
  #df2 <- merge(df_TG, df_WTmean, by = c("time","isoform")) %>% mutate(diff = abs(value - mean_exp))
  
  #p1 <- ggplot(df2, aes(x = time, y = diff, colour = isoform)) + geom_point() +
  #  stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
  #  scale_y_continuous(trans = 'log10') + mytheme + labs(x = "Age (months)", y = "Fold Change of Isoform Expression \n (TG - WT)") + theme(legend.position = "right")
  
  p <- ggplot(df, aes(x = time, y = value, colour = Isoform)) + geom_point() + 
    facet_grid(~group,scales = "free", space = "free") +
    stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
    mytheme + labs(x = "Age (months)", y = "Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
    theme(strip.background = element_blank(), legend.position = "right",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines")) 
  return(p)
}


