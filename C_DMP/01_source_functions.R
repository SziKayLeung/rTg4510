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

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(ggvenn))
suppressMessages(library(ggplot2))
suppressMessages(library(wesanderson))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(cowplot))
suppressMessages(library(stringr))
suppressMessages(library(extrafont))
suppressMessages(loadfonts())

simple_diff_stats <- function(trans, transExp, type){
  
  df <- transExp %>% filter(isoform == trans) %>% mutate(Exp = value)
  
  meanexp = df %>% group_by(time, group) %>% summarise_at(vars(value), funs(mean(., na.rm=TRUE))) %>% as.data.frame() %>% 
    mutate(groupings = paste0(time,group)) %>% .[,c("value","groupings")] %>% spread(., groupings, value) 
  
  if(nrow(meanexp != 0)){meanexp = meanexp %>% 
    mutate(log2fc_pheno = log2(`8Case`/`8Control`),
           log2fc_age = log2(`8Case`/`2Case`))} 
  
  if(type == "isoseq"){
    meanexp = meanexp[,c(5,2,4,1,3)]
  }
  
  return(meanexp)
}

#InputGene = "Ank1"
#Norm_transcounts = annotated$ont$Norm_transcounts
plot_transexp_overtime_meth <- function(InputGene, Norm_transcounts, transsig_output, name){
  print(InputGene)
  
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  plot_title <- paste0(InputGene,"\n",name,"\n\n")
  
  transsig_output <- transsig_output[transsig_output$associated_gene == InputGene,] %>% arrange(-desc(`p-value`)) 
  sigiso = transsig_output$isoform[1:5]
  df = df %>% filter(isoform %in% sigiso)
  colours = c("#F8766D",alpha("grey",0.6))
  
  p <- ggplot(df, aes(x = time, y = value, colour = Isoform)) + geom_point(size = 3) + 
    facet_grid(~group,scales = "free", space = "free") +
    stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted", size = 1.5) +
    mytheme + labs(x = "Age (months)", y = "Normalised Counts",title = plot_title) +
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines"),
          legend.position = c(0.2,0.8)) + scale_colour_discrete(name = "Isoform")
  
  if(name == "ONT Expression"){
    p <- p + theme(legend.position="bottom", legend.direction="vertical")
  }
  
  return(p)
}


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

### Methylation #################################################################
# Identify common genes with DMP/DMR and DIU/DIE/DFI
# DMP/DMR = Differentially methylated position/differentially methylated region 
# DIU/DIE = Differential isoform usage/differential isoform expression
# Use RNASeq as expression 
# Model == Genotype or Interaction (DMP/DMR)
# Output 
# 1/ Venn Diagram of overlap of DMP with DIU/DIE/DFI based on model
# 2/ Output Table of location of overlap DMP
Methylation_Integration <- function(Model){
  # Genes identified as DFI with major isoform switching (more interesting)
  # Genes with major isoform switch using the fold change for filtering minor isoforms
  # Genes with differential transcript expression with Rsqaured > 0.7
  #DFI_yes = DFI %>% filter(DFI_Status == "DFI", MajorIsoformSwitching == "YES")
  DIU_yes = tappasDIU$rnaseq %>% filter(DIU == "DIU") %>% filter(podiumChange == "YES")
  DIE_yes = tappassigtrans$WholeRNA_Transexp[tappassigtrans$WholeRNA_Transexp$`R-squared` > 0.7,]
  
  # Select DMP/DMR model  
  if(Model == "DMP_Genotype"){
    DMP = Whole_DMP$genotype
    DMP <- DMP %>% dplyr::rename(Location = betaResultsClean....Position..)
  }else if(Model == "DMP_Interaction"){
    DMP = Whole_DMP$interaction
    DMP <- DMP %>% dplyr::rename(Location = betaResultsClean....Position..)
  }else if(
    Model == "DMR_Genotype"){
    # less stringent for DMR genotype 
    #DFI_yes = DFI %>% filter(DFI_Status == "DFI")
    DIE_yes = tappassigtrans$WholeRNA_Transexp
    DMP = DMRs
    DMP <- DMP %>% rename(SYMBOL = Chipseeker_SYMBOL)
  }else{
    print("FAIL")
  }
  
  # create venn diagram 
  x <- list(DIU = DIU_yes$gene, DIE = DIE_yes$associated_gene, `DMP Genotype` = DMP$SYMBOL)
  p = ggvenn(x, fill_color = c("#EFC000FF", "#868686FF", "#CD534CFF"), stroke_size = 0.5, set_name_size = 4)
  
  # Interesting sectors 
  #DIU_DIE_DFI_DMP = intersect(intersect(intersect(DFI_yes$Gene,DIU_yes$gene),DIE_yes$associated_gene),DMP$SYMBOL)
  DIU_DIE_DMP = intersect(intersect(DIU_yes$gene,DMP$SYMBOL),DIE_yes$associated_gene)
  #DFI_DMP = intersect(DFI_yes$Gene,DMP$SYMBOL)
  DIU_DMP = intersect(DIU_yes$gene,DMP$SYMBOL)
  DIE_DMP = intersect(DIE_yes$associated_gene,DMP$SYMBOL)
  #cat("Overlap between DFI, DIU, DIU, DMP:", DIU_DIE_DFI_DMP, "\n")
  #cat("Overlap between DFI and DMP :", DFI_DMP,"\n")
  cat("Overlap between DIU, DIE and DMP:", DIU_DIE_DMP,"\n")
  cat("Overlap between DIU and DMP:", DIU_DMP,"\n")
  cat("Overlap between DIE and DMP:", DIE_DMP,"\n")
  
  All_Overlap = unique(c(DIU_DIE_DMP,DIU_DMP,DIE_DMP))
  cat("All Genes with Overlap", All_Overlap,"\n")
  
  # DMP/DMR positions of overlap 
  Location = DMP[DMP$SYMBOL %in% All_Overlap,] %>% 
    mutate(chr = word(Location,c(1),sep = fixed(":")), start = word(Location,c(2),sep = fixed(":")), end = word(Location,c(2),sep = fixed(":")))
  
  # No FDR for DMR Genotype
  if(Model == "DMR_Genotype"){
    Location = Location %>% select(chr, start, end, SYMBOL, median.meth.group1, median.meth.group2) 
  }else{
    Location = Location %>% select(chr, start, end, SYMBOL, meth.group1.WT, meth.group2.TG,FDR_adj_genotype) 
  }
  
  output = list(p,Location)
  names(output) = c("Venn","Positions")
  return(output)
}

# subset the mean values of WT and TG across ages from the original dataset
# input: list of chromosome positions to subset: i.e. chr1:XXX
subset_mean <- function(positions){
  
  print(positions)
  
  # subset the RRBS data by the chromosome of the gene and the DMP position  
  # Use the mean of all the methylation values
  dat = RRBS_completebetas[RRBS_completebetas$position %in% positions, ] %>% reshape2::melt() %>%
    group_by(variable, position) %>% summarise_at(vars(value), list(Methylation = mean))
  
  # merge subsetted data with the phenotype data; datawrangle for plot
  df = merge(dat,RRBS_Phenotype[,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% 
    mutate(Age_months = factor(Age_months),
           Genotype = factor(Genotype, levels = c("WT","TG")), 
           Sample = word(variable,c(1),sep = fixed("_")))
  
  return(df)
}

# Plot the Methylation of DMP of interest, and the correlation of isoform expression with DMP 
# gene = common gene from the DMPgene doc
# DMPposition = position of the DMP of interest 
# isoformdiff = PB.ID of the isoform of interest 
Methylation_Integration_plots <- function(gene, isoformdiff, norm_matrix, phenotype){
  #norm_matrix=loaded$rna$input_normalized_matrix
  #phenotype=phenotype$rna
  
  # Rename column label 
  Whole_DMP$genotype <- Whole_DMP$genotype %>% dplyr::rename(Location = betaResultsClean....Position..)
  Whole_DMP$interaction <- Whole_DMP$interaction %>% dplyr::rename(Location = betaResultsClean....Position..)
  Whole_DMP$pathology <- Whole_DMP$pathology %>% dplyr::rename(Location = betaResultsClean....Position..)
  
  # Find the DMP associated to the gene from Isabel's results
  Diff_Meth = rbind(Whole_DMP$genotype[,c("Location","SYMBOL")], # Genotype DMP
                    Whole_DMP$interaction[,c("Location","SYMBOL")],
                    Whole_DMP$pathology[,c("Location","SYMBOL")]) # Interaction DMP
  
  # remove duplicated DMP (from different analysis of genotype and interaction)
  Diff_Meth = Diff_Meth[!duplicated(Diff_Meth), ]
  DMPposition = Diff_Meth %>% filter(SYMBOL == gene) 
  print(DMPposition)
  
  # subset on original data to determine mean 
  df <- subset_mean(DMPposition$Location)
  
  # plot the methylation of DMP 
  age_names <- c(`2` = "2 mos",`4` = "4 mos", `6` = "6 mos",`8` = "8 mos")
  p1 = ggplot(df,aes(x = Genotype, y = Methylation, colour = Genotype)) + geom_boxplot() + 
    geom_point(size = 3, position = position_jitter(w = 0.1, h = 0)) + 
    labs(x = " ", y = "Methylation (%)") + mytheme +
    theme(legend.position = "bottom") + scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) +  
    scale_y_continuous(labels = function(x) x*100) + theme(legend.position = "none") +
    facet_grid(~Age_months,labeller = as_labeller(age_names)) + 
    theme(strip.background = element_blank())
  
  p1b <- ggplot(df, aes(x = Genotype, y = Methylation)) + geom_boxplot() + 
    mytheme + scale_y_continuous(labels = function(x) x*100) + 
    labs(x = " ", y = "Methylation (%)") 
  
  # isoform expression 
  df2 = norm_matrix %>% 
    tibble::rownames_to_column(., "isoform") %>% filter(isoform == isoformdiff) %>% reshape2::melt(id = "isoform") %>%
    dplyr::rename(sample = variable)
  df2 = merge(df2, phenotype, by = "sample") %>% dplyr::rename(Expression = value)
  
  # merge the isoform expression and methylation 
  merged = merge(df,df2,by.x = "Sample", by.y = "sample")
  
  # plot interaction 
  p2 = ggplot(merged,aes(x = Methylation, y = Expression)) + geom_point(size = 3, aes(colour = Genotype)) + 
    labs(y = "Isoform Expression (normalised)", x = "Methylation (%)") + mytheme + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) + 
    theme(legend.position = "bottom") +  scale_x_continuous(labels = function(x) x*100) + facet_wrap(~Age_months, ncol = 1) 
  
  # correlation 
  print(cor.test(merged$Methylation,merged$Expression))
  
  # linear regression 
  lr = data.frame(Exp = merged$Expression,Meth = merged$Methylation,Genotype = merged$Genotype, Age = merged$Age_months)
  lr$Age = as.numeric(as.character(lr$Age))
  lr$Genotype = as.numeric(as.character(recode(lr$Genotype, "WT" = "0","TG" = "1")))
  lr$InteractionMethGeno = lr$Meth * lr$Genotype
  lr$InteractionMethAge = lr$Meth * lr$Age
  lr$InteractionAgeGeno = lr$Age * lr$Genotype
  lr$InteractionAgeGenoMeth = lr$Age * lr$Genotype * lr$Meth
  
  
  res = lm(Exp ~ Meth + Genotype + InteractionMethGeno + InteractionMethAge + InteractionAgeGenoMeth + Age, data = lr)
  print(summary(res))
  
  output = list(p1,p2, df, p1b)
  return(output)
}




Methylation_Integration_stats <- function(gene, sig_trans_df, exp_matrix){
  # sig_trans_df = tappassigtrans$WholeRNA_Transexp
  # exp_matrix = annotated tappas file 
  print(gene)
  dat = DMP_Combine %>% filter(SYMBOL == gene) %>% mutate(Location = paste0(chr,":", start))
  nDMP = nrow(dat)
  
  # subset on original data to determine mean 
  df <- subset_mean(dat$Location)
  
  Loc = dat$Location
  WT = signif(mean(df[df$Genotype == "WT", "Methylation"], na.rm = T),digits = 3)
  TG = signif(mean(df[df$Genotype == "TG", "Methylation"], na.rm = T),digits = 3)
  Diff = log2(TG) - log2(WT)
  FDR_meth = signif(mean(dat$FDR_adj_genotype),digits = 3)
  
  Iso = sig_trans_df[sig_trans_df$associated_gene == gene,] 
  if(nrow(Iso) >= 1){
    Iso = Iso %>% arrange(desc(`p-value`)) %>% .[1,]
    type = "DTE"
    id = Iso$isoform
    trans = Iso$associated_transcript
    R = Iso$`R-squared`
    FDR_exp = signif(Iso$`p-value`,digits = 3)
    iso_diff = simple_diff_stats(id,exp_matrix$Norm_transcounts,"rnaseq") 
    fc_age = iso_diff$log2fc_age
    fc_pheno = iso_diff$log2fc_pheno
  }else{
    Iso = tappasDIU$rnaseq[tappasDIU$rnaseq$gene == gene,]
    type = "DTU"
    id = "NA"
    trans = "NA"
    R = "NA"
    FDR_exp = signif(Iso$`adjPValue`,digits = 3)
    fc_age = "NA"
    fc_pheno = "NA"
  }
  
  output = data.frame(Gene = gene,Analysis = type,Isoform = trans,Isoform_id = id,fc_age = fc_age, 
                      fc_pheno = fc_pheno, FDR = FDR_exp,
                      R = R, Num = nDMP,location = Loc,Diff = Diff, FDR = FDR_meth)
  return(output)
}

as3mt_dmr_figure <- function(){
  # DMP position within AS3MT region
  as3mt_DMR = c("chr19:46730273","chr19:46730291","chr19:46730292","chr19:46730304","chr19:46730305","chr19:46730332","chr19:46730333","chr19:46730349","chr19:46730758")
  dat = RRBS_completebetas[RRBS_completebetas$position %in% as3mt_DMR,] 
  
  # subset from RRBS and datawrangle for plot
  subsetted_dat = dat %>% mutate(location = as.numeric(as.character(word(position,c(2),sep = fixed(":"))))) %>% reshape2::melt(id = "location") %>% merge(.,RRBS_Phenotype[,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% mutate(Methylation = as.numeric(as.character(value)) * 100) 
  
  # mean of methylation across the DMR
  df = subsetted_dat %>% filter(location != "46730758") %>% 
    group_by(Genotype,Age_months,variable) %>% summarise_at(vars(Methylation), list(Methylation = mean)) %>%
    dplyr::rename(sample = variable)
  
  df2 = loaded$rna$input_normalized_matrix %>% tibble::rownames_to_column(., "isoform") %>% 
    filter(isoform == "PB.8363.2") %>% reshape2::melt() %>% dplyr::rename(Expression = value, sample = variable) 
  
  # merge the isoform expression and methylation 
  merged = merge(df,df2, by = "sample")
  
  # correlation 
  print(cor.test(merged$Methylation,merged$Expression))
  
  # linear regression 
  lr = data.frame(Exp = merged$Expression,Meth = merged$Methylation,Genotype = merged$Genotype, Age = merged$Age_months)
  lr$Age = as.numeric(as.character(lr$Age))
  lr$Genotype = as.numeric(as.character(recode(lr$Genotype, "WT" = "0","TG" = "1")))
  lr$InteractionMethGeno = lr$Meth * lr$Genotype
  lr$InteractionMethAge = lr$Meth * lr$Age
  lr$InteractionAgeGeno = lr$Age * lr$Genotype
  
  lr$InteractionAgeGenoMeth = lr$Age * lr$Genotype * lr$Meth
  
  res = lm(Exp ~ Meth + Genotype + InteractionMethGeno + InteractionMethAge + InteractionAgeGenoMeth + Age, data = lr)
  print(summary(res))
  
  # plot
  p = ggplot(subsetted_dat,aes(x = location, y = Methylation, colour = Genotype)) + 
    annotate("rect", xmin = 46730270, xmax = 46730350, ymin = 0, ymax = 100, alpha = 0.1,fill = wes_palette("Darjeeling1")[2]) +
    geom_point(size = 2) + ylim(0,100) + xlim(46730250,46730800) + stat_summary(data=subsetted_dat, aes(x=location, y=Methylation, group=Genotype), fun ="mean", geom="line", linetype = "dotted") + scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), "Genotype", labels = c("WT","TG")) + 
    mytheme + labs(y = "Methylation (%)", x = "Chromosome 19") + theme(legend.position = "bottom")
  
  # correlation plot 
  p2 = ggplot(merged,aes(x = Methylation, y = Expression)) + 
    geom_point(size = 3, aes(colour = Genotype)) + 
    labs(y = "Isoform Expression (normalised)", x = "Methylation (%)") + mytheme + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) + 
    theme(legend.position = "bottom") + facet_wrap(~Age_months,ncol=1)
  
  
  return(list(p,p2))
}

prnp_dmr_figure <- function(){
  # DMP position within Prnp
  prnp_DMR = c("chr2:131909818","chr2:131909823","chr2:131909854","chr2:131909856","chr2:131909863","chr2:131909865","chr2:131909895",
               "chr2:131909902","chr2:131909918","chr2:131909927","chr2:131909942","chr2:131909943","chr2:131909951","chr2:131909953",
               "chr2:131909959","chr2:131909972","chr2:131909982","chr2:131909987","chr2:131909990","chr2:131910163","chr2:131910165",
               "chr2:131910181","chr2:131910202")
  
  dat = RRBS_completebetas[RRBS_completebetas$position %in% prnp_DMR,] 
  
  # subset from RRBS and datawrangle for plot
  subsetted_dat = dat %>% mutate(location = as.numeric(as.character(word(position,c(2),sep = fixed(":"))))) %>% reshape2::melt(id = "location") %>% merge(.,RRBS_Phenotype[,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% mutate(Methylation = as.numeric(as.character(value)) * 100) 
  
  
  p1 = plot_transexp_overtime("Prnp",annotated$rna$Norm_transcounts,"isoseq","RNA-Seq Expression")
  
  p2 = ggplot(subsetted_dat,aes(x = location, y = Methylation, colour = Genotype)) + 
    annotate("rect", xmin = 131909854, xmax = 131910202, ymin = 0, ymax = 100, alpha = 0.1,fill = wes_palette("Darjeeling1")[2]) +
    geom_point(size = 2) + ylim(0,100) + xlim(131909818,131910210) + stat_summary(data=subsetted_dat, aes(x=location, y=Methylation, group=Genotype), fun ="mean", geom="line", linetype = "dotted") + scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), "Genotype", labels = c("WT","TG")) + 
    mytheme + labs(y = "Methylation (%)", x = "Chromosome 2") + theme(legend.position = "bottom")
  
  p3 = plot_transexp_overtime("Prnp",annotated$iso$Norm_transcounts,"isoseq","Iso-Seq Expression")
  
  return(list(p1,p2,p3))
}
