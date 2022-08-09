# Szi Kay Leung
# Functions script for whole differential analysis 

# packages
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
suppressMessages(library(data.table))
suppressMessages(library(readxl))
suppressMessages(library("xlsx"))
suppressMessages(library(ggvenn))
suppressMessages(library(pheatmap))

detach("package:plyr")
suppressMessages(library(extrafont))
#font_install('fontcm')
suppressMessages(loadfonts())

# do not output log files for venn diagrams
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis/GeneralPlots.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis/DifferentialPlots.R")

#

# segregate tappas output into the different models by filtering on the regression coefficients
segregate_tappasresults <- function(output_tappas, type){
  # tappas ouput with p values for each regression coefficient
  dat = output_tappas
  
  # apply Rsquared threshold at 0.5 
  dat = dat %>% filter(`p-value` < 0.05) %>% filter(`R-squared` > 0.5)
  
  if(type == "IsoSeq"){
    ### models
    # 1: casevscontrol yes, time no, timexcase no = genotype effect
    # 2. casevscontrol yes, time yes, timexcase no = genotype+age effect
    # 3. casevscontrol no, time yes, timexcase no = age effect
    # 4. casevscontrol no, time yes, timexcase yes = interaction effect
    # 5. casevscontrol no, time no, timexcase yes = interaction effect
    # 6. casevscontrol yes, time no, timexcase yes = interaction effect
    # 7. casevscontrol yes, time yes, timexcase yes = interaction effect
    ## yes = significant regression coefficient i.e != NA
    
    models = list(
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "TRUE" & is.na(p.valor_TimexCASE) == "TRUE"), # model 1
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"), # model 2
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"), # model 3
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "FALSE"), # model 4
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "TRUE" & is.na(p.valor_TimexCASE) == "FALSE"), # model 5
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "TRUE" & is.na(p.valor_TimexCASE) == "FALSE"), # model 6
      dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "FALSE") # model 7
    )
    names(models) = c("Model 1 Genotype","Model 2 Genotype+Age","Model 3 Age","Model 4 Interaction", "Model 5 Interaction", "Model 6 Interaction","Model 7 Interaction")
    models = lapply(models, function(x) x %>% arrange(-desc(`p-value`)))
    
    
  }else if(type == "RNASeq"){
    # Model 1
    model1 = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "TRUE" & 
                              is.na(p.valor_TimexCASE) == "TRUE" & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") 
    
    # Model 2
    model2a = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"
                             & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") # model 2
    
    model2b = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time2) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"
                             & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") # model 2
    
    model2c = dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "FALSE" & is.na(p.valor_Time3) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE"
                             & is.na(p.valor_Time2xCASE) == "TRUE" & is.na(p.valor_Time3xCASE) == "TRUE") # model 2
    
    model2 = bind_rows(list(model2a,model2b,model2c)) %>% distinct()
    
    # Model 3
    model3a =  dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE", 
                              is.na(p.valor_Time2xCASE) == "TRUE", is.na(p.valor_Time3xCASE) == "TRUE")
    model3b =  dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time2) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE",
                              is.na(p.valor_Time2xCASE) == "TRUE", is.na(p.valor_Time3xCASE) == "TRUE")
    model3c =  dat %>% filter(is.na(p.valor_CASEvsCONTROL) == "TRUE" & is.na(p.valor_Time3) == "FALSE" & is.na(p.valor_TimexCASE) == "TRUE",
                              is.na(p.valor_Time2xCASE) == "TRUE", is.na(p.valor_Time3xCASE) == "TRUE")
    model3 = bind_rows(list(model3a,model3b,model3c)) %>% distinct()
    
    non_interaction = bind_rows(model1,model2,model3)
    interaction = dat[dat$...1 %!in% non_interaction$...1,]
    
    models = list(model1,model2,model3,interaction)
    names(models) = c("Model 1 Genotype","Model 2 Genotype+Age","Model 3 Age","Model 4 - 7 Interaction")
    models = lapply(models, function(x) x %>% arrange(-desc(`p-value`)))
    
  }else{
    print("IsoSeq or RNASeq required as argument")
  }
  
  p <- sapply(models, nrow) %>% reshape2::melt() %>% rownames_to_column(var = "model") %>% mutate(num = c(1:nrow(.))) %>%
    ggplot(., aes(x = reorder(model, -num), y = value)) + geom_bar(stat = "identity") +
    labs(y = "Number of Differentially Expressed Genes",x = "") + mytheme + coord_flip() 
  
  df = sapply(models, nrow) %>% reshape2::melt() %>% as.data.frame() %>% rownames_to_column(var = "Model") %>% mutate(value = as.numeric(as.character(value)))
  print(df)
  cat("Total of Number differentially expressed genes:", sum(df$value),"\n")
  cat("Number differentially expressed under Interaction effect:", sum(df[grepl("Interaction",df$Model),"value"]),"\n")
  cat("Number differentially expressed under Genotype effect:", sum(df[df$Model == "Model 1 Genotype","value"]),"\n")
  
  output = list(models,p)
  names(output) = c("models","p")
  return(output)
}

### General Description #################################################################
summary_dataset <- function(){
  output = data.frame()
  count = 1
  for(d in unique(all.group.class.files$Dataset)){
    dat = all.group.class.files %>% filter(Dataset == d)
    
    total_genes = length(unique(dat$associated_gene)) 
    annotated_genes = length(unique(dat[!grepl("novelGene",dat$associated_gene),"associated_gene"]))
    novel_genes = length(unique(dat[grepl("novelGene",dat$associated_gene),"associated_gene"]))
    
    # number of unique genes, novel genes, annotated genes, isoforms
    output[1,count] =  total_genes
    output[2,count] =  paste0(annotated_genes," (", round(annotated_genes/total_genes*100,2),"%)") 
    output[3,count] =  paste0(novel_genes," (", round(novel_genes/total_genes*100,2),"%)") 
    output[4,count] =  nrow(dat)                          
    
    # 9 levels of structural cateogory 
    struct <- vector("numeric", 9)
    for(num in 1:length(levels(dat$structural_category))){
      cate <- nrow(dat[dat$structural_category == levels(dat$structural_category)[num],])
      struct[num] <- paste0(cate, " (", round(cate/dim(dat)[1]*100,2),"%)") 
    }
    output[5:13,count] <- struct
    
    # median,min,max length of isoforms, median number of exons
    output[14,count] <- paste0("Median: ", median(dat$length),", Range: ",min(dat$length), "-", max(dat$length))
    output[15,count] <- paste0("Median: ", median(dat$exons),", Range: ",min(dat$exons), "-", max(dat$exons))
    
    withincage = nrow(dat[abs(dat$dist_to_cage_peak) <= 50,])
    output[16,count] <- paste0(withincage," (", round(withincage/dim(dat)[1]*100,2),"%)")
    
    output[17,count] <- nrow(dat[dat$structural_category == "FSM",])
    
    colnames(output)[count] <- d
    count = count + 1
  }
  
  row.names(output) <- c("Total Number of Genes", "Annotated Genes","Novel Genes","Total Number of Isoforms","FSM",
                         "ISM","NIC","NNC","Genic\nGenomic","Antisense","Fusion","Intergenic","Genic\nIntron","Isoform Length (bp)","Number of Exons","Number of Ifsoforms with 50bp CAGE","FSM_count")
  
  
  output <- output[,c("WT_sq","TG_sq","WT_2m","WT_8m","TG_2m","TG_8m")]
  colnames(output) <- c("WT","TG","WT 2months","WT 8months", "TG 2months","TG 8months")
  return(output)
}

#### WGCNA
venn_WGCNA <- function(){
  genesigs = c(gene_sigs_WholeIso_lst$models$`Model 1 Genotype`$...1,gene_sigs_WholeIso_lst$models$`Model 2 Genotype+Age`$...1,
               gene_sigs_WholeIso_lst$models$`Model 4 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 5 Interaction`$...1,
               gene_sigs_WholeIso_lst$models$`Model 6 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 7 Interaction`$...1)
  all_modules = c(WGCNA$`Turquoise `$Gene,WGCNA$Red$Gene,WGCNA$Yellow$Gene)
  
  # RNA-Seq genes
  #Turquoise = length(intersect(gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1,WGCNA$`Turquoise `$Gene))
  #Red = length(intersect(gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1,WGCNA$Red$Gene))
  #Yellow = length(intersect(gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1,WGCNA$Yellow$Gene))
  #None = length(setdiff(gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1,all_modules))
  
  # Iso-Seq genes
  Turquoise = length(intersect(genesigs,WGCNA$`Turquoise `$Gene))
  Red = length(intersect(genesigs,WGCNA$Red$Gene))
  Yellow = length(intersect(genesigs,WGCNA$Yellow$Gene))
  None = length(setdiff(genesigs,all_modules))
  
  cat(length(genesigs))
  
  
  data <- data.frame(
    category=c("Turquoise Module", "Red Module", "Yellow Module","Not in Module"),
    count=c(Turquoise, Red, Yellow,None)
  )
  
  
  # Compute percentages
  data$fraction <- data$count / sum(data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax <- cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin <- c(0, head(data$ymax, n=-1))
  
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) /2
  
  # Compute a good label
  data$label <- paste0(data$count)
  
  # Make the plot
  p = ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
    scale_fill_manual(values = c(wes_palette("Royal1")[[1]],wes_palette("Zissou1")[[5]],
                                 wes_palette("Zissou1")[[1]],wes_palette("Zissou1")[3])) +
    scale_colour_manual(values = c(wes_palette("Royal1")[[1]],wes_palette("Zissou1")[[5]],
                                   wes_palette("Zissou1")[[1]],wes_palette("Zissou1")[3])) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_blank(),legend.text = element_text(size=16),legend.direction = "vertical")
  
  return(p)
  
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
Methylation_Integration<- function(Model){
  # Genes identified as DFI with major isoform switching (more interesting)
  # Genes with major isoform switch using the fold change for filtering minor isoforms
  # Genes with differential transcript expression with Rsqaured > 0.7
  #DFI_yes = DFI %>% filter(DFI_Status == "DFI", MajorIsoformSwitching == "YES")
  DIU_yes = tappasDIU$rnaseq %>% filter(DIU == "DIU") %>% filter(podiumChange == "YES")
  DIE_yes = tappassigtrans$WholeRNA_Transexp[tappassigtrans$WholeRNA_Transexp$`R-squared` > 0.7,]
  
  # Select DMP/DMR model  
  if(Model == "DMP_Genotype"){DMP = Whole_DMP$genotype}else if(
    Model == "DMP_Interaction"){DMP = Whole_DMP$interaction}else if(
      Model == "DMR_Genotype"){
      # less stringent for DMR genotype 
      #DFI_yes = DFI %>% filter(DFI_Status == "DFI")
      DIE_yes = tappassigtrans$WholeRNA_Transexp
      DMP = DMRs}else{print("FAIL")}
  
  # create venn diagram 
  x <- list(DIU = DIU_yes$gene, DIE = DIE_yes$associated_gene, `DMP Genotype` = DMP$Chipseeker_SYMBOL)
  p = ggvenn(x, fill_color = c("#EFC000FF", "#868686FF", "#CD534CFF"), stroke_size = 0.5, set_name_size = 4)
  
  # Interesting sectors 
  #DIU_DIE_DFI_DMP = intersect(intersect(intersect(DFI_yes$Gene,DIU_yes$gene),DIE_yes$associated_gene),DMP$Chipseeker_SYMBOL)
  DIU_DIE_DMP = intersect(intersect(DIU_yes$gene,DMP$Chipseeker_SYMBOL),DIE_yes$associated_gene)
  #DFI_DMP = intersect(DFI_yes$Gene,DMP$Chipseeker_SYMBOL)
  DIU_DMP = intersect(DIU_yes$gene,DMP$Chipseeker_SYMBOL)
  DIE_DMP = intersect(DIE_yes$associated_gene,DMP$Chipseeker_SYMBOL)
  #cat("Overlap between DFI, DIU, DIU, DMP:", DIU_DIE_DFI_DMP, "\n")
  #cat("Overlap between DFI and DMP :", DFI_DMP,"\n")
  cat("Overlap between DIU, DIE and DMP:", DIU_DIE_DMP,"\n")
  cat("Overlap between DIU and DMP:", DIU_DMP,"\n")
  cat("Overlap between DIE and DMP:", DIE_DMP,"\n")
  
  All_Overlap = unique(c(DIU_DIE_DMP,DIU_DMP,DIE_DMP))
  cat("All Genes with Overlap", All_Overlap,"\n")
  
  # DMP/DMR positions of overlap 
  Location = DMP[DMP$Chipseeker_SYMBOL %in% All_Overlap,] %>% 
    mutate(chr = word(Location,c(1),sep = fixed(":")), start = word(Location,c(2),sep = fixed(":")), end = word(Location,c(2),sep = fixed(":")))
  
  # No FDR for DMR Genotype
  if(Model == "DMR_Genotype"){
    Location = Location %>% select(chr, start, end, Chipseeker_SYMBOL, Chipseeker_annotation, meth.group1.WT, meth.group2.TG) 
  }else{
    Location = Location %>% select(chr, start, end, Chipseeker_SYMBOL, Chipseeker_annotation, meth.group1.WT, meth.group2.TG,FDR_adj_genotype) 
  }
  
  output = list(p,Location)
  names(output) = c("Venn","Positions")
  return(output)
}

# Plot the Methylation of DMP of interest, and the correlation of isoform expression with DMP 
# gene = common gene from the DMPgene doc
# DMPposition = position of the DMP of interest 
# isoformdiff = PB.ID of the isoform of interest 
Methylation_Integration_plots <- function(gene, isoformdiff){
  # Find the DMP associated to the gene from Isabel's results
  Diff_Meth = rbind(Whole_DMP$genotype[,c(c("Location","Chipseeker_SYMBOL"))], # Genotype DMP
                    Whole_DMP$interaction[,c("Location","Chipseeker_SYMBOL")]) # Interaction DMP
  
  # remove duplicated DMP (from different analysis of genotype and interaction)
  Diff_Meth = Diff_Meth[!duplicated(Diff_Meth), ]
  DMPposition = Diff_Meth %>% filter(Chipseeker_SYMBOL == gene) 
  print(DMPposition)
  
  # subset the RRBS data by the chromosome of the gene and the DMP position  
  # Use the mean of all the methylation values
  dat = subset(RRBS_completebetas, rownames(RRBS_completebetas) %in% DMPposition$Location) %>% reshape2::melt() %>%
    group_by(variable) %>% summarise_at(vars(value), list(Methylation = mean))
  
  # merge subsetted data with the phenotype data; datawrangle for plot
  df = merge(dat,RRBS_Phenotype[,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% 
    mutate(Age_months = factor(Age_months),Genotype = factor(Genotype, levels = c("WT","TG")), Sample = word(variable,c(1),sep = fixed("_")))
  
  # plot the methylation of DMP 
  age_names <- c(`2` = "2 mos",`4` = "4 mos", `6` = "6 mos",`8` = "8 mos")
  p1 = ggplot(df,aes(x = Genotype, y = Methylation, colour = Genotype)) + geom_boxplot() + 
    geom_point(size = 3, position = position_jitter(w = 0.1, h = 0)) + 
    labs(x = " ", y = "Methylation (%)") + mytheme +
    theme(legend.position = "bottom") + scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) +  
    scale_y_continuous(labels = function(x) x*100) + theme(legend.position = "none") +
    facet_grid(~Age_months,labeller = as_labeller(age_names)) + 
    theme(strip.background = element_blank())

  ggplot(df, aes(x = Genotype, y = Methylation)) + geom_boxplot() + 
    mytheme + scale_y_continuous(labels = function(x) x*100) + 
    labs(x = " ", y = "Methylation (%)") 
  
  
  # isoform expression (RNASeq expression normalised counts)
  df2 = tappasrna$input_norm %>% 
    tibble::rownames_to_column(., "isoform") %>% filter(isoform == isoformdiff) %>% reshape2::melt() #%>% 
    #mutate(sample = word(word(variable,c(1), sep = fixed("_")),c(2), sep = fixed(".")))
  df2 = merge(df2, tappasrna_phenotype, by = "variable") %>% dplyr::rename(Expression = value)
  
  # merge the isoform expression and methylation 
  merged = merge(df,df2,by.x = "Sample", by.y = "sample")
  
  # plot interaction 
  p2 = ggplot(merged,aes(x = Methylation, y = Expression)) + geom_point(size = 3, aes(colour = Genotype)) + labs(y = "Isoform Expression (normalised)", x = "Methylation (%)") + mytheme + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) + theme(legend.position = "bottom") +  scale_x_continuous(labels = function(x) x*100)
  
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
  
  output = list(p1,p2)
  return(output)
}

Methylation_Integration_stats <- function(gene){
  print(gene)
  dat = DMP_Combine %>% filter(Chipseeker_SYMBOL == gene)
  nDMP = nrow(dat)
  Loc = dat$Chipseeker_annotation[1]
  WT = signif(mean(dat$meth.group1.WT),digits = 3)
  TG = signif(mean(dat$meth.group2.TG),digits = 3)
  Diff = TG - WT 
  FDR_meth = signif(mean(dat$FDR_adj_genotype),digits = 3)
  
  Iso = tappassigtrans$WholeRNA_Transexp[tappassigtrans$WholeRNA_Transexp$associated_gene == gene,] 
  if(nrow(Iso) >= 1){
    Iso = Iso %>% arrange(desc(`p-value`)) %>% .[1,]
    type = "DTE"
    id = Iso$isoform
    trans = Iso$associated_transcript
    R = Iso$`R-squared`
    FDR_exp = signif(Iso$`p-value`,digits = 3)
    iso_diff = simple_diff_stats(id,RNAExp$Norm_transcounts,"rnaseq") 
    fc = signif(iso_diff$log2fc,digits = 3)
  }else{
    Iso = tappasDIU$rnaseq[tappasDIU$rnaseq$gene == gene,]
    type = "DTU"
    id = "NA"
    trans = "NA"
    R = "NA"
    FDR_exp = signif(Iso$`adjPValue`,digits = 3)
    fc = "NA"
  }
  
  
  output = data.frame(Gene = gene,Analysis = type,Isoform = trans,Isoform_id = id,fc = fc,FDR = FDR_exp,
                      R = R, Num = nDMP,location = Loc,Diff = Diff, FDR = FDR_meth)
  return(output)
}

as3mt_dmr_figure <- function(){
  # DMP position within AS3MT region
  as3mt_DMR = c("chr19:46730273","chr19:46730291","chr19:46730292","chr19:46730304","chr19:46730305","chr19:46730332","chr19:46730333","chr19:46730349","chr19:46730758")
  dat = subset(RRBS_completebetas, rownames(RRBS_completebetas) %in% as3mt_DMR) 
  
  # subset from RRBS and datawrangle for plot
  subsetted_dat = dat %>% rownames_to_column(var = "position") %>% mutate(location = as.numeric(as.character(word(position,c(2),sep = fixed(":"))))) %>% reshape2::melt(id = "location") %>% merge(.,RRBS_Phenotype[,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% mutate(Methylation = as.numeric(as.character(value)) * 100) 
  
  # mean of methylation across the DMR
  df = subsetted_dat %>% filter(location != "46730758") %>% 
    group_by(Genotype,Age_months,variable) %>% summarise_at(vars(Methylation), list(Methylation = mean))
  
  df2 = tappasrna$input_norm %>% tibble::rownames_to_column(., "isoform") %>% 
    filter(isoform == "PB.8363.2") %>% reshape2::melt() %>% dplyr::rename(Expression = value) %>%
    mutate(sample = word(word(variable,c(1), sep = fixed("_")),c(2), sep = fixed("."))) 
  
  # merge the isoform expression and methylation 
  merged = merge(df,df2, by.x = "variable", by.y = "sample")
  
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
  p2 = ggplot(merged,aes(x = Methylation, y = Expression)) + geom_point(size = 3, aes(colour = Genotype)) + 
    labs(y = "Isoform Expression (normalised)", x = "Methylation (%)") + mytheme + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) + 
    theme(legend.position = "bottom")
  
  
  return(list(p,p2))
}

prnp_dmr_figure <- function(){
  # DMP position within Prnp
  prnp_DMR = c("chr2:131909818","chr2:131909823","chr2:131909854","chr2:131909856","chr2:131909863","chr2:131909865","chr2:131909895",
               "chr2:131909902","chr2:131909918","chr2:131909927","chr2:131909942","chr2:131909943","chr2:131909951","chr2:131909953",
               "chr2:131909959","chr2:131909972","chr2:131909982","chr2:131909987","chr2:131909990","chr2:131910163","chr2:131910165",
               "chr2:131910181","chr2:131910202")
  
  dat = subset(RRBS_completebetas, rownames(RRBS_completebetas) %in% prnp_DMR) 
  
  # subset from RRBS and datawrangle for plot
  subsetted_dat = dat %>% rownames_to_column(var = "position") %>% mutate(location = as.numeric(as.character(word(position,c(2),sep = fixed(":"))))) %>% reshape2::melt(id = "location") %>% merge(.,RRBS_Phenotype[,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% mutate(Methylation = as.numeric(as.character(value)) * 100) 
  
  
  p1 = plot_transexp_overtime("Prnp",RNAExp$Norm_transcounts,"isoseq","RNA-Seq Expression")
  
  p2 = ggplot(subsetted_dat,aes(x = location, y = Methylation, colour = Genotype)) + 
    annotate("rect", xmin = 131909854, xmax = 131910202, ymin = 0, ymax = 100, alpha = 0.1,fill = wes_palette("Darjeeling1")[2]) +
    geom_point(size = 2) + ylim(0,100) + xlim(131909818,131910210) + stat_summary(data=subsetted_dat, aes(x=location, y=Methylation, group=Genotype), fun ="mean", geom="line", linetype = "dotted") + scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), "Genotype", labels = c("WT","TG")) + 
    mytheme + labs(y = "Methylation (%)", x = "Chromosome 2") + theme(legend.position = "bottom")
  
  p3 = plot_transexp_overtime("Prnp",IsoExp$Norm_transcounts,"isoseq","Iso-Seq Expression")
  
  return(list(p1,p2,p3))
}

