
## ---------- AD gene list -----------------

# Most recent human AD GWAS 2022
# convert human genes to mouse homologues 
ADGWAS <- read.table(paste0(dirnames$glob_metadata,"/ADGWASBellenguez.txt"), col.names = c("human"))
mouse2humanconversion <- read.csv(paste0(dirnames$utils, "mousehumangeneconversion.csv"))
ADGWAS <- merge(ADGWAS, mouse2humanconversion, by = "human")

# Full list of AD genes including Target Gene List 
ADRelatedGenes <- c(as.character(ADGWAS$mouse), TargetGene)

# List of housekeeping genes from Review
housekeepingGenes = c("Las1l","Rrp1","Gusb","Polr2b","Cyc1","Tbp","Xpnpep1","Gapdh","Actb","Rpl13a","Sdha")


## ---------- functions -----------------
meanTargetGeneExpression <- GlobalDESeq$resGeneAnno$wald$norm_counts_all %>% 
  filter(associated_gene %in% TargetGene) %>% 
  group_by(associated_gene) %>% 
  dplyr::summarise(mean = mean(normalised_counts)) 

meanOtherADGeneExpression <- GlobalDESeq$resGeneAnno$wald$norm_counts_all %>% 
  group_by(associated_gene) %>% 
  dplyr::summarise(mean = mean(normalised_counts)) %>%
  filter(min(meanTargetGeneExpression$mean) < mean, mean < max(meanTargetGeneExpression$mean))

rbind(meanTargetGeneExpression %>% mutate(Dataset = "Target"),
      meanOtherADGeneExpression %>% mutate(Dataset = "Others")) %>%
  ggplot(., aes(x = mean, fill = Dataset)) + geom_density(alpha = 0.2) +
  labs(x = paste0("Number of isoform"), y = "Density") +
  theme_classic() + scale_x_continuous(trans='log10') 

dat <- rbind(class.files$glob_iso %>% filter(associated_gene %in% TargetGene) %>% group_by(associated_gene) %>% tally() %>%
               mutate(Genes = "Target"),
             class.files$glob_iso %>% filter(!associated_gene %in% c(TargetGene)) %>% group_by(associated_gene) %>% tally() %>%
               mutate(Genes = "Other protein-coding")) 
ggplot(dat, aes(x = n, fill = Genes)) + geom_density(alpha = 0.4) +
  labs(x = paste0("Number of isoforms"), y = "Density") +
  theme_classic() + scale_x_continuous(trans='log10') +
  scale_fill_manual(values = c("green","orange","blue"))

dat %>% filter(Dataset %in% c("Target","Other protein-coding")) %>% 
  t.test(n ~ Dataset, data = .)

dat %>% filter(Dataset %in% c("AD GWAS","Other protein-coding")) %>% 
  t.test(n ~ Dataset, data = .)



expressionADGene <- function(){
  
  meanADGeneExpression <- GlobalDESeq$resGeneAnno$wald$norm_counts_all %>% 
    filter(associated_gene %in% TargetGene) %>% 
    group_by(associated_gene) %>% 
    dplyr::summarise(mean = mean(normalised_counts)) 
  
  meanOtherADGeneExpression <- GlobalDESeq$resGeneAnno$wald$norm_counts_all %>% 
    group_by(associated_gene) %>% 
    dplyr::summarise(mean = mean(normalised_counts)) %>%
    filter(min(meanADGeneExpression$mean) < mean, mean < max(meanADGeneExpression$mean))
  
  return(meanOtherADGeneExpression)
}

tallyADgene <- function(phenotype, similarLevel = NULL, all = TRUE, texpression = FALSE){
  
  geneAnnotation <- read.table(paste0(dirnames$utils,"gencode.vM22.annotation.geneannotation.txt"), header = T)
  proteinCodingGenes <- as.character(geneAnnotation[geneAnnotation$Class == "protein_coding", "GeneSymbol"])
  
  if(isTRUE(texpression)){
    message("Transcript expression")
    if(isTRUE(all)){
      dat <- group_class.files.diff[group_class.files.diff$Dataset %in% c(phenotype, "Both"),] 
    }else{
      dat <- group_class.files.diff[group_class.files.diff$Dataset == phenotype,] 
    }
    
  }else{
    # isoform diversity
    if(isTRUE(all)){
      dat <- class.files$glob_iso %>% 
        filter(isoform %in% group_class.files.diff[group_class.files.diff$Dataset %in% c(phenotype, "Both"),"isoform"]) 
    }else{
      dat <- class.files$glob_iso %>% 
        filter(isoform %in% group_class.files.diff[group_class.files.diff$Dataset == phenotype,"isoform"]) 
    }
  }
  
  # keep  only protein coding genes
  dat <- dat %>% filter(associated_gene %in% proteinCodingGenes)
    
  
  if(!is.null(similarLevel)){
    meanOtherADGeneExpression <- expressionADGene()
    dat <- dat %>% filter(associated_gene %in% c(as.character(ADRelatedGenes), as.character(meanOtherADGeneExpression$associated_gene)))
  }
  
  if(isFALSE(texpression)){
    dat <- dat %>% group_by(associated_gene) %>% tally() 
  }

  dat <- dat %>% 
    mutate(ADRelated = ifelse(associated_gene %in% c(as.character(ADRelatedGenes)),"AD","Non AD"),
           Dataset = phenotype) 
 
  return(dat)
  
}


plotIsoformDiversity <- function(common = TRUE, gexpression = FALSE, texpression = FALSE){
  
  if(isTRUE(common)){
    axistitle <- "all"
  }else{
    axistitle <- "unique"
  }
  
  # tally the number of unique isoforms in WT and TG mice
  isoDiversity_WT <- tallyADgene("WT", similarLevel = TRUE, all = common, texpression = FALSE)
  isoDiversity_TG <- tallyADgene("TG", similarLevel = TRUE, all = common, texpression = FALSE)
  
  
  if(isTRUE(texpression)){
    message("Plotting transcript expression only")
    
    isoExp_WT <- tallyADgene("WT", similarLevel = TRUE, all = common, texpression=TRUE) %>% mutate(FL = WTFL)
    isoExp_TG <- tallyADgene("TG", similarLevel = TRUE, all = common, texpression=TRUE) %>% mutate(FL = TGFL)
    
    # plot distribution
    p <- rbind(isoExp_TG[,c("associated_gene","FL", "Dataset", "ADRelated")], 
               isoExp_WT[,c("associated_gene","FL", "Dataset", "ADRelated")]) %>% 
      ggplot(., aes(x = FL, fill = Dataset)) + geom_density(alpha = 0.2) +
      labs(x = paste0("FL read count of ", axistitle, " isoforms"), y = "Density") +
      theme_classic()
    
    # stats
    resAD <- rbind(isoExp_TG[,c("associated_gene","FL", "Dataset", "ADRelated")], 
          isoExp_WT[,c("associated_gene","FL", "Dataset", "ADRelated")]) %>% 
      filter(ADRelated == "AD") %>% 
      t.test(FL ~ Dataset, data = .)

    message("AD only")
    print(resAD)

  }else if(isTRUE(gexpression)){
    
    meanOtherADGeneExpression <- expressionADGene() %>% merge(.,rbind(isoDiversity_TG, isoDiversity_WT),by = "associated_gene")
    meanOtherADGeneExpression <- meanOtherADGeneExpression %>% mutate(ADRelated = ifelse(associated_gene %in% housekeepingGenes, "HouseKeeping", ADRelated))
    
    p <- ggplot(meanOtherADGeneExpression, aes(y = mean, x = n, colour = Dataset)) + geom_point() + 
      scale_y_continuous(trans='log10') +
      facet_grid(ADRelated ~ .) +  
      labs(x = paste0("Number of ", axistitle, " isoforms"), y = "Mean gene expression") 
    
    message("AD only")
    geneExpressionIsoAD <- meanOtherADGeneExpression[meanOtherADGeneExpression$ADRelated == "AD",]
    resAD <- cor.test(geneExpressionIsoAD$mean, geneExpressionIsoAD$n)
    print(resAD)
    
    message("HK only")
    geneExpressionIsoHK <- meanOtherADGeneExpression[meanOtherADGeneExpression$ADRelated == "HouseKeeping",]
    resHK <- cor.test(geneExpressionIsoHK$mean, geneExpressionIsoHK$n)
    print(resHK)
    
    message("Other only")
    geneExpressionIsoHK <- meanOtherADGeneExpression[meanOtherADGeneExpression$ADRelated == "Non AD",]
    resHK <- cor.test(geneExpressionIsoHK$mean, geneExpressionIsoHK$n)
    print(resHK)
    
  }else{
    
    # plot distribution
    p <- rbind(isoDiversity_TG, isoDiversity_WT) %>% 
      ggplot(., aes(x = n, fill = Dataset)) + geom_density(alpha = 0.2) +
      labs(x = paste0("Number of ", axistitle, " isoforms"), y = "Density") +
      theme_classic()
    
    # t.test
    print(rbind(isoDiversity_TG, isoDiversity_WT) %>% filter(ADRelated == "AD") %>% t.test(n ~ Dataset, data = .))
  }
  
    p <- p + 
      facet_grid(ADRelated~.) + 
      scale_x_continuous(trans='log10') + 
      scale_fill_manual(values = c("red",label_colour("WT")), name = NULL) +
      scale_colour_manual(values = c("red",label_colour("WT")), name = NULL) +
      theme(strip.background = element_blank(), legend.position = c(0.9,0.9))
  
  return(p)
}


Allplots <- list(
  
  ## ---------- isoform diversity in AD-associated genes -----------------
  pTDiversityAll = plotIsoformDiversity(common = TRUE),
  pTDiversityUnique = plotIsoformDiversity(common = FALSE),
  
  ## ---------- isoform expression in AD-associated genes -----------------
  pTExpressionAll = plotIsoformDiversity(common = TRUE, texpression = TRUE),
  pTExpressionUnique = plotIsoformDiversity(common = FALSE, texpression = TRUE),

  ## ---------- gene expression of AD-associated genes -----------------
  pGExpressionAll = plotIsoformDiversity(common = TRUE, gexpression = TRUE),
  pGExpressionUnique = plotIsoformDiversity(common = FALSE, gexpression = TRUE)
  
)

plot_grid(plotlist = Allplots, ncol = 2, labels = c("A","B","C","D","E","F"), rel_heights = c(0.3,0.3,0.4))
