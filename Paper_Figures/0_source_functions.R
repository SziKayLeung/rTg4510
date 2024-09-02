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
suppressMessages(library(rtracklayer))


## ----------Functions-----------------

# load all the functions
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "aesthetics_basics_plots/draw_venn.R"))
source(paste0(LOGEN, "aesthetics_basics_plots/draw_density.R"))
source(paste0(LOGEN, "transcriptome_stats/plot_basic_stats.R"))
source(paste0(LOGEN,"longread_QC/plot_cupcake_collapse_sensitivty.R"))
source(paste0(LOGEN, "compare_datasets/base_comparison.R"))
source(paste0(LOGEN, "compare_datasets/whole_vs_targeted.R"))
source(paste0(LOGEN, "differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN, "differential_analysis/plot_usage.R"))
source(paste0(LOGEN, "merge_characterise_dataset/run_ggtranscript.R"))
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)

## ----------Theme-----------------

label_colour <- function(genotype){
  if(genotype %in% c("WT","Control","CONTROL")){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "WT_2mos"){colour = alpha(wes_palette("Royal1")[2],0.5)}else{
      if(genotype %in% c("TG","Case","CASE")){colour =wes_palette("IsleofDogs1")[4]}else{
        if(genotype == "TG_2mos"){colour = alpha(wes_palette("Royal1")[1],0.5)}else{
          if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
            if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
              if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
                if(genotype %in% c("AD")){colour = wes_palette("Royal1")[2]}else{
              }}}}}}}}
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
                     point.padding = 0.1,
                     force         = 100,
                     #nudge_y = 0.2,
                     nudge_x = 0.1,
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
    dplyr::select(isoform, log2normalised, sample) %>%
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


# recapitulate RNA-Seq and Iso-Seq at the gene level (whole transcriptome dataset)
recapitulate_gene_level <- function(){
  
  # effect size up and down in global DESeq2 output at gene level 
  GlobalDESeq$RresGeneAnno$wald$anno_res <- GlobalDESeq$RresGeneAnno$wald$anno_res %>% mutate(direction = ifelse(log2FoldChange < 0, "down","up"))
  GlobalDESeq$resGeneAnno$wald$res_wald_rResGeneSig <- GlobalDESeq$resGeneAnno$wald$res_Wald %>% 
    filter(isoform %in% GlobalDESeq$RresGeneAnno$lrt$anno_res$isoform) %>% 
    mutate(direction = ifelse(log2FoldChange < 0, "down","up"))
  
  # merge gene results from rnaseq and isoseq 
  # differentiate log2FC and direction betwen two datasets 
  
  mergedGeneRep <- merge(GlobalDESeq$RresGeneAnno$wald$anno_res %>% 
                           select(isoform, associated_gene, log2FoldChange, direction) %>% dplyr::rename("RresLog2FC" = "log2FoldChange", "RresDirection" = "direction"),
                         GlobalDESeq$resGeneAnno$wald$res_wald_rResGeneSig %>% 
                           select(isoform, log2FoldChange, direction) %>% dplyr::rename("resLog2FC" = "log2FoldChange", "resDirection" = "direction"), 
                         by = "isoform"
  )
  
  # create new column determining if direction is consistent 
  mergedGeneRep <- mergedGeneRep %>% mutate(effectSize = ifelse(RresDirection == resDirection, TRUE, FALSE))
  p <- ggplot(mergedGeneRep, aes(x = RresLog2FC, y = resLog2FC)) + geom_point() + mytheme + 
    labs(x = "RNA-Seq dataset: gene Log2FC", y = "Iso-Seq dataset: gene Log2FC")
  
  # binomial p-value 
  message("Binomial test of number of genes with consisitent effect size")
  consistentNum = nrow(subset(mergedGeneRep, effectSize == TRUE))
  totalNum = nrow(GlobalDESeq$RresGeneAnno$wald$anno_res)
  consistentNum/totalNum
  res <- binom.test(consistentNum, totalNum, alternative = c("two.sided"))
  print(res)
  print(res$p.value)
  
  # correlation 
  message("correlation test of gene Log2FC between RNA-Seq and Iso-Seq")
  res <- cor.test(mergedGeneRep$RresLog2FC,mergedGeneRep$resLog2FC, method = "pearson")
  print(res)
  print(res$p.value)
  
  return(p)
}

twovenndiagrams <- function(set1, set2, name1, name2){
  p <- venn.diagram(x = list(set1,set2), 
                    label_alpha = 0, category.names = c(name1,name2),filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', 
                    fill = c("#B3E2CD", "#FDCDAC"), main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",
                    cat.cex = 1,  cat.default.pos = "outer",  cat.col = c("#60756c", "#ba7443"),
                    cat.pos = c(-145, 200), cat.dist = c(-0.15,-0.03),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n"
                    print.mode = "raw")
  return(p)
}

tabulateIF <- function(classf, countcol){
  
  Counts <- classf %>% select(isoform,contains(countcol))
  rownames(Counts) <- Counts$isoform
  Counts <- Counts %>% select(-isoform)
  
  
  # Calculate the mean of normalised expression across all the samples per isoform
  meandf <- data.frame(meanvalues = apply(Counts,1,mean)) %>%
    rownames_to_column("isoform") %>% 
    # annotate isoforms with associated_gene and structural category
    left_join(., classf[,c("isoform","associated_gene","structural_category")], by = "isoform")  
  
  # Group meandf by associated_gene and calculate the sum of mean values for each group
  grouped <- aggregate(meandf$meanvalues, by=list(associated_gene=meandf$associated_gene), FUN=sum)
  
  # Calculate the proportion by merging back, and divide the meanvalues by the grouped values (x)
  merged <- meandf %>% 
    left_join(grouped, by = "associated_gene") %>%
    mutate(perc = meanvalues / x * 100) 
  return(merged)
}


plot_boxplot_SCN <- function(normCounts, iso, ageDiv = TRUE, genotypeDiv = TRUE){
  
  dat <- normCounts %>% dplyr::filter(isoform %in% iso) %>% dplyr::mutate(genotype = factor(genotype, levels = c("WT","TG")), cell = factor(ifelse(cell == "DN", "NeuN-", "NeuN+"), levels = c("NeuN+","NeuN-"))) 
  
  dat <<- dat
  
  p <- ggplot(dat, aes(x = cell, y  = TPM, colour = genotype)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(group = genotype), size = 3, position = position_jitterdodge()) + 
    labs(x = "Nuclei population", y = "TPM") + mytheme +
    scale_colour_manual(values = c(label_colour("WT"),label_colour("TG")), labels = c("WT","TG"), name = NULL) 
  
  print(summary(lm(TPM ~ cell + genotype + cell * genotype, data = dat)))
  
  plog2FC <- mean(log2(dat[dat$genotype  == "TG", "TPM"] + 1e-10)) -
    mean(log2(dat[dat$genotype  == "WT", "TPM"] + 1e-10))
  
  clog2FC <- mean(log2(dat[dat$cell  == "NeuN-", "TPM"] + 1e-10)) -
    mean(log2(dat[dat$cell == "NeuN+", "TPM"] + 1e-10))
  
  message("log2FC for TG vs WT: ", plog2FC)
  message("log2FC for NeuN- vs NeuN+: ", clog2FC)
  
  print(summary(lm(TPM ~ cell + genotype, data = dat)))
  
  if(length(iso) >= 2){
    
    p <- p + facet_grid(~isoform) +
      theme(strip.background = element_blank(), panel.spacing = unit(2, "lines"))
    
  }else{
    p <- p
    
    if(isTRUE(ageDiv)){
      p <- p +
        geom_point(aes(group = genotype), size = 3, position = position_jitterdodge()) + 
        facet_grid(~time, labeller = labeller(time = as_labeller(c("2m" = "2 months", "8m" = "8 months")))) +
        theme(strip.background = element_blank(), panel.spacing = unit(2, "lines"))
    }
    
    if(!isTRUE(ageDiv) & !isTRUE(genotypeDiv)){
      p <- ggplot(dat, aes(x = cell, y  = TPM)) + 
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 3) + 
        labs(x = "Nuclei population", y = "TPM") + mytheme 
      
    }
  }
  

  return(p)
  
}


BDR_plot <- function(norm_counts, iso = NULL, gene = NULL, sampleExclude = NULL, IF = NULL, Braak = FALSE){
  if(!is.null(iso) & is.null(IF) & is.null(gene)){
    print("Isoform")
    print(iso)
    dat <- norm_counts %>% filter(isoform == iso)
  }else if(!is.null(gene)){
    dat <- norm_counts %>% filter(grepl(gene,isoform)) %>%
      group_by(associated_gene, sample) %>% tally(normalised_counts, name = "normalised_counts")
  }else if(!is.null(IF)){
    print("IF")
    dat <- norm_counts %>% 
      map_if(is.numeric, ~./sum(.) * 100) %>%
      as_data_frame() %>% 
      filter(isoform == iso ) %>% 
      reshape2::melt(variable.name = "sample", value.name = "normalised_counts") 
  }else{
    print("error")
  }
  
  
  dat <- dat %>%
    mutate(sample = str_remove(sample, "B2.")) %>% 
    left_join(., phenotype, by = "sample") %>%
    mutate(BraakTangle_numeric = as.factor(BraakTangle_numeric))%>%
    filter(BraakTangle_numeric %in% c(0,1,2,5,6))  %>%
    mutate(phenotype = ifelse(BraakTangle_numeric %in% c(0,1,2),"Control","AD")) %>%
    mutate(phenotype = factor(phenotype, levels = c("Control","AD")))
  
  if(isFALSE(Braak)){
    p <- ggplot(dat, aes(x = phenotype, y = normalised_counts, fill = phenotype)) + geom_boxplot() + mytheme
  }else{
    p <- ggplot(dat, aes(x = BraakTangle_numeric, y = normalised_counts, fill = phenotype)) + geom_boxplot() + mytheme
  }
  
  if(!is.null(IF) & isFALSE(Braak)){
    p <- p + labs(x = "Phenotype", y = "Isoform fraction (%)")
  }else if(is.null(IF) & isTRUE(Braak)){
    p <- p + labs(x = "Braak Stage", y = "Normalized counts")
  }else if(!is.null(IF) & !isFALSE(Braak)){
    p <- p + labs(x = "Braak Stage", y = "Isoform fraction (%)")
  }else{
    p <- p + labs(x = "Phenotype", y = "Normalized counts")
  }
  
  # stats
  if(isFALSE(Braak)){
    sdat <- dat %>% filter(!sample %in% sampleExclude) %>% as.data.frame()
    sdat <<- sdat
    res <- t.test(normalised_counts ~ phenotype, data = sdat)
    log2FC <- mean(log2(sdat[sdat$phenotype == "AD", "normalised_counts"] + 1e-10)) -
      mean(log2(sdat[sdat$phenotype == "Control", "normalised_counts"] + 1e-10)) 
    
    print(res)
    message("Log2FC:", log2FC)
  }
  
  return(p)
}

# pclassfile with columns: <corrected_acc> from refined proteogenomics pipeline
plot_protein_general <- function(tclassfile, pclassfile, correctedCollapsedID = NULL){
  
  # count the number of transcripts in the transcript classification file
  RNATranscript <- tclassfile %>% group_by(associated_gene) %>% dplyr::summarize(RNATranscript = n())
  
  # count the number of transcripts in the protein classification file
  RNAIsoform <- pclassfile %>% dplyr::select(corrected_acc, associated_gene) %>% dplyr::filter(!is.na(corrected_acc)) %>% distinct() %>% 
    group_by(associated_gene) %>% dplyr::summarize(RNAIsoform = n())
  NumTranscriptIsoform <- merge(RNATranscript, RNAIsoform, by = "associated_gene") 
  
  if(is.null(correctedCollapsedID)){
    UniqueORF <- tclassfile %>% filter(!is.na(corrected_acc)) %>% 
      filter(!isoform %in% pclassfile$corrected_acc) %>% 
      group_by(structural_category, subcategory) %>% tally() 
  }else{
    print("Using corrected CollapsedID")
    UniqueORF <- tclassfile %>%  
      filter(isoform %in% correctedCollapsedID) %>% 
      group_by(structural_category, subcategory) %>% tally() 
  }
  
  p1 <- ggplot(NumTranscriptIsoform, aes(x = RNAIsoform, y = RNATranscript, label = associated_gene)) + 
    geom_point() + geom_text_repel() +
    geom_abline(intercept = 0, linetype = "dashed") + mytheme + 
    labs(x = "Number of protein isoforms", y = "Number of RNA transcripts")
  
  cortest <- cor.test(NumTranscriptIsoform$RNAIsoform, NumTranscriptIsoform$RNATranscript)
  print(cortest)
  
  p2 <- ggplot(UniqueORF, aes(y = n, x = subcategory)) + geom_bar(stat = "identity") + 
    facet_grid(rows = vars(structural_category), scales = "free", space = "free") + 
    coord_flip() + 
    labs(y = "Number of redundant RNA transcripts", x = "Subcategory") + 
    mytheme + theme(strip.background = element_blank())
  
  return(list(p1,p2))
}


generalggTranPlots <- function(isolist, inputgtf, classfiles, gene, cpat = NULL, species = NULL, squish = FALSE, pfam = NULL){
  IsoDf <- data.frame(
    Isoform = unlist(IsoDf <- isolist),
    Category = rep(names(IsoDf), lengths(IsoDf))
  )
  IsoDf$colour <- c(rep(NA,length(IsoDf$Category[IsoDf$Category != "DTE"])))
  IsoDf <<- IsoDf
  p <- ggTranPlots(inputgtf=inputgtf,classfiles=classfiles,
                   isoList = c(as.character(IsoDf$Isoform)),
                   selfDf = IsoDf, gene = gene, inputCpat = cpat, cpatSpecies = species, squish = squish, inputPfam = pfam) +
    labs(y = "")
  return(p)
}


# visualisation of Trem2 transcripts with same ORF LR.Trem2.54

visualise_ORFs <- function(refgtf,tgtf, pgtf, tclassfiles, pclassfiles, gene, transcript = NULL, cpat = NULL, species = NULL){
  
  refIDs <- unique(refgtf[refgtf$gene_name == gene & !is.na(refgtf$transcript_id), "transcript_id"])
  if(!is.null(transcript)){
    
    pID <- pclassfiles %>% dplyr::filter(corrected_acc %in% transcript) %>% .[["isoform"]]
    datIDs <- list(
      Reference = refIDs,
      `RNA Transcript` = pclassfiles %>% dplyr::filter(corrected_acc == transcript) %>% .[["isoform"]],
      `Protein Isoform` = unique(pgtf %>% dplyr::filter(transcript %in% pID) %>% .[["gene_id"]])
    )
      
      p <- generalggTranPlots(datIDs, tgtf, tclassfiles, gene, cpat, species) 
  }else{
    
    tID <- unique(pclassfiles[pclassfiles$associated_gene == gene,"corrected_acc"])
    
    tORFID <- unique(pgtf %>% dplyr::filter(transcript %in% tID) %>% .[["gene_id"]])
    transcriptIDs <- list(Reference = refIDs, `RNA Transcript` = tID) 
    proteinIDs <- list(Reference = refIDs,`Protein Isoform` = tORFID) 
    
    pT <- generalggTranPlots(transcriptIDs, tgtf, tclassfiles, gene)
    pP <- generalggTranPlots(proteinIDs, pgtf, pclassfiles, gene, cpat, species)
    
    p <- plot_grid(pT, pP, labels = c("i","ii"))
    
  }
  
  return(p)
}

visualise_ORFs_NMD <- function(transcripts, gene){
  transcriptList <- list(
    Reference =  RefIsoforms[[gene]][1:2],
    `NMD` = as.character(transcripts)
  )
  proteinList <- list(
    `NMD Protein` = unique(gtf$ptarg_merged %>% filter(transcript %in% as.character(transcripts)) %>% .[["gene_id"]])
  )
  
  bothList <- c(transcriptList, proteinList)
  print(bothList)
  p <- generalggTranPlots(bothList, gtf$targ_merged, class.files$targ_filtered, gene) 
  return(p)
}

plot_expression_summed <- function(TList, TName = NULL, AgeDiv = FALSE){
  
  if(length(TList) > 0){
    gene <- class.files$targ_filtered[class.files$targ_filtered$isoform %in% TList, "associated_gene"]
    message("************************ Gene:", (unique(gene)))
    
    dat <- subset(Exp$targ_ont$normAll, row.names(Exp$targ_ont$normAll) %in% TList) %>% dplyr::select(-associated_gene) %>%
      apply(.,2,sum) %>% 
      reshape2::melt(value.name = "sumReads") %>% 
      tibble::rownames_to_column(., var = "sample") %>% 
      mutate(sample = word(sample, c(2), sep = fixed("_"))) %>% 
      left_join(., phenotype$targ_ont, by = "sample") %>%
      mutate(group = factor(ifelse(group == "CASE","TG","WT"), levels = c("WT","TG")))
    
    dat <<- dat
    
    
    if(isFALSE(AgeDiv)){
      p <- ggplot(dat, aes(x = group, y = sumReads)) + geom_boxplot(outlier.shape = NA) + 
        geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), aes(colour = factor(time)), size = 3) + 
        mytheme + labs(x = "", y = "Normalized counts", subtitle = TName) + 
        scale_colour_manual(values = c("grey","azure4","black","red"), name = "Age (months)") +
        theme(legend.position = "top") 
    }else{
      p <- ggplot(dat, aes(x = group, y = sumReads, colour = as.factor(time))) + geom_boxplot() +
        geom_point(position=position_jitterdodge(), aes(colour = factor(time)), size = 3) +
        mytheme + labs(x = "", y = "Normalized counts", subtitle = TName) + 
        scale_colour_manual(values = c("grey","azure4","black","red"), name = "Age (months)") +
        theme(legend.position = "top") 
    }
    
    message("linear regression:")
    res <- lm(sumReads ~ group + time, data = dat)
    print(summary(res))
    
    ##transform our data into log2 base.
    # add 1 to deal with 0
    dat$sumReads <- dat$sumReads + 1
    dat <- dat %>% mutate(log2Reads = log2(sumReads))
    control <- mean(dat[dat$group == "WT","log2Reads"])
    TG <- mean(dat[dat$group == "TG","log2Reads"])
    foldchange <- TG - control
    message("log2FC of TG relative to control")
    print(foldchange)
    
    return(p)
    
  }else{
    message("********************* No Transcripts")
  }
  
  
}

