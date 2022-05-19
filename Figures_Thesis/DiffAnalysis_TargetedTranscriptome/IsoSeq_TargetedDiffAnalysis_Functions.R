# Szi Kay Leung
# Functions script for Thesis Chapter on Differential Analysis using Whole Transcriptome Dataset

'%!in%' <- function(x,y)!('%in%'(x,y))
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


reannotate_tamamerged_output <- function(tama_merged_transcript_file){
  # merge the Norm_transcounts from RNA-Seq Expression with tama merge annotation file, with the corresponding Pbid from targeted transcriptome 
  colnames(tama_merged_transcript_file)[[1]] = "tama_id"  # initially called "isoform" 
  dat = merge(tama_merged_transcript_file, merged_targetedid[,c("transcript_id","isoform")], by.x = "tama_id", by.y = "transcript_id", all.x = T)
  
  # note there are some transcripts that were not detected in targeted transcriptome but were in whole transcriptome and were aligned with RNA-Seq
  targeted_detected = dat[!is.na(dat$isoform),]  # these are the isoforms that were detected 
  missing_isoforms = unique(dat[is.na(dat$isoform),"tama_id"]) # these are the isoforms that were not detected
  
  # find the corresponding pbid from the whole transcriptome 
  wholeco_missing_isoforms = merged_wholeid[merged_wholeid$transcript_id %in% missing_isoforms,c("transcript_id","isoform")]
  
  # remerge the Norm_transcounts with the list of reannotated pbid from whole transcriptome
  targeted_not_detected = merge(dat %>% select(-c("isoform")),wholeco_missing_isoforms,by.x = "tama_id", by.y = "transcript_id", all.y = T) %>% 
    mutate(isoform = paste0(isoform,"_Whole"))
  
  combined = rbind(targeted_detected,targeted_not_detected)    
  combined = combined %>% mutate(Isoform = paste0(isoform,"_",structural_category))
  return(combined)
}


### TappAS Differential Analysis #################################################################
input_tappasfiles <- function(tappas_input_dir){
  # read in files generated from TAPPAS
  tappasfiles <- list(paste0(tappas_input_dir,"/Data/gene_matrix.tsv"),
                      paste0(tappas_input_dir,"/Data/gene_transcripts.tsv"),
                      paste0(tappas_input_dir,"/InputData/input_normalized_matrix.tsv"),
                      paste0(tappas_input_dir,"/Data/result_gene_trans.tsv"),
                      paste0(tappas_input_dir,"/Data/transcript_matrix.tsv"))
  
  tappasfiles <- lapply(tappasfiles, function(x) read.table(x, sep = "\t", header = T))
  names(tappasfiles) <- c("gene_matrix.tsv","gene_transcripts.tsv","input_normalized_matrix.tsv","result_gene_trans.tsv","transcript_matrix.tsv")
  
  # InputExpression Matrix from TAPPAS of which transcripts are filtered during normalisation 
  # match the filtered transcripts with the associated gene uing the isoform id from classification file
  #tappasfiles$tappAS_Transcripts_InputExpressionMatrix.tsv <- 
  #  merge(class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], 
  #        tappasfiles$tappAS_Transcripts_InputExpressionMatrix.tsv, by.x = "isoform", by.y = "Id")
  
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
  
  #classification_file = class.files
  #tappas_normalised_expmatrix = tappasrna$input_normalized_matrix.tsv
  #phenotype = tappasrna_phenotype
  # Annotate the normalised expression matrix from tappas with the associated gene and sample phenotypes for plots 
  Norm_transcounts = 
    # annotate the tappas output normalised expression matrix of transcripts by the associated gene name and transcript name
    merge(classification_file [,c("isoform","associated_gene","associated_transcript","structural_category")], 
          tappas_normalised_expmatrix, by.x = "isoform", by.y = 0) %>% reshape2::melt() %>% 
    # annotate the samples to the phenotype 
    left_join(., phenotype, by = c("variable" = "variable")) %>% 
    # change factor levels for plots
    mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("WT", "TG")),
           structural_category=recode(structural_category, `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM",`novel_in_catalog`="NIC"),
           Isoform = paste0(isoform,"_",structural_category))
  
  # Deduce gene expression from the sum of normalised transcript counts 
  GeneExp = Norm_transcounts %>% group_by(associated_gene,variable) %>% dplyr::summarise(Exp = sum(value)) %>%
    left_join(., phenotype, by = c("variable" = "variable"))
  
  output <- list(Norm_transcounts,GeneExp)
  names(output) <- c("Norm_transcounts","GeneExp")
  return(output)
}

# plot gene expression or by isoformID 
plot_mergedexp <- function(InputGene,IsoformID,GeneExp,Norm_transcounts,type){
  if (InputGene != "NA"){
    df <- GeneExp %>% filter(associated_gene == InputGene) 
    df$group <- factor(df$group, levels = c("CONTROL","CASE"))
    plot_title <- paste0(InputGene,"\n",type)
    #cat("Mean Expression across groups for", InputGene,"\n")
    
    #meanexp = df %>% group_by(time, group) %>% summarise_at(vars(Exp), funs(mean(., na.rm=TRUE))) %>% as.data.frame() %>% 
    #  mutate(groupings = paste0(time,group)) %>% .[,c("Exp","groupings")] %>% spread(., groupings, Exp) %>% mutate(log2fc = log2(`8CASE`/`2CASE`)) %>% .[,c(5,2,4,1,3)]
    
    # 3 significant factor 
    #meanexp = signif(meanexp, digits = 3)
    
    #print(meanexp)
    #meanexp_output[[InputGene]] <- meanexp
    #meanexp_output <<- meanexp_output
    
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
    labs(y = "Normalised Counts", x = "Age (months)", title = paste0(plot_title,"\n\n")) + mytheme +
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
# Input: tappassig with the all the sheets
# Plots: 
# P1: venn diagram of genes that are differentially expressed between RNA+RNA(Isabel),Iso+RNA,Iso+Iso
tappas_genesig <- function(R){
  
  ## Genes that are already filtered by significance (p <0.05 and R > 0.5)
  # Genes associated with interaction effects (results from using RNA-Seq or Iso-Seq as expression)
  WholeRNA_Interaction = c(gene_sigs_WholeRNA_lst$models$`Model 4 Interaction`$...1,gene_sigs_WholeRNA_lst$models$`Model 5 Interaction`$...1,
                           gene_sigs_WholeRNA_lst$models$`Model 6 Interaction`$...1,gene_sigs_WholeRNA_lst$models$`Model 7 Interaction`$...1)
  
  WholeIso_Interaction = c(gene_sigs_WholeIso_lst$models$`Model 4 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 5 Interaction`$...1,
                           gene_sigs_WholeIso_lst$models$`Model 6 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 7 Interaction`$...1)
  
  # Genes associated with genotype effects
  WholeRNA_Genotype = gene_sigs_WholeRNA_lst$models$`Model 1 Genotype`$...1
  WholeIso_Genotype = gene_sigs_WholeIso_lst$models$`Model 1 Genotype`$...1
  
  
  cat("Nummber of DEG from Isabel's supp, Genotype Effect:", nrow(Isabel_gene_Tg4510AgeGenotypeDEG),"\n")
  cat("Nummber of DEG from Isabel's supp, Genotype & Age Effect:", nrow(Isabel_gene_Tg4510GenotypeDEG),"\n")
  
  # first column of the input table = gene list 
  p1 <- venn.diagram(x = list(Isabel_gene_Tg4510AgeGenotypeDEG$Gene, WholeRNA_Interaction, WholeIso_Interaction), 
                     category.names = c("Reference genome, \n RNA-Seq Expression","Iso-Seq transcriptome, \n RNA-Seq Expression","Iso-Seq transcriptome, \n Iso-Seq Expression"), 
                     filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC","#CBD5E8"), main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",
                     cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n"
                     print.mode = "raw")
  
  p2 <- venn.diagram(x = list(Isabel_gene_Tg4510GenotypeDEG$Gene, WholeRNA_Genotype, WholeIso_Genotype), 
                     category.names = c("Reference genome, \n RNA-Seq Expression",
                                        "Iso-Seq transcriptome, \n RNA-Seq Expression","Iso-Seq transcriptome, \n Iso-Seq Expression"), 
                     filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC","#CBD5E8"), 
                     main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",
                     cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-30, 0, 30), cat.dist = c(0.055, 0.055, 0.065),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n"
                     print.mode = "raw")
  
  output <- list(p1,p2)
  names(output) <- c("p1","p2")
  return(output)
  
}

RNASeq_IsoSEQ_Geneexp <- function(){
  DEA_uIiso = setdiff(WholeIso_Interaction,WholeRNA_Interaction)
  DEA_uRNA = setdiff(WholeRNA_Interaction,WholeIso_Interaction)
  DEA_commonRNAIso = intersect(WholeRNA_Interaction,WholeIso_Interaction)
  
  plot_mergedexp("Cd34","NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")
  plot_mergedexp("Cd34","NA",wholetappas_rnaexp$GeneExp,wholetappas_rnaexp$Norm_transcounts,"Iso-Seq Expression")
  plot_mergedexp("Thra","NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")
  plot_mergedexp("Thra","NA",wholetappas_rnaexp$GeneExp,wholetappas_rnaexp$Norm_transcounts,"Iso-Seq Expression")
  
  tappasiso$gene_matrix_mean <- apply(tappasiso$gene_matrix.tsv,1,mean) %>% reshape2::melt() %>% `colnames<-`(c("meanIsoGeneExp"))
  wholetappas_isoexp$GeneExp <- merge(wholetappas_isoexp$GeneExp,tappasiso$gene_matrix_mean %>% rownames_to_column(var = "associated_gene"), by = "associated_gene")
  
  tappasrna$gene_matrix_mean <- apply(tappasrna$gene_matrix.tsv,1,mean) %>% reshape2::melt() %>% `colnames<-`(c("meanRNAGeneExp"))
  wholetappas_rnaexp$GeneExp <- merge(wholetappas_rnaexp$GeneExp,tappasrna$gene_matrix_mean %>% rownames_to_column(var = "associated_gene"), by = "associated_gene")
  
  meanGeneExp = merge(distinct(wholetappas_rnaexp$GeneExp[,c("associated_gene","meanRNAGeneExp")]),
                      distinct(wholetappas_isoexp$GeneExp[,c("associated_gene","meanIsoGeneExp")]), by = "associated_gene")
  for(i in 1:nrow(meanGeneExp)){
    meanGeneExp$type[[i]] = if(meanGeneExp$associated_gene[[i]] %in% DEA_uIiso){"IsoSeq_Only"
    }else if(meanGeneExp$associated_gene[[i]] %in% DEA_uRNA){"RNASeq Only"
    }else if(meanGeneExp$associated_gene[[i]] %in% DEA_commonRNAIso){"Common DEA"}else{"Not DEA"}}
  
  meanGeneExp %>% filter(type != "Not DEA") %>% 
    ggplot(., aes(x = meanRNAGeneExp, y = meanIsoGeneExp, colour = type)) + geom_point() + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10") 
  
  rbind(tappassiggene$WholeRNA_Genexp[tappassiggene$WholeRNA_Genexp$...1 %in% DEA_uRNA,c("R-squared")] %>% mutate(type = "RNASeq"),
        tappassiggene$WholeRNA_Genexp[tappassiggene$WholeRNA_Genexp$...1 %in% DEA_commonRNAIso,c("R-squared")] %>% mutate(type = "Both")) %>% 
    ggplot(., aes(`R-squared`, fill = type)) + geom_density(alpha = 0.2)
  
  rbind(tappassiggene$WholeIso_Genexp[tappassiggene$WholeIso_Genexp$...1 %in% DEA_uIiso,c("R-squared")] %>% mutate(type = "IsoSeq"),
        tappassiggene$WholeIso_Genexp[tappassiggene$WholeIso_Genexp$...1 %in% DEA_commonRNAIso,c("R-squared")] %>% mutate(type = "Both")) %>% 
    ggplot(., aes(`R-squared`, fill = type)) + geom_density(alpha = 0.2)
  
}

plot_transexp <- function(InputGene,Norm_transcounts,type, name){
  plot_title <- paste0(InputGene,"\n",name,"\n\n")
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene)
  df$grouping1 <- paste0(df$group, "_", df$time)
  df <- df[order(match(df$grouping1, c("WT_2","WT_8","TG_2","TG_8"))),]
  #df <- df %>% arrange(desc(group)) %>% arrange(desc(time)) 
  df$time <- as.factor(df$time)
  # df$group <- factor(df$group, level = c("WT","TG"))
  #df$time <- factor(paste(df$time,"months"), levels = c("2 months","8 months"))
  df$groupings <- paste0(df$isoform,"_", df$group, "_", df$time)
  
  #df$groupings <- paste0(df$isoform,"_", df$group)
  #p <-
  #  ggplot(df, aes(x = reorder(isoform,-value), y = value)) + geom_line(aes(group = groupings, colour = group), position = position_dodge(width = 0.3)) + 
  #geom_jitter(size = 3, position = position_jitterdodge()) +
  #  geom_point(aes(group = groupings, colour = group), position = position_dodge(width = 0.3)) +
  #  mytheme + labs(x = "Isoforms", y = "Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
  #  theme(strip.background = element_blank(), 
  #        plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
  #        panel.spacing = unit(2, "lines"), 
  #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #  legend_theme +
  #  facet_grid(~time,  scales = "free", space = "free") +
  #  scale_y_continuous(trans='log10') +
  #  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  
  # show only the first highest expressed isoforms
  top10_isoforms = df %>% group_by(isoform) %>% tally(value) %>% arrange(-n) %>% as.data.frame() %>% .[1:10,1]
  
  p <- df %>% filter(isoform %in% top10_isoforms) %>%
    ggplot(., aes(x = reorder(isoform,-value), y = value)) + geom_line(aes(group = groupings, colour = group), position = position_dodge(width = 0.3)) + 
    geom_point(size = 3, aes(group = groupings, shape = time, colour = group), position = position_dodge(width = 0.3)) +
    mytheme + labs(x = "Isoforms", y = "Normalised Counts",title = plot_title) +
    theme(strip.background = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    legend_theme +
    #scale_y_continuous(trans='log10') +
    guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  
  if(type == "rnaseq"){
    p = df %>% filter(isoform %in% top10_isoforms) %>%
      ggplot(., aes(x = reorder(isoform,-value), y = value)) + geom_boxplot(aes(group = groupings, colour = group),position = position_dodge(width = 0.7)) + 
      #geom_point(size = 3, aes(group = groupings, shape = time, colour = group),position = position_dodge(width = 0.7))+
      mytheme + labs(x = "Isoforms", y = "Normalised Counts",title = plot_title) +
      theme(strip.background = element_blank(), 
            plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
            panel.spacing = unit(2, "lines"), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      legend_theme +
      #scale_y_continuous(trans='log10') +
      guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  }
 
  
  #p <- df %>% mutate(grouping = paste0(group,"_",time)) %>% group_by(grouping, isoform) %>% dplyr::summarise(structural_category,time,group, mean_exp = mean(value), .groups = 'drop') %>% ggplot(., aes(x = reorder(isoform,-mean_exp), y = mean_exp)) + geom_point(aes(colour = group, shape = time),size = 3) + mytheme + 
  #  labs(x = "", y = "Mean Normalised Isoform Expression",title = plot_title) +
  #  theme(strip.background = element_blank(), 
  #        plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
  #        panel.spacing = unit(2, "lines"), 
  #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #  legend_theme +
  #  facet_grid(~structural_category,  scales = "free", space = "free") +
  #  scale_y_continuous(trans='log10') +
  #  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  

  p <- p + scale_shape_manual(name = "Age (months)", values=c(1, 16, 2, 17)) + 
      scale_colour_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT")))  
  
  return(p)
}

plot_transexp_overtime <- function(InputGene,Norm_transcounts, name){
  
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  plot_title <- paste0(InputGene,"\n",name,"\n\n")
  
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
    mytheme + labs(x = "Age (months)", y = "Normalised Counts",title = plot_title) +
    theme(strip.background = element_blank(), legend.position = "bottom",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines")) +
    guides(colour=guide_legend(ncol=3,bycol=TRUE))
  
  if(name == "ONT Expression"){
    p = p + theme(legend.position = "None")
  }
  return(p)
}

# redundant
plot_transexp_overtime_specific <- function(InputGene,Norm_transcounts, name){
  
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  plot_title <- paste0(InputGene,"\n",name,"\n\n")
  
  df = df %>% mutate(Isoform = ifelse(associated_transcript != "Novel", associated_transcript, paste0(isoform,"_",structural_category)))
  
  p <- ggplot(df, aes(x = time, y = value, colour = Isoform)) + geom_point() + 
    facet_grid(~group,scales = "free", space = "free") +
    stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
    mytheme + labs(x = "Age (months)", y = "Normalised Counts",title = plot_title) +
    theme(strip.background = element_blank(), legend.position = "bottom",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines")) +
    guides(colour=guide_legend(ncol=3,bycol=TRUE))
  return(p)
}



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
  cat("Number differentially expressed under Interaction effect:", sum(df[grepl("Interaction",df$Model),"value"]),"\n")
  cat("Number differentially expressed under Genotype effect:", sum(df[df$Model == "Model 1 Genotype","value"]),"\n")
  
  output = list(models,p)
  names(output) = c("models","p")
  return(output)
}


twovenndiagrams <- function(set1, set2, name1, name2){
  p <- venn.diagram(x = list(set1,set2), 
                    category.names = c(name1,name2),filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', 
                    fill = c("#B3E2CD", "#FDCDAC"), main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",#
                    cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-147, 147), cat.dist = c(0.055, 0.055),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n"
                    print.mode = "raw")
  return(p)
}



DIU_analysis_output_venn <- function(){
  
  # whole trancriptome: Iso-Seq scaffold + Iso-Seq expression 
  # venn diagram of the genes identified as differentially isoform usage 
  # methods for removing lowly expressed isoforms:
  # proportion: threshold at 0.2 (20%)
  # fold change: threshold at 2.5 (relative fold change)
  v1 = twovenndiagrams(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_isoseq_fc$gene,"DIU Genes: Iso-Seq Expression \n filtered by proportion","DIU Genes: Iso-Seq Expression \n filtered by fold change")
  v2 = twovenndiagrams(tappasDIU$DIU_rnaseq_prop$gene,tappasDIU$DIU_rnaseq_fc$gene,"DIU Genes: RNA-Seq Expression \n filtered by proportion","DIU Genes: RNA-Seq Expression \n filtered by fold change")
  v3 = twovenndiagrams(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_rnaseq_prop$gene,"DIU Genes: Iso-Seq Expression \n filtered by proportion","DIU Genes: RNA-Seq Expression \n filtered by proportion")
  v4 = twovenndiagrams(tappasDIU$DIU_isoseq_fc$gene,tappasDIU$DIU_rnaseq_fc$gene,"DIU Genes: Iso-Seq Expression \n filtered by fold change","DIU Genes: RNA-Seq Expression \n filtered by fold change")
  
  output <- list(v1,v2,v3,v4)
  names(output) <- c("v1","v2","v3","v4")
  return(output)
}

IF_plot <- function(gene, gene_transcripts_input, normexp_input, type){
  normexp = normexp_input
  genetrans = gene_transcripts_input
  
  # subset expression by all the isoforms associated with the gene 
  iso = subset(genetrans, geneName == gene) 
  isoexp = normexp[which(rownames(normexp) %in% iso$transcript),]
  
  # calculate mean expression of WT vs Control 
  meanisoexp = cbind(isoexp %>% dplyr::select(contains("WT")) %>% apply(.,1,mean), isoexp %>% dplyr::select(contains("TG")) %>% apply(.,1,mean)) %>% `colnames<-`(c("WT_mean", "TG_mean")) 
  
  # determine isoform fraction by mean expression of isoform over sum of mean expression of all isoforms
  IF = apply(meanisoexp, 2, function(x) x/sum(x) * 100)
  
  p = IF %>% reshape2::melt() %>% ggplot(., aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Isoform", y = "Isoform Fraction (%)") + mytheme + scale_fill_manual(values = c(label_colour("TG"),label_colour("WT")), "Genotype", labels = c("WT","TG")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = "bottom")
  
  
  ### further subset by age 
  if(type == "isoseq"){ages = c("2","8")}else{ages = c("2","4","6","8")}
  meanisoexp_ages = list()
  for(a in 1:length(ages)){
    meanisoexp_ages[[a]] <- cbind(isoexp %>% dplyr::select(contains(paste0("WT_",ages[a]))) %>% apply(.,1,mean), 
                                  isoexp %>% dplyr::select(contains(paste0("TG_",ages[a]))) %>% apply(.,1,mean)) %>% `colnames<-`(c(paste0("WT_mean_",ages[a]),paste0("TG_mean_",ages[a])))
  }
  
  meanisoexp_ages = do.call(cbind, meanisoexp_ages)
  IF_ages = apply(meanisoexp_ages, 2, function(x) x/sum(x) * 100)
  
  IF_ages_plots = IF_ages %>% reshape2::melt() %>% mutate(age = factor(word(Var2, c(3), sep = fixed("_"))), group = factor(word(Var2, c(1), sep = fixed("_")),levels = c("WT","TG"))) 
  
  p1 = ggplot(IF_ages_plots, aes(x = group, y = value, fill = age)) + geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(~ Var1) + labs(x = "Genotype", y = "Isoform Fraction (%)") + mytheme +  theme(legend.position = "bottom") 
  
  p2 = ggplot(IF_ages_plots, aes(x = group, y = value, fill = Var1)) + geom_bar(stat = "identity") + facet_grid(~ age) +
    labs(x = "Genotype", y = "Isoform Fraction (%)", fill = "Isoforms") + mytheme + theme(legend.position = "bottom", legend.direction="vertical")
  
  
  if(type == "isoseq"){
    p1 = p1 + scale_fill_manual(values = c(wes_palette("IsleofDogs2")[1], wes_palette("IsleofDogs2")[4]), name = "Age (months)")
  }else{
    p1 = p1 + scale_fill_manual(values = c(wes_palette("IsleofDogs2")[1], wes_palette("IsleofDogs2")[2], wes_palette("IsleofDogs2")[3], 
                                           wes_palette("IsleofDogs2")[4]), name = "Age (months)")
  }
  
  output <- list(p,p1,p2)
  names(output) <- c("IFgen","IFage1","IFage2")
  return(output)
}

# plot expression of DIU genes
# DIU genes identified from using Iso-Seq expression as abundance, and commonly identified from both filtering methods (proportion and fold change)
# Gene expression = sum of normalised expression or FL read counts from associated isoforms that were prefiltered 
# prefiltered isoforms = filtered lowly-expressed isoforms as part of normalisation (refer to prefiltering prior to TMM)
# mean gene expression/median gene expression = mean across all the samples 
# 1. Determined normalised gene expression/read counts across each sample by summing normalised isoform expression/read counts 
# 2. Determined mean or median of normalised gene expression across all samples (n = 12)
DIU_genes_exp <- function(){
  # common genes identified by both proporition and fold change
  # commonDIU_isoiso: whole transcriptome: Iso-Seq scaffold + RNA_Seq expresion 
  commonDIUisoiso = intersect(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_isoseq_fc$gene)
  
  #commonDIUisorna_prop = intersect(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_rnaseq_prop$gene)
  #commonDIUisorna_fc = intersect(tappasDIU$DIU_isoseq_fc$gene,tappasDIU$DIU_rnaseq_fc$gene)
  #commonDIUisorna = intersect(commonDIUisorna_fc,commonDIUisorna_prop)
  
  ###1. Using classification files (whole transcriptome), obtain the total gene FL counts associated with the isoforms that were retained from pre-filtering 
  # normalised_matrix.tsv = isoforms that were retained after prefiltering
  # subset the classification files on retained isoforms, and obtain total gene FL counts (across all samples)
  # genes_iso = number of prefiltered isoforms associated per gene
  retained_transcripts = row.names(tappasiso$input_normalized_matrix.tsv)
  FL_genes = class.files %>% filter(isoform %in% retained_transcripts) %>% group_by(associated_gene) %>% tally(FL) 
  genes_iso = class.files %>% filter(isoform %in% retained_transcripts) %>% group_by(associated_gene) %>% tally()
  
  ###2. Obtain the median and mean of gene FL counts across the samples 
  # sum_dat = sum of the FL counts per gene per sample (FL counts of retained isoforms)
  # med_sumdat = median of the gene across all samples 
  # mean_sumdat mean of the gene across all samples
  dat = class.files %>% filter(isoform %in% retained_transcripts) %>% select(associated_gene,isoform, contains("FL.")) 
  rownames(dat) <- paste(dat$isoform,"_",dat$associated_gene)
  sum_dat = dat %>% select(-isoform) %>% group_by(associated_gene) %>% summarise(across(everything(), ~ sum(., is.na(.), 0))) %>% column_to_rownames(., var = "associated_gene")
  med_sumdat = data.frame(apply(sum_dat,1,median)) %>% rownames_to_column(var = "gene") %>% `colnames<-`(c("gene", "median_FL_expression"))
  mean_sumdat = data.frame(apply(sum_dat,1,mean)) %>% rownames_to_column(var = "gene") %>% `colnames<-`(c("gene", "mean_FL_expression"))
  
  # Merge all the counts for DIU genes 
  ###3. Use only the commonly identified genes from lowly-expressed filtering methods: proportion and fold change 
  merged_counts = merge(merge(merge(tappasDIU$DIU_isoseq_prop, FL_genes, by.x = "gene", by.y = "associated_gene"),mean_sumdat, by = "gene"),med_sumdat, by = "gene") %>% filter(gene %in% commonDIUisoiso)
  
  ## check normalised gene expression make sense using example of 1110032F04Rik
  ## iso = isoforms associated with 1110032F04Rik, note can not just use pb.id as wrong assignment 
  ## iso_exp = obtain the normalised expression of those isoforms 
  ## sum the iso_exp per column across the samples --> gene expression --> determine median and mean across samples
  ## median and mean == merged_counts$mean_expression, and merged_counts$median_expression
  #iso = tappasiso$result_gene_trans.tsv[tappasiso$result_gene_trans.tsv$gene == "1110032F04Rik","isoform"]
  #iso_exp = tappasiso$input_normalized_matrix.tsv %>% rownames_to_column(var = "isoform") %>% filter(isoform %in% iso) 
  #median(apply(iso_exp[,-1],2,sum))
  #mean(apply(iso_exp[,-1],2,sum))
  #head(merged_counts)
  
  # plots of gene expression of commonly identified DTU genes, with the normalised gene expression, and FL reads
  p1 = density_plot(merged_counts,"sum_expression","n","Normalised Gene Counts","Gene FL Reads","") + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")
  p2 = density_plot(merged_counts,"median_expression","median_FL_expression","Median Normalised Gene Counts","Median Gene FL Reads","") + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")
  p3 = density_plot(merged_counts,"mean_expression","mean_FL_expression","Mean Normalised Gene Counts","Mean Gene FL Reads","") + scale_y_continuous(trans = "log10", labels = scales::number_format(accuracy = 1)) + scale_x_continuous(trans = "log10",labels = scales::number_format(accuracy = 1))
  
  output = list(p1,p2,p3)
  names(output) = c("p1","p2","p3")
  return(output)
}

# Using Fold Change from RNASeq differential isoform usage
# Podium Change = TRUE or FALSE for whether differential isoform usage
DIU_RNASEQ_results <- function(){
  cat("Number of Genes with DIU:", nrow(tappasDIU$DIU_rnaseq_fc %>% filter(adjPValue < 0.05)),"\n")
  cat("Number of Genes with DIU with no major switching:", nrow(tappasDIU$DIU_rnaseq_fc %>% filter(adjPValue < 0.05 & podiumChange == "FALSE")),"\n")
  cat("Number of Genes with DIU with no switching:", nrow(tappasDIU$DIU_rnaseq_fc %>% filter(adjPValue < 0.05 & podiumChange == "TRUE")),"\n")
  
  # significant DEG with interaction effects using Iso-Seq as expression
  genesigs_interaction_wholeiso = c(gene_sigs_WholeIso_lst$models$`Model 4 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 5 Interaction`$...1,
                                    gene_sigs_WholeIso_lst$models$`Model 6 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 7 Interaction`$...1)
  
  # significant DEG with interaction effects using RNA-Seq as expression
  genesigs_interaction_wholerna = c(gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1)
  
  # Cases of Differentially Expressed Genes with Differential Isoform Usage
  DEG_DIU = list(
    #case1 = DIU genes with major switching isoform, and DEG
    #case2 = DIU genes but no major switching, and DEG
    #case3 = DIU genes with no major switching, but not DEG
    #case4 = DIU genes with major switching, but not DEG
    tappasDIU$DIU_rnaseq_fc %>% filter(adjPValue < 0.05 & podiumChange == "TRUE" & gene %in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue)),
    tappasDIU$DIU_rnaseq_fc %>% filter(adjPValue < 0.05 & podiumChange == "FALSE" & gene %in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue)),
    tappasDIU$DIU_rnaseq_fc %>% filter(adjPValue < 0.05 & podiumChange == "FALSE" & gene %!in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue)),
    tappasDIU$DIU_rnaseq_fc %>% filter(adjPValue < 0.05 & podiumChange == "TRUE" & gene %!in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue))
  )
  names(DEG_DIU) = c("DIU_DEG_maj","DIU_DEG_nomaj","DIU_notDEG_nomaj","DIU_notDEG_maj")
  
  cat("Number of Genes with DIU with major switching, but genes are differentially expressed:", nrow(DEG_DIU$DIU_DEG_maj),"\n")
  cat("Number of Genes with DIU with no major switching, but genes are differentially expressed:", nrow(DEG_DIU$DIU_DEG_nomaj),"\n")
  cat("Number of Genes with DIU with no major switching, but genes are not differentially expressed:", nrow(DEG_DIU$DIU_notDEG_nomaj),"\n")
  cat("Number of Genes with DIU with major switching, but genes are not differentially expressed:", nrow(DEG_DIU$DIU_notDEG_maj),"\n")
  
  return(DEG_DIU)
}

diff_across_rnaseq <- function(gene){
  p <- list(plot_mergedexp(gene,"NA",wholetappas_rnaexp$GeneExp,wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Gene Expression"),
            plot_transexp_overtime(gene,wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Isoform Expression") + theme(legend.position = "bottom", legend.direction = "vertical"),
            IF_plot(gene,tappasrna$gene_transcripts.tsv, tappasrna$input_normalized_matrix.tsv, "rnaseq")[[1]],
            IF_plot(gene,tappasrna$gene_transcripts.tsv, tappasrna$input_normalized_matrix.tsv, "rnaseq")[[3]])
  return(p)
}


IR_ORF <- function(){
  dat = lapply(group.class.files, function(x) x %>% filter(subcategory == "intron_retention") %>% group_by(associated_gene) %>% tally())
  IR = bind_rows(dat, .id = "Sample") %>% spread(., Sample, n) %>% mutate(TG_WT_8mos = `TG_8mos` - `WT_8mos`,TG_WT_2mos = `TG_2mos` - `WT_2mos`)
  
  total_genes = sapply(group.class.files, function(x) length(unique(x$associated_gene))) %>% reshape2::melt(value.name = "total_genes") %>% rownames_to_column(var = "Sample")
  total_iso = sapply(group.class.files, function(x) nrow(x)) %>% reshape2::melt(value.name = "total_isoforms") %>% rownames_to_column(var = "Sample")
  IR_genes = bind_rows(dat, .id = "Sample") %>% group_by(Sample) %>% tally() %>% full_join(.,total_genes, by = "Sample") %>% mutate(perc = n/total_genes * 100)
  
  
  plot_Genes <- function(dataset,grp, age_name, y_name){
    x.var <- rlang::sym(quo_name(enquo(grp)))
    cols <- enquo(grp) 
    
    
    # difference in number of isoforms tally 
    dat2 = dataset %>% group_by_at(vars(!!cols)) %>% tally()
    dat2$sign_freq = dat2$n * ifelse(dat2[[grp]] > 0 , 1, -1)
    dat2[[grp]]  = reorder(dat2[[grp]] , dat2$sign_freq, sum)
    
    if(grp == "TG_WT_8mos"){
      dat2_agg = aggregate(sign_freq ~ TG_WT_8mos, FUN = sum, data = dat2)
    } else{
      dat2_agg = aggregate(sign_freq ~ TG_WT_2mos, FUN = sum, data = dat2) 
    }
    dat2_agg[[grp]] = as.numeric(as.character(dat2_agg[[grp]]))
    dat2_agg$col = ifelse(dat2_agg[[grp]] == "0", "neutral", ifelse(dat2_agg$sign_freq < 0, "negative","positive"))
    
    print(dat2_agg)
    p <- ggplot(dat2_agg, aes(x = !! x.var , y = abs(sign_freq), fill = factor(col))) +
      geom_bar(position = 'identity', stat = "identity") +
      guides(fill = FALSE) + mytheme +
      labs(y = y_name, x = paste0("Difference in number of IR isoforms per gene between WT and TG \n at ", age_name,"(TG - WT)"))
    
    return(p)
  }
  
  # number of ORFs
  ORF = lapply(group.class.files, function(x) x %>% filter(predicted_NMD == "TRUE") %>% group_by(associated_gene) %>% tally())
  ORF_genes = bind_rows(ORF, .id = "Sample") %>% group_by(Sample) %>% tally() %>% full_join(.,total_genes, by = "Sample") %>% mutate(perc = n/total_genes * 100)
  ORF_transcripts = bind_rows(ORF, .id = "Sample") %>% group_by(Sample) %>% tally(n) %>% full_join(.,total_iso, by = "Sample") %>% mutate(perc = n/total_isoforms * 100)
  
  ORF_counts = bind_rows(ORF, .id = "Sample") %>% spread(., Sample, n) %>% mutate(TG_WT_8mos = `TG_8mos` - `WT_8mos`,TG_WT_2mos = `TG_2mos` - `WT_2mos`)
  
  p1 = plot_Genes(IR,"TG_WT_8mos","8 months","Number of Genes with intron-retained isoforms")
  p2 = plot_Genes(IR,"TG_WT_2mos","2 months","Number of Genes with intron-retained isoforms")
  p3 = plot_Genes(ORF_counts,"TG_WT_8mos","8 months","Number of Genes with isoforms predicted for ORF")
  p4 = plot_Genes(ORF_counts,"TG_WT_2mos","2 months","Number of Genes with isoforms predicted for ORF")
  
  return(list(p1,p2,p3,p4))
}


### Summary Info #################################################################
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
    
    colnames(output)[count] <- d
    count = count + 1
  }
  
  row.names(output) <- c("Total Number of Genes", "Annotated Genes","Novel Genes","Total Number of Isoforms","FSM",
                         "ISM","NIC","NNC","Genic\nGenomic","Antisense","Fusion","Intergenic","Genic\nIntron","Isoform Length (bp)","Number of Exons","Number of Ifsoforms with 50bp CAGE")
  
  
  output <- output[,c("WT_sq","TG_sq","WT_2m","WT_8m","TG_2m","TG_8m")]
  colnames(output) <- c("WT","TG","WT 2months","WT 8months", "TG 2months","TG 8months")
  return(output)
}

DEI_genes_unique <- function(x){
  DEI_IsoSeq = unique(c(trans_sigs_WholeIso_lst$models$`Model 4 Interaction`$isoform,trans_sigs_WholeIso_lst$models$`Model 5 Interaction`$isoform,
                        trans_sigs_WholeIso_lst$models$`Model 6 Interaction`$isoform,trans_sigs_WholeIso_lst$models$`Model 7 Interaction`$isoform))
  
  DEI_IsoSeq_only = setdiff(DEI_IsoSeq,unique(trans_sig_WholeRNA_lst$models$`Model 4 - 7 Interaction`$isoform))
  DEI_IsoSeq_RNAseq = intersect(DEI_IsoSeq,unique(trans_sig_WholeRNA_lst$models$`Model 4 - 7 Interaction`$isoform))
  cat("Number of Differentially expressed isoforms (Iso-Seq) associated with genotype and interaction effects:", length(DEI_IsoSeq),"\n")
  cat("Number of Differentially expressed isoforms (Iso-Seq) not supported by RNA-Seq expression input:",length(DEI_IsoSeq_only),"\n")
  
  DEI_Isoseq_only_FL = class.files[class.files$isoform %in% DEI_IsoSeq_only,] %>% select(starts_with("FL.")) %>% apply(.,1,mean) %>% reshape2::melt() %>% merge(.,class.files[,c("isoform","associated_gene")],by.x = 0, by.y = "isoform", all.x = T) 
  
  DEI_IsoSeq_RNAseq_FL = class.files[class.files$isoform %in% DEI_IsoSeq_RNAseq,] %>% select(starts_with("FL.")) %>% apply(.,1,mean) %>% reshape2::melt() %>% merge(.,class.files[,c("isoform","associated_gene")],by.x = 0, by.y = "isoform", all.x = T) 
  cat("Number of DIE (IsoSeq) with median FL reads < 24):", nrow(DEI_Isoseq_only_FL %>% filter(value < 24)))
  
  p1 = ggplot(DEI_Isoseq_only_FL, aes(x = value)) + geom_density(colour = "red") + mytheme + labs(x = "Mean Number of FL Reads", y= "Density")
  p2 = plot_transexp_overtime("Slc1a3",wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression") + theme(legend.position = "bottom",legend.direction="vertical")
  p3 = plot_transexp_overtime("Slc1a3",wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Expression") + theme(legend.position = "bottom",legend.direction="vertical")
  p4 = plot_transexp_overtime("Gja1",wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression") + theme(legend.position = "bottom",legend.direction="vertical")
  p5 = plot_transexp_overtime("Gja1",wholetappas_rnaexp$Norm_transcounts,"RNA-Seq Expression") + theme(legend.position = "bottom",legend.direction="vertical")
  
  output = list(p1,p2,p3,p4,p5)
  names(output) = list("p1_density","p2_slc1a3_isoseq","p3_slc1a3_rnaseq","p4_Gja1_isoseq","p5_Gja1_rnaseq")
  return(output)
}
