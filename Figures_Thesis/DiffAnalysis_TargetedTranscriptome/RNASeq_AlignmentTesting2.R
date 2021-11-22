library(VennDiagram)
root_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/RNASeq_SQANTI3/TESTING/"
ERCC_sq= read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/ERCC/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt", header = T, as.is = T, sep = "\t")

ERCCfiles = lapply(list(paste0(root_dir,"K24_ERCCReadsNotCollapsed/abundance.tsv"),
             paste0(root_dir,"K24_ERCCReadsCollapsed/abundance.tsv")), function(x) read.table(x, header = T))
names(ERCCfiles) = c("ERCCNotCollapsed","ERCCCollapsed")

ERCCfiles = lapply(ERCCfiles, function(x) merge(x, ERCC_sq[,c("isoform","chrom")], by.x = "target_id", by.y = "isoform"))
ggplot(ERCCfiles$ERCCNotCollapsed,aes(x = chrom, y = est_counts)) + geom_point() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_ERCC <- function(dat, title){
  ggplot(dat,aes(x = chrom, y = est_counts)) + geom_point() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(y = "Estimated Counts", title = title)
}

plot_grid(plot_ERCC(ERCCfiles$ERCCNotCollapsed,"ERCC Not Collapsed"),plot_ERCC(ERCCfiles$ERCCCollapsed,"ERCC Collapsed"))
                 
length(unique(ERCCfiles$ERCCNotCollapsed$chrom))


## correlation of FL reads and length


mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 14,family="ArialMT"),
                 axis.text.x= element_text(size=12,  family="ArialMT"),
                 axis.text.y= element_text(size=12,  family="ArialMT"))

venn_diagram_plot_twocircles <- function(set1, set2, label_set1, label_set2){
  
  p <- venn.diagram(
    x = list(set1, set2),
    category.names = c(paste(label_set1,"Cortex"),paste(label_set2,"Cortex")),
    filename = NULL,
    output=TRUE,
    
    # Circles
    lwd = 0.2,
    lty = 'blank',
    fill = c("#FDCDAC","#B3E2CD"),
    
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "ArialMT",
    
    # Set names
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "ArialMT",
    #rotation = 1, 
    main = "\n\n\n\n",
    
    print.mode = "raw"
  )
  
  return(p)
  
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
                            x = 0.05, y = 0.97, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 20, fontface = "italic")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  print(paste0("corr.value", corr.value))
  print(paste0("p.value", p.value))
  
  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    #stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point() +
    #scale_fill_distiller(palette=4, direction=1, name = "Density") +
    theme_bw() +
    labs(x = x_lab, y = y_lab, title = paste(title)) + 
    #geom_smooth(method=lm, colour = "black") + 
    mytheme + 
    theme(legend.position = "none") + guides(fill=FALSE)
  
  return(p)
}



ERCCNotCollapsed = ERCCfiles$ERCCNotCollapsed[log10(ERCCfiles$ERCCNotCollapsed$est_counts) > 2,c("target_id")]
ERCC_Collapsed = ERCCfiles$ERCCCollapsed[log10(ERCCfiles$ERCCCollapsed$est_counts) > 2,c("target_id")]

# ERCC Not Collapsed picked up the same ERCCs as collapsed 
plot_grid(venn_diagram_plot_twocircles(ERCCNotCollapsed,ERCC_Collapsed,"ERCC Not Collapsed","ERCC Collapsed"))


## correlation of FL reads and length
corr_plot <- function(ERCC){
  p <- ERCCfiles$ERCCNotCollapsed %>% filter(chrom == ERCC) %>% 
  density_plot(.,"length","est_counts","Length","Estimated Counts",ERCC)
  
  return(p)
}

multipledetected_ERCC = ERCCfiles$ERCCNotCollapsed %>% group_by(chrom) %>% tally() %>% filter(n > 3)

#for(i in multipledetected_ERCC$chrom){print(corr_plot(i))}



