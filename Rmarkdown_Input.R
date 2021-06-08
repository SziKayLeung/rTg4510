# packages
suppressMessages(library(ggplot2))
suppressMessages(library(wesanderson))
suppressMessages(library(dplyr))
suppressMessages(library(knitr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(grid))
suppressMessages(library(VennDiagram))
#suppressMessages(library(magicfor))
#suppressMessages(library(extrafont))
#magic_for(print, silent = TRUE) # call magic_for()
#magic_for(put) # call magic_for()

#library(extrafont)
#font_install('fontcm')
#loadfonts()

# plot theme
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=20,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.8, colour = "black"),
                 axis.title.y = element_text(vjust=0.8, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 #legend.position = c(.90, 0.90),
                 #legend.justification = c(1,1),
                 #legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 20,family="CM Roman"),
                 axis.text.x= element_text(size=16,family="CM Roman"),
                 axis.text.y= element_text(size=16,family="CM Roman"),
                 plot.margin = unit(c(2,2,2,2), "cm"))

# plot label colour
label_colour <- function(genotype){
  if(genotype == "WT"){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "TG"){colour = wes_palette("Royal1")[2]}else{
      if(genotype == "Merged"){colour = wes_palette("Royal1")[4]
    }
    }
  }
  return(colour)
}

<<<<<<< HEAD:Rmarkdown_Input.R
<<<<<<< HEAD:Rmarkdown_Input.R
# To scale axis into 1000s
ks <- function(x){ format(x/1000, big.mark=",")} 

=======
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:Rmarkdown_Input.R
=======
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:Rmarkdown_Input.R
#
input_multiple_csv <- function(input_dir_path, str_pattern){
  dat_files <- list.files(paste0(input_dir_path), pattern = str_pattern, full.names = T)
  dat <- lapply(dat_files, function(x) read.csv(x))
  names(dat) <- list.files(paste0(input_dir_path), pattern = str_pattern)

  for(i in 1:length(names(dat))){
    dat[[i]]$Sample <- paste(word(names(dat)[[i]], c(1), sep = fixed ('.')))
  }

  return(dat)
}

# different types of plots
plot_boxplot <- function(dat, x_var, y_var, xlabel, ylabel, group_var){

  x_var <- rlang::sym(quo_name(enquo(x_var)))
  y_var <- rlang::sym(quo_name(enquo(y_var)))
  group_var <- rlang::sym(quo_name(enquo(group_var)))

  p <- ggplot(dat, aes(x = reorder( !! x_var, - !!y_var), y = !!y_var)) +
    geom_boxplot() +
    geom_jitter(aes(color = !!group_var), shape=16, position=position_jitter(0.2)) +
    theme_bw() +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = comma) +
    labs(x = xlabel, y = ylabel) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(margin = margin(t = 0, r = 15 , b = 0, l = 0)))


  return(p)

<<<<<<< HEAD:Rmarkdown_Input.R
<<<<<<< HEAD:Rmarkdown_Input.R
}

# read_merge_gtf
# Aim: Extract coordinates with output sqanti2 filtered gtf file and merge with sqanti2 classification.txt using PBID
# Input: sqanti2.classification.gtf and sqanti2_classification file (already read)
# Output: dataframe with 3 columns: V9 from input gtf file, genome coordinates, isoform (pacbio id)
read_merge_gtf <- function(gtf_input, sqanti2_class){
  # read in gff
  input <- read.delim(gtf_input, header=F, comment.char="#") %>%
    # filter only transcripts in column 3
    filter(V3 == "transcript") %>%
    # take only chromosome (V1), start coordinates (V4), end coordinates (v5)
    mutate(gtf_coordinates = paste(V1,":", V4,"-", V5)) %>%
    .[,-c(1:8)]

  # create separate column for merging
  input$isoform <- word(input$V9,2, sep = ";")
  input$isoform <- word(input$isoform,3, sep = " ")
  input$isoform <- gsub('"', '', input$isoform)

  final_merge <- merge(input, sqanti2_class, by = "isoform", all.y = TRUE)
  return(final_merge)
}

read_gtf <- function(gtf_input){
  input <- read.delim(gtf_input, header=F, comment.char="#") %>%
    # filter only transcripts in column 3
    filter(V3 == "transcript") %>%
    # take only chromosome (V1), start coordinates (V4), end coordinates (v5)
    mutate(gtf_coordinates = paste(V1,":", V4,"-", V5)) %>%
    .[,-c(1:8)]

  # create separate column for merging
  input$isoform <- word(input$V9,2, sep = ";")
  input$isoform <- word(input$isoform,3, sep = " ")
  input$isoform <- gsub('"', '', input$isoform)

  return(input)
}

=======
}

# read_merge_gtf
# Aim: Extract coordinates with output sqanti2 filtered gtf file and merge with sqanti2 classification.txt using PBID
# Input: sqanti2.classification.gtf and sqanti2_classification file (already read)
# Output: dataframe with 3 columns: V9 from input gtf file, genome coordinates, isoform (pacbio id)
read_merge_gtf <- function(gtf_input, sqanti2_class){
  # read in gff
  input <- read.delim(gtf_input, header=F, comment.char="#") %>%
    # filter only transcripts in column 3
    filter(V3 == "transcript") %>%
    # take only chromosome (V1), start coordinates (V4), end coordinates (v5)
    mutate(gtf_coordinates = paste(V1,":", V4,"-", V5)) %>%
    .[,-c(1:8)]

  # create separate column for merging
  input$isoform <- word(input$V9,2, sep = ";")
  input$isoform <- word(input$isoform,3, sep = " ")
  input$isoform <- gsub('"', '', input$isoform)

  final_merge <- merge(input, sqanti2_class, by = "isoform", all.y = TRUE)
  return(final_merge)
}

read_gtf <- function(gtf_input){
  input <- read.delim(gtf_input, header=F, comment.char="#") %>%
    # filter only transcripts in column 3
    filter(V3 == "transcript") %>%
    # take only chromosome (V1), start coordinates (V4), end coordinates (v5)
    mutate(gtf_coordinates = paste(V1,":", V4,"-", V5)) %>%
    .[,-c(1:8)]

  # create separate column for merging
  input$isoform <- word(input$V9,2, sep = ";")
  input$isoform <- word(input$isoform,3, sep = " ")
  input$isoform <- gsub('"', '', input$isoform)

  return(input)
}

>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:Rmarkdown_Input.R
=======
}

# read_merge_gtf
# Aim: Extract coordinates with output sqanti2 filtered gtf file and merge with sqanti2 classification.txt using PBID
# Input: sqanti2.classification.gtf and sqanti2_classification file (already read)
# Output: dataframe with 3 columns: V9 from input gtf file, genome coordinates, isoform (pacbio id)
read_merge_gtf <- function(gtf_input, sqanti2_class){
  # read in gff
  input <- read.delim(gtf_input, header=F, comment.char="#") %>%
    # filter only transcripts in column 3
    filter(V3 == "transcript") %>%
    # take only chromosome (V1), start coordinates (V4), end coordinates (v5)
    mutate(gtf_coordinates = paste(V1,":", V4,"-", V5)) %>%
    .[,-c(1:8)]

  # create separate column for merging
  input$isoform <- word(input$V9,2, sep = ";")
  input$isoform <- word(input$isoform,3, sep = " ")
  input$isoform <- gsub('"', '', input$isoform)

  final_merge <- merge(input, sqanti2_class, by = "isoform", all.y = TRUE)
  return(final_merge)
}

read_gtf <- function(gtf_input){
  input <- read.delim(gtf_input, header=F, comment.char="#") %>%
    # filter only transcripts in column 3
    filter(V3 == "transcript") %>%
    # take only chromosome (V1), start coordinates (V4), end coordinates (v5)
    mutate(gtf_coordinates = paste(V1,":", V4,"-", V5)) %>%
    .[,-c(1:8)]

  # create separate column for merging
  input$isoform <- word(input$V9,2, sep = ";")
  input$isoform <- word(input$isoform,3, sep = " ")
  input$isoform <- gsub('"', '', input$isoform)

  return(input)
}

>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:Rmarkdown_Input.R

plot_line <- function(dat, x_var, y_var, xlabel, ylabel, group_var){

  x_var <- rlang::sym(quo_name(enquo(x_var)))
  y_var <- rlang::sym(quo_name(enquo(y_var)))
  group_var <- rlang::sym(quo_name(enquo(group_var)))

  p <- ggplot(dat, aes(x = reorder( !! x_var, - !!y_var), y = !!y_var, group = !! group_var)) +
    geom_line(aes(linetype=!!group_var, color = !! group_var)) +
    geom_point(shape=20) +
    theme_bw() +
    #scale_y_continuous(labels = comma) +
    labs(x = xlabel, y = ylabel) +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(margin = margin(t = 0, r = 15 , b = 0, l = 0)))

  return(p)

}

plot_line_no_reorder <- function(dat, x_var, y_var, xlabel, ylabel, group_var){

  x_var <- rlang::sym(quo_name(enquo(x_var)))
  y_var <- rlang::sym(quo_name(enquo(y_var)))
  group_var <- rlang::sym(quo_name(enquo(group_var)))

  p <- ggplot(dat, aes(x = !! x_var, y = !!y_var, group = !! group_var)) +
    geom_line(aes(linetype=!!group_var, color = !! group_var)) +
    geom_point(shape=20) +
    theme_bw() +
    #scale_y_continuous(labels = comma) +
    labs(x = xlabel, y = ylabel) +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(margin = margin(t = 0, r = 15 , b = 0, l = 0)))

  return(p)

}

# density plot
density_plot <- function(dat,x.var,y.var, x_lab, y_lab,title){

  print(cor.test(dat[[x.var]],dat[[y.var]]))
  cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")

  corr.value <- cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  p.value <- cor.test(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")$p.value


  # corr.value <- cor(FSM_TPM$ISOSEQ_TPM_Normalised,FSM_TPM$RNASeq_TPM) # normalised ISOSEQ FL counts to length
  corr <- grobTree(textGrob(paste("r = ", round(corr.value, 2)),
                            x = 0.05, y = 0.97, hjust = 0,
                            gp = gpar(col = "black", fontsize = 20, fontface = "italic", family="CM Roman")))

  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  print(paste0(title))
  print(paste0("Correlation between", x.var, "and", y.var))
  print(paste0("corr.value", corr.value))
  print(paste0("p.value", p.value))


  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1, name = "Density") +
    theme_bw() +
    labs(x = x_lab, y = y_lab, title = paste(title,"\n\n")) +
    geom_smooth(method=lm, colour = "black") +
    mytheme +
    theme(legend.position = "none")

  return(p)
}

# TOFU_read <file_suffix> <output>
# Aim: Read in the multiple files generated from Cupcake ToFU from testing CCS Parameters
# <file_suffix> = ending of files generated from Cupcake TOFU
# <output> = name of dataframe to be saved to environment
TOFU_read <- function(file_suffix, output){

  # list and name the input files based on the file suffix
  input <- list.files(working_tofu_dir, pattern = paste0(file_suffix,"*$"), full.names = T)
  ifelse(file_suffix == "filtered.abundance.txt" || file_suffix == "abundance.txt",
         input <- lapply(input, function(x) read.table(x, as.is = T, sep = "\t", header = T)),
         input <- lapply(input, function(x) read.table(x)))
  names(input) <- list.files(working_tofu_dir, pattern = paste0(file_suffix,"*$"))

  # bind all the lists into one dataframe
  input <- bind_rows(input, .id = "variable") %>%
    mutate(variable = paste0(word(variable, c(1), sep = fixed(".")),".", word(variable, c(2), sep = fixed("."))))

  if(file_suffix == "collapsed.group.txt"){
    colnames(input)[2] <- "PB.ID"
    #input$count_collapsed <-count.fields(textConnection(input$V2), sep = ",")
  } else if(file_suffix == "collapsed.filtered.gff" || file_suffix == "collapsed.gff"){
    colnames(input)[14] <- "PB.ID"
  } else if(file_suffix == "filtered.abundance.txt" || file_suffix == "abundance.txt"){
    input <- input %>% select("variable","pbid","count_fl")
    colnames(input)[2] <- "PB.ID"
  } else{
    print("file_suffix: collapsed.group.txt OR collapsed.filtered.gff OR filtered.abundance.txt")
  }

  assign(output, input, envir=.GlobalEnv)
}

# heatmap
draw_heatmap <- function(stage){

  ## Heatmap of the number of FL reads associated per ERCC across the samples
  # Create matrix of the heatmap for input
  reads_heatmap <- collapsed_TOFU_reads %>% filter(V3 == "transcript") %>% group_by(variable, V1) %>% tally(count_fl) %>%
    spread(., variable, n) %>%
    column_to_rownames(., var = "V1") %>%
    as.matrix()

  # Create a metadata dataframe for additional information to annotate the heatmap
  reads_metadata <- as.data.frame(unique(collapsed_TOFU_reads$variable)) %>%
    `colnames<-`(c("variable")) %>%
    mutate(num_passes = word(variable, c(2), sep = fixed("_"))) %>%
    mutate(RQ = word(variable, c(4), sep = fixed("_"))) %>%
    mutate(Sample = word(variable, c(1), sep = fixed("_")))

  # only keep the relevant columns
  if(stage == "Part1"){
    ann <- data.frame(reads_metadata$num_passes, reads_metadata$RQ)
    colnames(ann) <- c('Passes', 'RQ')
    colours <- list('Passes' = c('1' = alpha(wes_palette("Darjeeling2")[4],0.7), '2' = wes_palette("Zissou1")[2],'3' =wes_palette("Darjeeling2")[2] ),
                    'RQ' = c('0.8' = wes_palette("Chevalier1")[3], '0.9' = wes_palette("Darjeeling1")[2],
                             '0.95' = wes_palette("Darjeeling1")[1],'0.99' = wes_palette("Darjeeling2")[5]))
  } else if(stage == "Part2"){
    # for ERCC not detected - replace NA as 1, as required value for heatmap clustering
    # note replace with 1 rather than 0 as log10 1 = 0
    reads_heatmap[is.na(reads_heatmap)] <- 1

    ann <- data.frame(reads_metadata$num_passes, reads_metadata$RQ, reads_metadata$Sample)
    # determine the colours and format for the annotations
    colnames(ann) <- c('Passes', 'RQ','Sample')
    colours <- list('Passes' = c('1' = alpha(wes_palette("Darjeeling2")[4],0.7), '2' = wes_palette("Zissou1")[2],'3' =wes_palette("Darjeeling2")[2]),
                    'RQ' = c('0.8' = wes_palette("Chevalier1")[3], '0.9' = wes_palette("Darjeeling1")[2],
                             '0.95' = wes_palette("Darjeeling1")[1],'0.99' = wes_palette("Darjeeling2")[5]),
                    'Sample' = c('K18' = wes_palette("GrandBudapest1")[2],'K23' = wes_palette("GrandBudapest1")[3], 'M21' = wes_palette("GrandBudapest1")[4]))
  } else {
    print("Part1 or Part2 for function input")
  }

  colAnn <- HeatmapAnnotation(df = ann,
                              which = 'col',
                              col = colours,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'))



  ht_list <- Heatmap(log10(reads_heatmap),
                     top_annotation = colAnn,
                     show_row_dend = FALSE, show_column_names = FALSE, col=plasma(100),
                     name = 'FL Count (Log10)',
                     row_names_gp = gpar(fontsize = 10),
                     heatmap_legend_param = list(direction = "horizontal"))

  return(draw(ht_list, heatmap_legend_side = "top", annotation_legend_side = "top",merge_legend = TRUE))
}

venn_diagram_plot_twocircles <- function(set1, set2, label_set1, label_set2){

  # not to generate log file
  flog.threshold(ERROR)


  p <- venn.diagram(
    x = list(set1, set2),
    category.names = c(label_set1,label_set2),
    filename = NULL,
    output=TRUE,

    # Circles
    lwd = 0.2,
    lty = 'blank',
    fill = c("#B3E2CD","#FDCDAC"),

    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "ArialMT",

    # Set names
    cat.cex = 2,
    cat.default.pos = "text",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "ArialMT",
    #rotation = 1,
    main = "\n\n\n\n",

    print.mode = "raw"
  )

  return(p)

}

<<<<<<< HEAD:Rmarkdown_Input.R
<<<<<<< HEAD:Rmarkdown_Input.R
venn_diagram_plot_threecircles <- function(set1, set2, set3, label_set1, label_set2, label_set3){
  
  # not to generate log file
  flog.threshold(ERROR)
  
  
  p <- venn.diagram(
    x = list(set1, set2, set3),
    category.names = c(label_set1,label_set2, label_set3),
    filename = NULL,
    output=FALSE,
    
    # Circles
    lwd = 0.2,
    lty = 'blank',
    fill = c("#B3E2CD","#FDCDAC","#CBD5E8"),
    
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "ArialMT",
    
    # Set names
    cat.cex = 1,
    cat.default.pos = "text",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.055),
    cat.fontfamily = "ArialMT",
    rotation = 1,
    main = "\n\n\n\n",
    
    print.mode = "raw"
  )
  
  return(p)
  
}



=======
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:Rmarkdown_Input.R
=======
>>>>>>> c34b526e4f39c14e6ca6579cf2dd5726cda490a5:Rmarkdown_Input.R
# output directory
output_plot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/QC/Rmarkdown"
