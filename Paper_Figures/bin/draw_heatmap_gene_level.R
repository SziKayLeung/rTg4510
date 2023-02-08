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
