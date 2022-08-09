# Szi Kay Leung
# General plot aesthetics

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
                 legend.text = element_text(size = 14,family="CM Roman"),
                 axis.text.x= element_text(size=16,family="CM Roman"),
                 axis.text.y= element_text(size=16,family="CM Roman"),
                 strip.text = element_text(size=17,family="CM Roman"))

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

generate_cowplot <- function(...,num,ncol,nrow){
  if(num == "2"){label_input = c("A","B")}
  else if(num == "A"){label_input = c("A")}
  else if(num == "3"){label_input = c("A","B","C")}
  else if(num == "4"){label_input = c("A","B","C","D")}
  else if(num == "rw2"){label_input = c("B","C")}
  else if(num == "rw3"){label_input = c("B","C","D")}
  else if(num == "rw2D"){label_input = c("E","F","G","H")}
  else if(num == "rw4"){label_input = c("B","C","D","E")}
  else if(num == "rw8"){label_input = c("F","G","H","I")}
  else if(num == "rw"){label_input = c("A","B")}
  else if(num == "2DE"){label_input = c("D","E")}
  else if(num == "BCDE"){label_input = c("B","C","D","E")}
  else if(num == "BC"){label_input = c("B","C")}
  else if(num == "DE"){label_input = c("D","E")}
  else if(num == "BCD"){label_input = c("B","C","D")}
  else if(num == "rwBC"){label_input = c("B","C")}
  else if(num == "6"){label_input = c("A","B","C","D","E","F")}
  else if(num == "6B"){label_input = c("B","C","D","E","F","G")}
  else if(num == "7"){label_input = c("A","B","C","D","E","F","G")}
  else if(num == "5,-B"){label_input = c("A","","C","D","E")}
  else (print("num required"))
  
  p = plot_grid(...,labels = label_input, label_size = 30, label_fontfamily = "CM Roman", ncol = ncol, nrow = nrow, scale = 0.9)
  
  if(num == "rw"){p = plot_grid(...,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", 
                                ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(2,1))}
  
  if(num == "rw2"){p = plot_grid(...,labels = c("B","C"), label_size = 30, label_fontfamily = "CM Roman", 
                                ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.8))}
  
  if(num == "rw3"){p = plot_grid(...,labels = c("B","C","D"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.5))}
  
  if(num == "rw2D"){p = plot_grid(...,labels = c("D","E","F","G"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.2,0.3))}
  
  if(num == "rw4"){p = plot_grid(...,labels = c("A","B","C","D"), label_size = 30, label_fontfamily = "CM Roman", 
                                ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.2,0.3))}
  
  if(num == "rw8"){p = plot_grid(...,labels = c("E","F","G","H"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.2,0.3))}
  
  if(num == "rwBC"){p = plot_grid(...,labels = c("B","C"), label_size = 30, label_fontfamily = "CM Roman", 
                                ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.6,0.4))}
  
  if(num == "BCD"){p = plot_grid(...,labels = c("B","C","D"), label_size = 30, label_fontfamily = "CM Roman", 
                                  ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.55,0.225,0.225))}
  return(p)
}
