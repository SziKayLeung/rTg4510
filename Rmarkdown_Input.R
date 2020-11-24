# packages
suppressMessages(library(ggplot2))
suppressMessages(library(wesanderson))
suppressMessages(library(dplyr))
suppressMessages(library(knitr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(magicfor))
suppressMessages(library(extrafont))
magic_for(print, silent = TRUE) # call magic_for()
magic_for(put) # call magic_for()

#library(extrafont)
#font_install('fontcm')
loadfonts()

# plot theme
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 #legend.position = c(.90, 0.90),
                 #legend.justification = c(1,1),
                 #legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 14,family="CM Roman"),
                 axis.text.x= element_text(size=12,family="CM Roman"),
                 axis.text.y= element_text(size=12,family="CM Roman"))

# plot label colour
label_colour <- function(genotype){
  if(genotype == "WT"){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "TG"){colour = wes_palette("Royal1")[2]}
  }
  return(colour)
}



# output directory
output_plot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Rmarkdown"