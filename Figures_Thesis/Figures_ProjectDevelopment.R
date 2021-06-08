# Szi Kay Leung
# Script for Introduction and Project Dev Chapter 

# plot theme
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


suppressMessages(library("ggplot2"))
suppressMessages(library("cowplot"))
suppressMessages(library("extrafont"))
suppressMessages(loadfonts())

output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis"

### Functions ########################

# plot ERCC Length 
plot_ERCCLength <- function(){
  ERCC_length <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/ERCC/ERCC_Stats.csv", as.is = T)
  p <- ggplot(ERCC_length, aes(x = Length ))+ geom_histogram(color="black", fill="lightblue") + mytheme + labs(x = "ERCC Length (nucleotides)", y = "Frequency") +  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 10, l = 0)))
  return(p)
}

pdf (paste0(output_plot_dir,"/ProjectDevelopment.pdf"), width = 10, height = 15)
plot_grid(plot_ERCCLength(),NULL,NULL,NULL,scale = 0.9)
dev.off()
