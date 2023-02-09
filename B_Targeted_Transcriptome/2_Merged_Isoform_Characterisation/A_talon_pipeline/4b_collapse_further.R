## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Date: Jan 2023
## 
## Explore the output from cupcake collapse after parameter optimisation of --max_5_diff and max_3_diff
## Output of 4_collapse_further.sh
## Format of group.txt, generated from cupcake collapse (v8.6)
##      PB.1.1      isoform1,isoform2,isoform3
##      PB.1.2      isoform4,isoform5 ...
## ----------------------------------


## ---------- library -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("cowplot"))
suppressMessages(library("VennDiagram"))


## ---------- input -----------------

# directory names
dirnames <- list(
  coll = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/4_characterise/CollapseMore/"
)

# input collapse output from parameter optimisation
input.collapse <- list(
  default = read.table(paste0(dirnames$coll, "1000_5diff_100_3diff/IsoSeqONT.collapsed.group.txt")),
  d50 = read.table(paste0(dirnames$coll, "50_5diff_50_3diff/IsoSeqONT.collapsed.group.txt")),
  d5 = read.table(paste0(dirnames$coll, "5_5diff_5_3diff/IsoSeqONT.collapsed.group.txt")),
  d10 = read.table(paste0(dirnames$coll, "10_5diff_10_3diff/IsoSeqONT.collapsed.group.txt")),
  d20 = read.table(paste0(dirnames$coll, "20_5diff_20_3diff/IsoSeqONT.collapsed.group.txt"))
)


## ---------- functions -----------------

# Aim: plot the number of transcripts collapsed 
# Input:  
  # param5 = str: integer of 5_diff denoted in _group.txt. file
  # param3 = str: integer of 3_diff denoted in _group.txt. file
# Output:
  # bar-plot
hist_n_transcripts <- function(param5, param3){
  
  # read collapsed_group file based on param5 and param3
  collapse <- read.table(paste0(dirnames$coll, param5,"_5diff_", param3,"_3diff/IsoSeqONT.collapsed.group.txt"))
  collapse$V2 <- as.character(collapse$V2)
  
  # determine the number of collapsed isoforms from delimiter
  collapse$n_transcript <- sapply(collapse$V2, function(x) ifelse(grepl(",", x), lengths(gregexpr(",", x)) + 1, 1))
  
  print(nrow(collapse))
  
  # bar plot of the number of transcrips collapsed
  p <- collapse %>% group_by(n_transcript) %>% tally() %>% 
    ggplot(., aes(x = n_transcript, y = log10(n))) + geom_bar(stat = "identity") +
    labs(x = "Number of transcripts collapsed", y = "Frequency", title = paste0("max_5_diff:", param5, "\nmax_3_diff:", param3)) + 
    theme_bw()

  return(p)
}


## ---------- plots -----------------

# plot the frequency of transcripts collapsed across the different parameters
pColl <- list(
  default = hist_n_transcripts(1000,100),
  d5 = hist_n_transcripts(5,5),
  d10 = hist_n_transcripts(10,10),
  d20 = hist_n_transcripts(20,20)
)

plot_grid(pColl$d5,pColl$d10, pColl$d20, pColl$default)

ggvenn(
  list(d20 = input.collapse$d20$V2, d5 = input.collapse$d5$V2, d50 = input.collapse$d50$V2, default = input.collapse$default$V2), 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)