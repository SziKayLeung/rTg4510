load("/gpfs/mrc0/projects/Research_Project-MRC148213/EmmaW/DementiaMouseRRBSmethylation/rTg4510_RRBSbetasComplete.RData")
#load("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP/RRBSmergedThatCorrelate.Rdata")

library("dplyr")
library("ggvenn")
library("ggplot2")
library("stringr")


RRBS_Phenotype= read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/DMP/Tg4510_phenotype_RRBS.csv")

# Targeted Transcriptome
targetgene <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/4_DMP/Mouse_TargeGenes_Coordinates.csv")
DMPgene <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/4_DMP/Mouse_DMPGenes_Coordinates.csv")
targetgene = targetgene %>% mutate(chr = word(coordinates, c(3), sep = fixed(" ")), 
                                   start = as.numeric(as.character(word(coordinates, c(6), sep = fixed(" ")))), 
                                   end = as.numeric(as.character(word(coordinates, c(7), sep = fixed(" ")))),
                                   strand = ifelse(start < end, "Forward","Backward"))

DMPgene = DMPgene %>% mutate(chr = word(coordinates, c(2), sep = fixed(" ")), 
                                   start = as.numeric(as.character(word(coordinates, c(4), sep = fixed(" ")))), 
                                   end = as.numeric(as.character(word(coordinates, c(6), sep = fixed(" ")))),
                                   strand = ifelse(start < end, "Forward","Backward"))


RRBS = RRBS_completebetas %>% tibble::rownames_to_column(., var = "Position") %>% mutate(Chromosome = word(Position,c(1), sep = fixed(":")),
                                                                                    Position = word(Position,c(2), sep = fixed(":")))

rTg4510_samples <- RRBS_Phenotype$Sample_ID
RRBS_Tg4510 = RRBS %>% select(Chromosome, Position, matches(paste(rTg4510_samples, collapse = "|")))

targetgene$gene = str_replace_all(targetgene$gene, fixed(" "), "")
num_gene_probes <- function(input_data,gene){
  gchr = paste0("chr",input_data[input_data$gene == gene, "chr"])
  # 100bp wobble 
  gstart = input_data[input_data$gene == gene, "start"] - 3000
  gend = input_data[input_data$gene == gene, "end"] + 3000
  
  print(gstart)
  print(gend)
  dat = RRBS_Tg4510 %>% filter(Chromosome == gchr) %>% filter(Position > gstart & Position < gend)
  return(dat)
}

numprobes = lapply(targetgene$gene, function(x) num_gene_probes(target_gene,x))
names(numprobes) = targetgene$gene                                        
uscsgenome = lapply(numprobes, function(x) x %>% mutate(Position2 = Position) %>% select(Chromosome,Position,Position2))

# DMP 
DMP_numprobes = lapply(DMPgene$gene, function(x) num_gene_probes(DMPgene,x))
names(DMP_numprobes) = DMPgene$gene  
DMP_uscsgenome = lapply(DMP_numprobes, function(x) x %>% mutate(Position2 = Position) %>% select(Chromosome,Position,Position2))
write.table(DMP_uscsgenome$Ank1,"ANK1_DMP.txt",quote = F, row.names = F, col.names = F)

# subset into chromosome 18 
RRBS_chr18 <- RRBS %>% filter(Chromosome == "chr18")
