# Script borowed from E.Walker 
# Mouse Array Data: Methylation from normalised beta 

suppressMessages(library("dplyr"))
suppressMessages(library("readr"))
suppressMessages(library("stringr"))

# normalised betas
load("/gpfs/mrc0/projects/Research_Project-MRC148213/Mouse_Array_Data/Normalised_Data_Sesame_2.rdat")
#add in annotation info
probe_mapping_file <- read_csv("/gpfs/mrc0/projects/Research_Project-MRC148213/Mouse_Array_Data/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #24048 probes 
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- species_probes$MusMusculus
species_probes$Bp <- str_split(species_probes$MusMusculus, ":", simplify = T)[,2]
species_probes <- species_probes[,-2]
species_probes$Position <- paste("chr", paste(species_probes$Chr, species_probes$Bp, sep = ":"), sep ="")
identical(species_probes$probeID, row.names(Normalised_Sesame_Betas)) #[1] TRUE
Normalised_Sesame_Betas_annotated <- cbind(Normalised_Sesame_Betas, species_probes)


# check for number of probes within target genes
targetgene <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/4_DMP/Mouse_TargeGenes_Coordinates.csv")
targetgene = targetgene %>% mutate(chr = word(coordinates, c(3), sep = fixed(" ")), 
                                   start = as.numeric(as.character(word(coordinates, c(5), sep = fixed(" ")))), 
                                   end = as.numeric(as.character(word(coordinates, c(7), sep = fixed(" ")))),
                                   strand = ifelse(start < end, "Forward","Backward"))

Normalised_Sesame_Betas_annotated = Normalised_Sesame_Betas_annotated %>% mutate(Position_Base = as.numeric(as.character(Bp)))

targetgene$gene = str_replace_all(targetgene$gene, fixed(" "), "")
num_gene_probes <- function(gene){
  gchr = targetgene[targetgene$gene == gene, "chr"]
  # 100bp wobble 
  gstart = targetgene[targetgene$gene == gene, "start"] - 100
  gend = targetgene[targetgene$gene == gene, "end"] + 100
  
  dat = Normalised_Sesame_Betas_annotated %>% filter(Chr == gchr) %>% filter(Position_Base > gstart & Position_Base < gend)
  return(dat)
}


numprobes = lapply(targetgene$gene, function(x) num_gene_probes(x))
names(numprobes) = targetgene$gene

all_numprobes = bind_rows(numprobes, .id = "gene") 
targetprobes = all_numprobes %>% group_by(gene) %>% tally()
probepos = all_numprobes %>% select(gene,probeID,Position)                   

