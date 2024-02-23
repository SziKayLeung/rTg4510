
# plot the number of transcripts unique and common to DN and NeuN
p3 <- class.files$targ_sorted %>% group_by(associated_gene, dataset) %>% tally() %>% 
  ggplot(., aes(x = associated_gene, y = n, fill = dataset)) + geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Gene", y = "Number of isoforms")


NeuNrawCounts <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/2_cutadapt_merge/NeuN/read_numbers.txt")
DNrawCounts <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/2_cutadapt_merge/DN/read_numbers.txt")
libraryPrep <- read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/0_metadata/libraryPrepMolarity.csv")
phenotype <- read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/0_metadata/SCNPhenotype.csv")

head(phenotype)
head(libraryPrep)

libraryPrep <- merge(phenotype, libraryPrep, by = "tissue")
RawCounts <- rbind(NeuNrawCounts %>% mutate(cell = "NeuN"), DNrawCounts %>% mutate(cell = "DN")) %>% mutate(barcode = word(V1,c(1),sep=fixed("_")))
RawCounts <- merge(RawCounts, libraryPrep, by = "barcode") %>% dplyr::rename("totalReads" = "V2")
library(scales)
p1 <- ggplot(RawCounts, aes(x = genotype, y = totalReads, colour = age)) + 
  geom_point(size = 3) + 
  scale_y_continuous(labels = label_comma()) + facet_grid(~cell) +
  labs(x = "Genotype", y = "Total ONT raw reads") 
  
p2 <- ggplot(RawCounts, aes(y = totalReads, x = fmol, colour = cell)) + geom_point(size = 3) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Amount in library preparation (fmol)", y = "Total ONT raw reads") 

plot_grid(plot_grid(p1,p2),p3, nrow = 2, labels = c("A","B","C"))
