#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Comparison of jones2024 dataset
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------------------------

dirnames$jones2024 <- "/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/0_otherStudies/Jones2024/"

## ---------- classification file -----------------

# sqanti classification file from jones
class.names.files$jones <- paste0(dirnames$jones2024, "sqanti/chr_annotations_classification.txt")
class.files$jones <- read.table(class.names.files$jones, as.is = T, sep = "\t", header = T)
class.files$jones_targeted <- class.files$jones %>% filter(associated_gene %in% TargetGene)
message("number of transcripts in jones dataset:", nrow(class.files$jones_targeted))

class.names.files$jonesfiltered <- paste0(dirnames$jones2024, "sqanti/chr_annotations_RulesFilter_result_classification.txt")
class.files$jonesfiltered <- SQANTI_class_preparation(class.names.files$jonesfiltered,"nstandard")
class.files$jonesfiltered_targeted <- class.files$jonesfiltered %>% filter(associated_gene %in% TargetGene)
message("number of transcripts in jones dataset:", nrow(class.files$jonesfiltered_targeted))

# non-canonical junctions
filteredReasons <- read.table(paste0(dirnames$jones2024, "sqanti/chr_annotations_filtering_reasons.txt"), sep = "\t", as.is = T, header = T)
filteredReasons %>% filter(isoform %in% class.files$jones_targeted$isoform)

phenotype$jones <- read.csv(paste0(dirnames$jones2024,"comparison/sample_collection_metadata.csv")) 
tCounts <- read.table(paste0(dirnames$jones2024, "data_minus_bam/nextflow/bambu/counts_transcript.txt"), header = T)
tCountsMelt <- tCounts %>% filter(TXNAME %in% class.files$jones_targeted$isoform) %>% reshape2::melt() %>% 
  mutate(sample_id = word(variable,c(1),sep=fixed("_")))
tCountsMelt <- merge(tCountsMelt,phenotype$jones[,c("sample_id","mouse_id","sex","tissue")], by = "sample_id")
tCountsMelt <- merge(tCountsMelt, class.files$jones_targeted[,c("isoform","associated_gene","associated_transcript","structural_category")], by.x = "TXNAME", by.y = "isoform")
ggplot(tCountsMelt, aes(x = TXNAME, y = value, colour = sex)) + geom_boxplot() + facet_grid(tissue~associated_gene, scales = "free", space = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



## ---------- gffCompare output -----------------

# rTg4510 as reference
rTg4510Ref <- read.table(paste0(dirnames$jones2024,"comparison/rTg4510Jones.chr_annotations_corrected.gtf.tmap"), header = T)
message("Number of detected transcripts: ", nrow(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$jones_targeted$isoform,]))
table(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$jones_targeted$isoform,"class_code"])

# jones dataset as reference
jonesRef <- read.table(paste0(dirnames$jones2024, "comparison/JonesrTg4510.all_iso_ont_collapsed.filtered_counts_filtered.gtf.tmap"), header = T)
detected <- table(jonesRef[jonesRef$ref_id %in% class.files$jones_targeted$isoform, "class_code"])

rTg4510Ref[rTg4510Ref$qry_id %in% c(as.character(unique(tCountsMelt[tCountsMelt$value > 10,"TXNAME"]))),]
