#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Comparison of Kumar2024 dataset
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------------------------

dirnames$kumar2024 <- "/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/0_otherStudies/Kumar2024/"

## ---------- classification file -----------------

# sqanti classification file from kumar
class.names.files$kumar <- paste0(dirnames$kumar2024, "/cupcake/Kumar2024_RulesFilter_result_classification.txt")
class.files$kumar <- SQANTI_class_preparation(class.names.files$kumar,"nstandard")
class.files$kumar_targeted <- class.files$kumar %>% filter(associated_gene %in% TargetGene)
message("number of transcripts in Kumar dataset:", nrow(class.files$kumar_targeted))


## ---------- gffCompare output -----------------

# rTg4510 as reference
rTg4510Ref <- read.table(paste0(dirnames$kumar2024,"/comparison/rTg4510Kumar.Kumar2024.filtered.gtf.tmap"), header = T)
message("Number of detected transcripts: ", nrow(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$kumar_targeted$isoform,]))
message("Number of exact match transcripts:", nrow(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$kumar_targeted$isoform & rTg4510Ref$class_code == "=",]))
message("Number of partial transcripts:", nrow(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$kumar_targeted$isoform & rTg4510Ref$class_code == "c",]))


# Kumar dataset as reference
KumarRef <- read.table(paste0(dirnames$kumar2024, "comparison/KumarrTg4510.all_iso_ont_collapsed.filtered_counts_filtered.gtf.tmap"), header = T)
detected <- KumarRef[KumarRef$ref_id %in% class.files$kumar_targeted$isoform & tmap$class_code == "=", "qry_id"]
(nrow(class.files$targ_filtered) - length(detected))/nrow(class.files$targ_filtered)
