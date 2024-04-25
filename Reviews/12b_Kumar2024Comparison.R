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


## ---------- gffCompare output - Targeted dataset -----------------

# rTg4510 as reference
rTg4510Ref <- read.table(paste0(dirnames$kumar2024,"/comparison/targeted/rTg4510Kumar.Kumar2024.filtered.gtf.tmap"), header = T)
message("Number of detected transcripts: ", nrow(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$kumar_targeted$isoform,]))
message("Number of exact match transcripts:", nrow(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$kumar_targeted$isoform & rTg4510Ref$class_code == "=",]))
message("Number of partial transcripts:", nrow(rTg4510Ref[rTg4510Ref$qry_id %in% class.files$kumar_targeted$isoform & rTg4510Ref$class_code == "c",]))


# Kumar dataset as reference
KumarRef <- read.table(paste0(dirnames$kumar2024, "comparison/targeted/KumarrTg4510.all_iso_ont_collapsed.filtered_counts_filtered.gtf.tmap"), header = T)
detected <- KumarRef[KumarRef$ref_id %in% class.files$kumar_targeted$isoform & tmap$class_code == "=", "qry_id"]
(nrow(class.files$targ_filtered) - length(detected))/nrow(class.files$targ_filtered)


## ---------- gffCompare output - Whole dataset -----------------

# rTg4510 as reference
rTg4510Ref <- read.table(paste0(dirnames$kumar2024,"/comparison/whole/rTg4510Kumar.Kumar2024.filtered.gtf.tmap"), header = T)
message("Number of exact match or partial transcripts:", length(unique(rTg4510Ref[rTg4510Ref$class_code %in% c("=","c"),"qry_id"])))
3659/nrow(class.files$glob_iso) * 100


# Kumar dataset as reference
KumarRef <- read.table(paste0(dirnames$kumar2024, "comparison/whole/KumarrTg4510.WholeIsoSeq.collapsed.filtered.gtf.tmap"), header = T)
length(unique(KumarRef$qry_id)) == nrow(class.files$glob_iso)
message("Number of exact match or partial transcripts:", length(unique(KumarRef[KumarRef$class_code %in% c("=","c"),"qry_id"])))
Detected <- KumarRef[KumarRef$class_code %in% c("=","c"),"qry_id"]
class.files$glob_iso <- class.files$glob_iso %>% mutate(Kumar = factor(ifelse(isoform %in% Detected,TRUE,FALSE), levels = c(TRUE,FALSE)))
class.files$glob_iso$meanFL <- class.files$glob_iso %>% select(contains("FL.")) %>% apply(., 1, mean)

p1 <- ggplot(class.files$glob_iso, aes(x = Kumar, y = log10(meanFL))) + geom_boxplot() + 
  theme_classic() + labs(y = "log10 mean FL reads", x = "Detected in Kumar2024 dataset") 

t.test(meanFL ~ Kumar, data = class.files$glob_iso)
nrow(class.files$kumar)
