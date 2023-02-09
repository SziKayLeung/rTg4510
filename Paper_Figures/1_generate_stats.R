## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##



## ---------- Source function and config files -----------------

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "rTg4510_config.R"))
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(SC_ROOT,"bin/draw_heatmap_gene_level.R"))


# mean number of genes
cat("Mean number of genes:", mean(sapply(sub_class.files, function(x) length(unique(x[["associated_gene"]])))))
cat("Mean number of genes in WT samples:", mean(sapply(sub_class.files[WT], function(x) length(unique(x[["associated_gene"]])))))
cat("Mean number of genes in TG samples:", mean(sapply(sub_class.files[TG], function(x) length(unique(x[["associated_gene"]])))))


# mean number of isoforms
cat("Mean number of isoforms:", mean(sapply(sub_class.files, function(x) nrow(x))))
cat("Mean number of isoforms in WT samples:", mean(sapply(sub_class.files[WT], function(x) nrow(x))))
cat("Mean number of isoforms in TG samples:", mean(sapply(sub_class.files[TG], function(x) nrow(x))))

