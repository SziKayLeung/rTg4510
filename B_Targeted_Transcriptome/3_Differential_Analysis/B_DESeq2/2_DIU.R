## ---------- Script -----------------
##
## Purpose: perform differential isoform usage analysis on mouse rTg4510 ONT and Iso-Seq targeted datasets 
## Adapted tappAS DIU analysis scripts
## https://github.com/ConesaLab/tappAS/blob/master/scripts/DIU.R
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- source functions -----------------

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_usage.R"))


## ---------- input -----------------

# directory names
root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
dirnames <- list(
  # targeted sequencing (Iso-Seq, ONT)
  targ_root = paste0(root_dir, "rTg4510/G_Merged_Targeted/B_cupcake_pipeline"),
  targ_iso_metadata = paste0(root_dir,"rTg4510/0_metadata/B_isoseq_targeted"),
  targ_ont_metadata = paste0(root_dir,"rTg4510/0_metadata/F_ont_targeted"),
  targ_output = paste0(root_dir, "/rTg4510/01_figures_tables/Targeted_Transcriptome")
)

# classification files
sqname="/3_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification"
class.files <- list(
  targ_all = SQANTI_class_preparation(paste0(dirnames$targ_root, sqname,".targetgenes_counts.txt"),"ns"),
  target_filtered = SQANTI_class_preparation(paste0(dirnames$targ_root, sqname,".targetgenes_counts_filtered.txt"),"ns")
)

# DESEQ results
TargetedDESeq <- list(
  ontResTranAnno = readRDS(file = paste0(dirnames$targ_output, "/Ont_DESeq2TranscriptLevel.RDS")),
  isoResTranAnno = readRDS(file = paste0(dirnames$targ_output, "/IsoSeq_DESeq2TranscriptLevel.RDS"))
) 

# expression
Exp <- list(
  targ_ont = list(
    raw = class.files$targ_all %>% dplyr::select(associated_gene, contains("ONT")) %>% select(!contains("sum_FL")),
    norm = TargetedDESeq$ontResTranAnno$wald$norm_counts %>% 
      dplyr::select(sample,associated_gene,normalised_counts, isoform) %>% spread(., sample, value = normalised_counts) %>% 
      remove_rownames %>% tibble::column_to_rownames(var="isoform"),
    normAll = TargetedDESeq$ontResTranAnno$wald$norm_counts_all %>% 
      dplyr::select(sample,associated_gene,normalised_counts, isoform) %>% spread(., sample, value = normalised_counts) %>% 
      remove_rownames %>% tibble::column_to_rownames(var="isoform")
  ),
  targ_iso = list(
    raw = class.files$targ_all %>% dplyr::select(associated_gene, contains("Iso.Seq")) %>% select(!contains("sum_FL")),
    norm = TargetedDESeq$isoResTranAnno$wald$norm_counts %>% 
      dplyr::select(sample,associated_gene,normalised_counts, isoform) %>% spread(., sample, value = normalised_counts) %>% 
      remove_rownames %>% tibble::column_to_rownames(var="isoform"),
    normAll = TargetedDESeq$isoResTranAnno$wald$norm_counts_all %>% 
      dplyr::select(sample,associated_gene,normalised_counts, isoform) %>% spread(., sample, value = normalised_counts) %>% 
      remove_rownames %>% tibble::column_to_rownames(var="isoform")
  )
)

# phenotype
phenotype <- list(
  targeted_rTg4510_ont = read.csv(paste0(dirnames$targ_ont_metadata, "/ONTBarcoded_Phenotype.csv")) %>% 
    mutate(group = Phenotype, time = Age,
           sample = paste0("ONT_",sample),col = paste0(sample,"_",Phenotype)),
  targ_iso = read.table(paste0(dirnames$targ_iso_metadata, "/TargetedMouse_PhenotypeTAPPAS.txt"), header = T) %>% mutate(variable = paste0(sample))
)

# factors
factorsInput <- list(
  targ_ont_WTTG = phenotype$targeted_rTg4510_ont %>% dplyr::select(sample, Phenotype) %>% magrittr::set_rownames(.$sample) %>% 
    dplyr::rename("Replicate" = "Phenotype") %>% 
    mutate(Replicate = ifelse(Replicate == "WT",1,2)),
  targ_ont_age = phenotype$targeted_rTg4510_ont %>% mutate(grouptime = paste0(group,time)) %>% 
    dplyr::select(sample, grouptime) %>% magrittr::set_rownames(.$sample) %>% 
    mutate(Replicate = match(.$grouptime, unique(.$grouptime))),
  
  targ_iso_WTTG = phenotype$targ_iso %>% mutate(sample = str_replace(.$sample, "FL.","Iso.Seq_")) %>%
    magrittr::set_rownames(.$sample) %>%
    mutate(Replicate = match(.$group, unique(.$group)))
) 
factorsInput <- lapply(factorsInput, function(x) x[order(x$Replicate),, drop = FALSE])



## ---------- runDIU -----------------


resultsDIU <- list(
  ontDIUGeno = runDIU(transMatrixRaw=Exp$targ_ont$raw,transMatrix=Exp$targ_ont$normAll,classf=class.files$targ_all,
                   myfactors=factorsInput$targ_ont_WTTG,filteringType="FOLD",filterFC=2),
  ontDIUPath = runDIU(transMatrixRaw=Exp$targ_ont$raw,transMatrix=Exp$targ_ont$normAll,classf=class.files$targ_all,
                   myfactors=factorsInput$targ_ont_age,filteringType="FOLD",filterFC=2),
  isoDIUGeno = runDIU(transMatrixRaw=Exp$targ_iso$raw,transMatrix=Exp$targ_iso$normAll,classf=class.files$targ_all,
         myfactors=factorsInput$targ_iso_WTTG,filteringType="FOLD",filterFC=2)
)

## ---------- write output -----------------

saveRDS(resultsDIU, file = paste0(dirnames$targ_output, "/resultsDIU.RDS"))

for(i in 1:3){
  write.csv(resultsDIU[[i]], paste0(dirnames$targ_output,"/",names(resultsDIU)[[i]],".csv"),row.names=F)
}
