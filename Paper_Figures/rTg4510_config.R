# Szi Kay Leung: sl693@exeter.ac.uk

## --------------------------- 
# variables
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")

samples <- c("K17","K18","K23","K24","L22","M21","O18","O23","Q20","Q21","S18","S23")
TG <- c("K18","K24", "L22","O18","Q20","S18")
WT <- setdiff(samples, TG)

## --------------------------- 
# directory names
root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
dirnames <- list(
  # global transcriptome (Iso-Seq, Iso-Seq + RNA-Seq)
  glob_metadata = paste0(root_dir, "rTg4510/0_metadata/A_isoseq_whole"),
  glob_root = paste0(root_dir, "rTg4510/A_IsoSeq_Whole"),
  glob_rnaseq = paste0(root_dir,"rTg4510/C_RNASeq/1_RNASeq_Isabel"),
  glob_tabroot = paste0(root_dir, "rTg4510/A_IsoSeq_Whole/3_differential/2_Results"),
  glob_output = paste0(root_dir, "/rTg4510/01_figures_tables/Whole_Transcriptome"),
  
  # targeted sequencing (Iso-Seq, ONT)
  targ_iso_metadata = paste0(root_dir,"rTg4510/0_metadata/B_isoseq_targeted"),
  targ_ont_metadata = paste0(root_dir,"rTg4510/0_metadata/F_ont_targeted"),
  targ_iso_root = paste0(root_dir, "/rTg4510/B_IsoSeq_Targeted/thesis_dump/DiffAnalysis_noRNASEQ"),
  targ_ont_root = paste0(root_dir, "/rTg4510/F_ONT_Targeted/thesis_dump/TALON"),
  targ_anno = paste0(root_dir,"/rTg4510/F_ONT_Targeted/thesis_dump/TALON/All/Merged/TargetGenes")
)


## --------------------------- 
# Final classification file
class.names.files <- list(
  glob_iso = paste0(dirnames$glob_root, "/2_post_isoseq3/9_sqanti3/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt"),
  targ_iso = paste0(dirnames$targ_iso_root, "/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"),
  targ_ont = paste0(dirnames$targ_ont_root, "/All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon_classification.txt")
) 
class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

# for downstream subsetting of the global transcriptome by WT and TG mice sample
sub_class.files <- lapply(samples, function(x) subset_class_by_sample(class.files$glob_iso,x))
names(sub_class.files) <- samples


## ---------------------------
# Isabel's supplementary table of differentially expressed genes in rTg4510 
rnaseq_results <- list(
  AgeGenotypeDEG = read.csv(paste0(dirnames$glob_rnaseq, "/Isabel_Supp4_Tg4510AgeGenotypeDEG.csv"), header = T, as.is = T),
  GenotypeDEG = read.csv(paste0(dirnames$glob_rnaseq, "/Isabel_Supp2_Tg4510GenotypeDEG.csv"), header = T, as.is = T)
)


## ---------------------------
# TAPPAS (Differential Analysis) 
tappas_dir <- list(
  glob_iso = paste0(dirnames$glob_tabroot, "/IsoSeq_Expression"), 
  glob_rna = paste0(dirnames$glob_tabroot, "/RNASeq_Expression"),
  targ_iso = paste0(dirnames$targ_iso_root, "/TAPPAS_OUTPUT/IsoSeq_Expression"),
  targ_ont = paste0(dirnames$targ_ont_root, "/TAPPAS_OUTPUT")
)

phenotype <- list(
  glob_iso = read.table(paste0(dirnames$glob_metadata, "/WholeIsoSeq_PhenotypeTAPPAS.txt"), header = T) %>% 
    mutate(variable = paste0(sample)),
  glob_rna = read.table(paste0(dirnames$glob_metadata,"/WholeAllMouse_PhenotypeTAPPAS.txt"), header = T) %>% 
    mutate(col_names = paste0(group,".",sample),
           age = paste0(time, "_mos"), 
           pheno = ifelse(group == "CONTROL", "WT", "TG"), 
           variable = paste0("FL.", sample,"_",pheno,"_",age)),
  
  targ_iso = read.table(paste0(dirnames$targ_iso_metadata, "/TargetedMouse_PhenotypeTAPPAS.txt"), header = T) %>% mutate(variable = paste0(sample)),
  targ_ont = read.table(paste0(dirnames$targ_ont_metadata, "/ONT_phenotype.txt"), header = T) %>% mutate(variable = paste0(sample))
  
)


# Output from Tappas_DEA.R
tappassiggene <- list(
  glob = read_dea_files(paste0(dirnames$glob_tabroot,"/DifferentialGeneExpression_Analysis.xlsx")),
  targ_iso = read_dea_files(paste0(dirnames$targ_iso_root,"/TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx")),
  targ_ont = read_dea_files(paste0(dirnames$targ_ont_root,"/TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx"))
)

tappassigtrans <- list(
  glob = read_dea_files(paste0(dirnames$glob_tabroot,"/DifferentialTransExpression_Analysis.xlsx")),
  targ_iso = read_dea_files(paste0(dirnames$targ_iso_root,"/TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx")),
  targ_ont = read_dea_files(paste0(dirnames$targ_ont_root,"/TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx"))
)
