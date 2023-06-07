# Szi Kay Leung: sl693@exeter.ac.uk

LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen"
source(paste0(LOGEN,"/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"/target_gene_annotation/summarise_gene_stats.R"))


## --------------------------- 
# variables
TargetGene <- c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf")

wholesamples <- c("K17","K18","K23","K24","L22","M21","O18","O23","Q20","Q21","S18","S23")
wholeTG <- c("K18","K24", "L22","O18","Q20","S18")
wholeWT <- setdiff(wholesamples, wholeTG)

targetedWT <- c("K19","K23","K21","K17","S19","M21","O23","P19","Q21","S23","Q17","Q23")
targetedTG <- c("K18","K20","K24","L22","O18","O22","T20","Q20","S18","Q18","L18","T18")

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
  targ_root = paste0(root_dir, "rTg4510/G_Merged_Targeted/B_cupcake_pipeline"),
  targ_iso_root = paste0(root_dir, "rTg4510/B_IsoSeq_Targeted"),
  targ_ont_root = paste0(root_dir, "rTg4510/F_ONT_Targeted"),
  #targ_anno = paste0(root_dir,"rTg4510/F_ONT_Targeted/thesis_dump/TALON/All/Merged/TargetGenes")
  targ_anno = paste0(root_dir,"rTg4510/G_Merged_Targeted/B_cupcake_pipeline/4_characterise/TargetGenes"),
  targ_output = paste0(root_dir, "/rTg4510/01_figures_tables/Targeted_Transcriptome")
)



## ------------- Phenotype files -------------------

phenotype <- list(
  glob_iso = read.table(paste0(dirnames$glob_metadata, "/WholeIsoSeq_PhenotypeTAPPAS.txt"), header = T) %>% 
    mutate(variable = paste0(sample)),
  glob_rna = read.table(paste0(dirnames$glob_metadata,"/WholeAllMouse_PhenotypeTAPPAS.txt"), header = T) %>% 
    mutate(col_names = paste0(group,".",sample),
           age = paste0(time, "_mos"), 
           pheno = ifelse(group == "CONTROL", "WT", "TG"), 
           variable = paste0("FL.", sample,"_",pheno,"_",age)),
  
  targ_iso = read.table(paste0(dirnames$targ_iso_metadata, "/TargetedMouse_PhenotypeTAPPAS.txt"), header = T) %>% mutate(variable = paste0(sample)),
  targ_ont = read.table(paste0(dirnames$targ_ont_metadata, "/ONT_phenotype.txt"), header = T) %>% mutate(variable = paste0(sample)),
  tg4510_samples = read.csv(paste0(dirnames$glob_metadata, "/Tg4510_fullsample.csv"))[,c("Genotype","Age_in_months", "Sample.ID","RIN","ng.ul")]
  
)
phenotype$whole_rTg4510_iso <- phenotype$tg4510_samples %>% filter(Sample.ID %in% wholesamples) %>% mutate(Phenotype = Genotype)
phenotype$targeted_rTg4510_iso <- phenotype$tg4510_samples %>% filter(Sample.ID %in% c(targetedTG, targetedWT)) %>% 
  mutate(group = Genotype, time = Age_in_months,sample = paste0("Iso.Seq_",Sample.ID),col = paste0(Sample.ID,"_",group)) %>% 
  mutate(group = factor(group, levels = c("WT","TG")))
phenotype$targeted_rTg4510_ont <- read.csv(paste0(dirnames$targ_ont_metadata, "/ONTBarcoded_Phenotype.csv"))
phenotype$targeted_rTg4510_ont <- phenotype$targeted_rTg4510_ont %>% mutate(group = Phenotype, time = Age) %>% mutate(group = factor(group, levels = c("WT","TG")))
phenotype$targeted_rTg4510_ont <- phenotype$targeted_rTg4510_ont %>% mutate(sample = paste0("ONT_",sample),col = paste0(sample,"_",Phenotype)) 


## --------------------------- 
# Final classification file
class.names.files <- list(
  glob_iso = paste0(dirnames$glob_root, "/2_post_isoseq3/9_sqanti3/WholeIsoSeq.collapsed_RulesFilter_result_classification.txt"),
  targ_all = paste0(dirnames$targ_root, "/3_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts.txt"),
  targ_filtered = paste0(dirnames$targ_root, "/3_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt"),
  iso_match = paste0(dirnames$targ_iso_root, "/7b_matched_only/8d_sqanti3/MatchedMouse_RulesFilter_result_classification.txt")#,
  #targ_iso = paste0(dirnames$targ_iso_root, "/thesis_dump/DiffAnalysis_noRNASEQ/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"),
  #targ_ont = paste0(dirnames$targ_ont_root, "/thesis_dump/TALON/All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon_classification.txt")
) 
class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

# for downstream subsetting of the global transcriptome by WT and TG mice sample
sub_class.files <- lapply(wholesamples, function(x) subset_class_by_sample(class.files$glob_iso,x))
names(sub_class.files) <- wholesamples

group_class.files <- bind_rows(
  subset_class_phenotype(class.files$glob_iso, phenotype$whole_rTg4510_iso, "WT") %>% mutate(Dataset = "WT"),
  subset_class_phenotype(class.files$glob_iso, phenotype$whole_rTg4510_iso, "TG") %>% mutate(Dataset = "TG")
)

# Expression 
rawExp <- list(
  targ_ont_all = class.files$targ_all %>% dplyr::select(associated_gene, contains("ONT")) %>% select(!contains("sum_FL"))
)

## ---------- DESeq2 results -----------------
cat("Input DESEQ2 results\n")
GlobalDESeq <- list(
  resTranAnno = readRDS(file = paste0(dirnames$glob_output, "/IsoSeq_DESeq2TranscriptLevel.RDS")),
  resGeneAnno = readRDS(file = paste0(dirnames$glob_output, "/IsoSeq_DESeq2GeneLevel.RDS")),
  RresTranAnno = readRDS(file = paste0(dirnames$glob_output, "/RNASeqHybrid_DESeq2TranscriptLevel.RDS")),
  RresGeneAnno = readRDS(file = paste0(dirnames$glob_output, "/RNASeqHybrid_DESeq2GeneLevel.RDS")),
  resGeneComparison = readRDS(file = paste0(dirnames$glob_output, "/Comparison_DESeq2GeneLevel.RDS"))
) 


TargetedDESeq <- list(
  ontResTranAnno = readRDS(file = paste0(dirnames$targ_output, "/Ont_DESeq2TranscriptLevel.RDS")),
  isoResTranAnno = readRDS(file = paste0(dirnames$targ_output, "/IsoSeq_DESeq2TranscriptLevel.RDS")),
  ontResGeneAnno = readRDS(file = paste0(dirnames$targ_output, "/Ont_DESeq2GeneLevel.RDS")),
  isoResGeneAnno = readRDS(file = paste0(dirnames$targ_output, "/IsoSeq_DESeq2GeneLevel.RDS"))
) 


## ---------- DIU results (EdgeR) -----------------
TargetedDIU <- readRDS(file = paste0(dirnames$targ_output, "/resultsDIU.RDS"))

## ---------- Expression -----------------
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

## ---------- merged targeted results overview -----------------

Targeted <- list(
  Genes = array(read.table(paste0(dirnames$targ_ont_metadata, "/TargetGenes.tsv"))[["V1"]]),
  
  # final transcript classifications of all merged samples
  Gene_class = lapply(list.files(path = dirnames$targ_anno, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE, full = T), 
                      function(x) read.csv(x)),
  
  # CPAT 
  cpat = read.table(paste0(dirnames$targ_root, "/4_characterise/CPAT/all_iso_ont.ORF_prob.best.tsv"), header = T) %>% 
    mutate(coding_status = ifelse(Coding_prob >= 0.44, "Coding","Non_Coding")),
  
  # noORF file from CPAT
  noORF = read.table(paste0(dirnames$targ_root, "/4_characterise/CPAT/all_iso_ont.no_ORF.txt")) %>% 
    mutate(coding_status = "No_ORF") %>% `colnames<-`(c("seq_ID", "coding_status")),
  
  # reference 
  ref_gencode = read.csv(paste0(dirnames$targ_root, "/4_characterise/TargetGenesRef/TargetGene_Reference_LengthNum.csv")),
  ref_altcon = read.csv(paste0(dirnames$targ_root, "/4_characterise/TargetGenesRef/TargetGene_Reference_AltConExons.csv")),
  
  # ont abundance 
  ont_abundance = class.files$targ_filtered %>% dplyr::select(isoform, contains("ONT")) %>% mutate(annot_transcript_id = isoform)
  
)
names(Targeted$Gene_class) = lapply(list.files(path = dirnames$targ_anno, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE), 
                                    function(x) word(x, c(1), sep = fixed("/")))
Merged_gene_class_df <- all_summarise_gene_stats(Targeted$Gene_class, class.files$targ_filtered, Targeted$cpat, Targeted$noORF, Targeted$Genes)


## ---------------------------
# Isabel's supplementary table of differentially expressed genes in rTg4510 
rnaseq_results <- list(
  AgeGenotypeDEG = read.csv(paste0(dirnames$glob_rnaseq, "/Isabel_Supp4_Tg4510AgeGenotypeDEG.csv"), header = T, as.is = T),
  GenotypeDEG = read.csv(paste0(dirnames$glob_rnaseq, "/Isabel_Supp2_Tg4510GenotypeDEG.csv"), header = T, as.is = T)
)


## ---------------------------
# TAPPAS (Differential Analysis) 
#tappas_dir <- list(
#  glob_iso = paste0(dirnames$glob_tabroot, "/IsoSeq_Expression"), 
#  glob_rna = paste0(dirnames$glob_tabroot, "/RNASeq_Expression"),
#  targ_iso = paste0(dirnames$targ_iso_root, "/thesis_dump/DiffAnalysis_noRNASEQ/TAPPAS_OUTPUT/IsoSeq_Expression"),
#  targ_ont = paste0(dirnames$targ_ont_root, "/thesis_dump/TALON/TAPPAS_OUTPUT")
#)



# Output from Tappas_DEA.R
#tappassiggene <- list(
#  glob = read_dea_files(paste0(dirnames$glob_tabroot,"/DifferentialGeneExpression_Analysis.xlsx")),
#  targ_iso = read_dea_files(paste0(dirnames$targ_iso_root,"/thesis_dump/DiffAnalysis_noRNASEQ/TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx")),
#  targ_ont = read_dea_files(paste0(dirnames$targ_ont_root,"/thesis_dump/TALON/TAPPAS_OUTPUT/DifferentialGeneExpression_Analysis.xlsx"))
#)

#tappassigtrans <- list(
#  glob = read_dea_files(paste0(dirnames$glob_tabroot,"/DifferentialTransExpression_Analysis.xlsx")),
#  targ_iso = read_dea_files(paste0(dirnames$targ_iso_root,"/thesis_dump/DiffAnalysis_noRNASEQ/TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx")),
#  targ_ont = read_dea_files(paste0(dirnames$targ_ont_root,"/thesis_dump/TALON/TAPPAS_OUTPUT/DifferentialTransExpression_Analysis.xlsx"))
#)
