## ---------- Script -----------------
##
## Purpose: input variables for ONT targeted mouse transcriptome datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)

suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))


## ---------- Directory and input files -----------------

dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/",
  meta = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/F_ont_targeted/",
  talon = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/thesis_dump/TALON/All/",
  ont_sq = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/7_sqanti3/",
  iso_sq = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/B_IsoSeq_Targeted/thesis_dump/DiffAnalysis_noRNASEQ/SQANTI3/",
  targetgenes = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/F_ONT_Targeted/thesis_dump/TALON/All/Merged/TargetGenes",
  ref = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/4_characterise/TargetGenesRef/"
)


# Miscellaneous input files 
misc_input <- list(
  
  # list of samples
  tg4510_samples = read.csv(paste0(dirnames$meta, "Tg4510_fullsample.csv"))[,c("Genotype","Sample.ID","RIN","ng.ul")],
  
  # list of target genes
  TargetGene = read.table(paste0(dirnames$meta, "TargetGenes.tsv"))[["V1"]],
  
  # list of barcodes for each sample
  BarcodedPhenotype = read.csv(paste0(dirnames$meta, "ONTBarcoded_Phenotype.csv")),
  
  # ONT retained list of isoforms from filtering and merging with Iso-Seq 
  ont_retained_id = read.table(paste0(dirnames$talon,"/Merged/IsoSeqONT_RetainedIDs.txt")),
  
  # abundance file of ONT retained IDs
  ONT_abundance = read.csv(paste0(dirnames$talon,"Merged/ONT_Retained_Abundance.csv")) %>% mutate(FL = rowSums(select_if(., is.numeric), na.rm = TRUE)),
  
  # path of directory with comprehensive annotation of genes 
  target_anno_dir = paste0(dirnames$talon, "Merged/TargetGenes"),
  
  # final transcript classifications of all merged samples
  Gene_class = lapply(list.files(path = dirnames$targetgenes, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE, full = T), 
                      function(x) read.csv(x)),
  
  # CPAT 
  cpat = read.table(paste0(dirnames$talon, "/Merged/SQANTI3/IsoSeqONT_cpat.ORF_prob.best.tsv"), header = T) %>% 
    mutate(coding_status = ifelse(Coding_prob >= 0.44, "Coding","Non_Coding")),
  
  # noORF file from CPAT
  noORF = read.table(paste0(dirnames$talon, "/Merged/SQANTI3/IsoSeqONT_cpat.no_ORF.txt")) %>% 
    mutate(coding_status = "No_ORF") %>% `colnames<-`(c("seq_ID", "coding_status")),
  
  
  # reference 
  ref_gencode = read.csv(paste0(dirnames$ref, "TargetGene_Reference_LengthNum.csv")),
  ref_altcon = read.csv(paste0(dirnames$ref, "TargetGene_Reference_AltConExons.csv"))
  
)

names(misc_input$Gene_class) = lapply(list.files(path = dirnames$targetgenes, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE), 
                           function(x) word(x, c(1), sep = fixed("/")))


## ---------- SQANTI files -----------------

# relevant columns for subsetting
sq_cols = c("isoform", "min_cov","associated_gene", "exons", "length", "dist_to_cage_peak", 
            "within_cage_peak", "dist_to_polya_site","within_polya_site","polyA_motif",
            "structural_category","subcategory","Dataset","diff_to_TTS","diff_to_TSS","RNASeq_supported")

input.class.files <- list(
  
  # ONT TALON unfiltered targeted dataset
  ont_unfiltered = SQANTI_class_preparation(paste0(dirnames$ont_sq,"Unfiltered/basic/ONTTargeted_unfiltered_talon_classification.txt"),"ONT") %>% 
    filter(associated_gene %in% misc_input$TargetGene),
  
  # ONT TALON filtered targeted dataset
  ont_filtered = SQANTI_class_preparation(paste0(dirnames$ont_sq,"Filtered/basic/ONTTargeted_filtered_talon_classification.txt"),"ONT") %>% 
    filter(associated_gene %in% misc_input$TargetGene) %>% mutate(Dataset = "ONT"),
  
  # Iso-Seq
  iso = SQANTI_class_preparation(paste0(dirnames$iso_sq,"AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"),"standard") 
    %>% filter(associated_gene %in% misc_input$TargetGene) %>% mutate(Dataset = "Iso-Seq"),
  
  # Merged ONT and Iso-Seq targeted dataset
  merged = SQANTI_class_preparation(paste0(dirnames$talon, "Merged/SQANTI3/IsoSeqONT_final_genename_classification.txt"), "not_standard")
  
)

# ONT retained targeted dataset (after merging with Iso-Seq)
# Note isoforms may be removed in TALON filtering but kept if also detected in Iso-Seq dataset
input.class.files$ont_retained = input.class.files$ont_unfiltered %>% filter(isoform %in% misc_input$ont_retained_id$V1) %>% mutate(Dataset = "ONT") %>% 
  # merge with abundance file of retained ONT isoforms
  dplyr::select(-FL) %>% merge(., misc_input$ONT_abundance, by.x = "isoform", by.y = "annot_transcript_id")

# Merge of original sqanti classification of the transcripts that were retained in the final curated dataset
input.class.files$concat_retained = annotate_class_binary(rbind(input.class.files$ont_retained[sq_cols], input.class.files$iso[sq_cols]))

# Merged ONT and Iso-Seq targeted dataset, with 3'ISM removed 
input.class.files$merged_noISM = SQANTI_remove_3prime(input.class.files$merged)  

# identify dataset in merged classification file based on isoform ID
input.class.files$merged_noISM <- wrapper_identify_dataset(input.class.files$merged_noISM)

# include abundance of Iso-Seq and ONT-derived isoforms
input.class.files$merged_noISM_exp <- split_iso_ont_characteristics(input.class.files$merged_noISM, input.class.files$iso, input.class.files$ont_retained)