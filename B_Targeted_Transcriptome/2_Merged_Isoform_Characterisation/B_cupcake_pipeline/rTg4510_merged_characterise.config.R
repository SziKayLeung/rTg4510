## ---------- Script -----------------
##
## Purpose: appending input variables merged datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ----------------------------------

# source rTg4510_ont_characterise.config.R

## ---------- Appending to existing lists -----------------

# directory paths
dirnames$merge_root <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/G_Merged_Targeted/"
dirnames$cupcake_merged <- paste0(dirnames$merge_root,"B_cupcake_pipeline/")
dirnames$targetgenes <- paste0(dirnames$cupcake_merged,"4_characterise/TargetGenes")

# input class.names and abundance files
input.class.names.files$cupcake_merged <- paste0(dirnames$cupcake_merged,"3_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.txt")
input.abundance.files$cupcake_merged <- read.csv(paste0(dirnames$cupcake_merged,"2_collapse/demux_fl_count.csv"), row.names = "isoform")


## ---------- Annotations -----------------

# SQANTI classification files
input.class.files$cupcake_merged <- SQANTI_class_preparation(input.class.names.files$cupcake_merged,"nstandard")
targeted.class.files$cupcake_merged <- subset_class_by_targets(input.class.files$cupcake_merged, misc_input$TargetGene)
targeted.class.files$cupcake_merged <- quantify_class_abundance(targeted.class.files$cupcake_merged, input.abundance.files$cupcake_merged)


# Miscellaneous input files 
merged_misc_input <- list(
  
  # final transcript classifications of all merged samples
  Gene_class = lapply(list.files(path = dirnames$targetgenes, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE, full = T), 
                      function(x) read.csv(x)),
  
  # CPAT 
  cpat = read.table(paste0(dirnames$cupcake_merged, "4_characterise/CPAT/all_iso_ont.ORF_prob.best.tsv"), header = T) %>% 
    mutate(coding_status = ifelse(Coding_prob >= 0.44, "Coding","Non_Coding")),
  
  # noORF file from CPAT
  noORF = read.table(paste0(dirnames$cupcake_merged, "4_characterise/CPAT/all_iso_ont.no_ORF.txt")) %>% 
    mutate(coding_status = "No_ORF") %>% `colnames<-`(c("seq_ID", "coding_status")),
  
  # reference 
  ref_gencode = read.csv(paste0(dirnames$cupcake_merged, "4_characterise/TargetGenesRef/TargetGene_Reference_LengthNum.csv")),
  ref_altcon = read.csv(paste0(dirnames$cupcake_merged, "4_characterise/TargetGenesRef/TargetGene_Reference_AltConExons.csv"))
  
)

names(merged_misc_input$Gene_class) = lapply(list.files(path = dirnames$targetgenes, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE), 
                                      function(x) word(x, c(1), sep = fixed("/")))
misc_input <- append(misc_input, merged_misc_input)
rm(merged_misc_input)

