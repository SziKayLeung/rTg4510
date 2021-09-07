# Szi Kay Leung
# Aim: To run differential transcript usage using IsoformSwitchAnalyze.R
# 16/12/2020: Using SQANTI2 data from Targeted Transcriptome 24 Samples (Using Merged data followed by demutiplexing)

########## READ IN NEWFILE WITH Q21
##### Analysis #####
# 1. Input SQANTI2 classification file, and prepare counts and abundance 
# 2. Importdata into IsoformSwitchAnalyzeR
# 3. Filtering on Gene Expression and Isoform Expression level
# 4. IsoformSwitchTestDEXSeq: Statiscal test for identifying Isoform Switching
#*******************


# Packages
library(dplyr)
library(IsoformSwitchAnalyzeR)
options(scipen=999) # removal of scientific notation 

input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/SQANTI_OLD/"
collapsed_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged/TOFU/"

##### 1. Input #####
##### Input SQANTI2 classification file, and prepare counts and abundance

# Read in SQANTI2 classification file of all merged data
class <- read.table(paste0(input_dir,"All_Targeted_Merged.collapsed.filtered_classification.filtered_lite_classification.txt"),header=T)

# Prepare Counts: keep only PacBio isoform ID and columns with isoform FL counts 
WT <- paste("FL", c("K20","K18", "K24","L22","O18","Q20","S18","O22","Q18","L18","T20","T18"),sep=".")
TG <- paste("FL", c("K23","K21", "K19","K17","O23","M21","S23","Q21","P19","Q17","S19"),sep=".")
count_cols <- do.call(c, list("isoform", WT, TG))

counts <- class[,which((names(class) %in% count_cols)==TRUE)]
colnames(counts)[1] <- c("isoform_id")

# Prepare Abundance: Determine total FL reads and calculate TPM from counts 
abundance <- cbind(counts[,1], as.data.frame(apply(counts[,-1],2, function(x) round(x/sum(x)*1000000,0))))
colnames(abundance)[1] <- c("isoform_id")

##### 2. ImportData #####
##### ImportData into IsoformSwitchAnalyzeR
# design dataframe input for IsoformSwitchAnalyzeR
myDesign <- data.frame(
  sampleID = colnames(counts)[-1],
  condition = c(rep("WT", 12), rep("TG",11))
)

# using SQANTI generated gtf and fasta, rather than gencode
#isoformExonAnnoation = paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite.gtf"),
aSwitchList <- importRdata(
  isoformCountMatrix   = counts,
  isoformRepExpression = abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = paste0(input_dir,"All_Targeted_Merged.collapsed.filtered_classification.filtered_lite.gtf"),
  isoformNtFasta       = paste0(input_dir,"All_Targeted_Merged.collapsed.filtered_classification.filtered_lite.fasta"),
)



##### 3. Filtering #####
##### Filtering on Gene Expression and Isoform Expression level
# isoformExpressionCutoff = removal of lowly-expressed isoforms with potential untrustworthy expression levels
# geneExpressionCutoff = removal of lowly-expressed genes, which cause untrustworthy calculations of the Isoform Fractions (IF) 
# removeSingleIsoformGenes = removal of single isoform genes that cannot have changes in isoform usage
# The filtering removed 36475 ( 95.98% of ) transcripts. There is now 1526 isoforms left
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 0,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)


##### 4. IsoformSwitch #####
##### IsoformSwitchTestDEXSeq: Statiscal test for identifying Isoform Switching
# reduceToSwitchingGenes = reduce/subset the switchAnalyzeRlist to genes which each contains at least one differential used isoform, as indicated by the alpha and dIFcutoff cutoffs, to speed up rest of workflow  
# reduceFurtherToGenesWithConsequencePotential = FALSE; i.e. do not further reduce/subset genes that only have differential used isoforms with consequences (i.e. protein domains) for curiousity 
# alpha = cutoff for fdr correct p-values to determine significant switches 
# dIFcutoff = cutoff for changes in absolute isoform usage to be considered switching
SwitchListAnalyzed_Sg <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  alpha = 0.05, # Default
  dIFcutoff = 0.01, # Default (10%)
  onlySigIsoforms = FALSE
)


extractSwitchSummary(SwitchListAnalyzed)
SwitchListAnalyzed <- analyzeORF(
  SwitchListAnalyzed,
  orfMethod = "longest",
  showProgress=FALSE
)

SwitchListAnalyzed <- extractSequence(
  SwitchListAnalyzed, 
  writeToFile=TRUE,
  alsoSplitFastaFile=TRUE,
  filterAALength=TRUE
)


#pfam
# sed 's/^.//' all_pfam.txt > all_pfam_mod.txt
SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = SwitchListAnalyzed,
  pathToPFAMresultFile = "all_pfam_mod.txt",
)

### Add SignalP analysis
SwitchListAnalyzed <- analyzeSignalP(
  switchAnalyzeRlist       = SwitchListAnalyzed,
  pathToSignalPresultFile  = "all_output_protein_type.txt"
)


SwitchListAnalyzed <- analyzeCPAT(
  switchAnalyzeRlist   = SwitchListAnalyzed,
  pathToCPATresultFile = "cpat.txt",
  codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

# alternative splicing
SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchListAnalyzed,
  quiet=TRUE
)

consequencesOfInterest <- c(
  'intron_retention',
  'coding_potential',
  'ORF_seq_similarity',
  'NMD_status',
  'domains_identified'
)

SwitchListAnalyzed <- analyzeSwitchConsequences(
  SwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest,
  showProgress=FALSE,
)

extractSwitchSummary(SwitchListAnalyzed, filterForConsequences = FALSE)
extractSwitchSummary(SwitchListAnalyzed, dIFcutoff = 0.4, filterForConsequences = FALSE)
extractTopSwitches(
  SwitchListAnalyzedSubset, 
  filterForConsequences = TRUE, 
  n = NA, 
  sortByQvals = FALSE
)

SwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
  SwitchListAnalyzed, 
  SwitchListAnalyzed$isoformFeatures$condition_1 == 'TG'
)

SwitchListAnalyzedSubset 
extractTopSwitches(
  SwitchListAnalyzedSubset, 
  filterForConsequences = TRUE, 
  n = 2, 
  sortByQvals = TRUE
)

switchPlotTopSwitches(
  switchAnalyzeRlist = SwitchListAnalyzed, 
  n = Inf,                                          
  filterForConsequences = FALSE,
  fileType = "pdf",                                 
  pathToOutput = "dtu-plots"
)

extractSwitchOverlap(
  switchAnalyzeRlist = SwitchListAnalyzed,
  filterForConsequences=TRUE
)

extractConsequenceSummary(
  SwitchListAnalyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)

extractConsequenceEnrichment(
  SwitchListAnalyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)

extractSplicingEnrichment(
  SwitchListAnalyzed,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)

ggplot(data=SwitchListAnalyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

### Switch vs Gene changes:
ggplot(data=SwitchListAnalyzed$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + 
  facet_wrap(~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()