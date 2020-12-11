# Szi Kay Leung
# Aim: To run differential transcript usage using IsoformSwitchAnalyze.R
# 23/11/2020: Using output from SQANTI3 (chained data as input) with FL read count for differential transcript usage
# Prequisite: Run sqanti_gtf_modify.py --> all_samples.chained.rep_classification.filtered_lite_mod.gtf

library(dplyr)
library(IsoformSwitchAnalyzeR)

options(scipen=999) # removal of scientific notation 

# read in SQANTI2 classification file of all merged data
input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples/SQANTI3/"
# sqanti classification file filtering removed monoexonic transcripts 
class <- read.table(paste0(input_dir,"all_samples.chained_classification.filtered_lite_classification.txt"),header=T, as.is = T, sep = "\t")

isoform_count <- class %>% group_by(associated_gene) %>% tally()


# PacBio isoform ID, keep only columns with isoform FL counts (and convert FL counts to TPM)
counts <- me_class[, c("isoform","FL.K17","FL.K18","FL.K24","FL.L22","FL.M21","FL.O18","FL.O23","FL.Q20","FL.Q21","FL.S18","FL.S23")] 
abundance <- cbind(counts[,1], as.data.frame(apply(counts[,-1],2, function(x) round(x/sum(x)*1000000,0))))
colnames(counts)[1] <- c("isoform_id")      # column name important for input into IsoformSwitch
colnames(abundance)[1] <- c("isoform_id")

# using SQANTI generated gtf and fasta, rather than gencode
#isoformExonAnnoation = paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite.gtf"),
aSwitchList <- importRdata(
  isoformCountMatrix   = counts,
  isoformRepExpression = abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = paste0(input_dir,"all_samples.chained_classification.filtered_lite.gtf"),
  isoformNtFasta       = paste0(input_dir,"all_samples.chained_classification.filtered_lite.fasta"),
)



############ from chaining 
input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples/TOFU/"
counts <- read.table(paste0(input_dir,"all_samples.chained_count.txt"),header=T, as.is = T, sep = "\t")
counts[is.na(counts)] <- 0
abundance <- cbind(counts[,1], as.data.frame(apply(counts[,-1],2, function(x) round(x/sum(x)*1000000,0))))
colnames(counts)[1] <- c("isoform_id")      # column name important for input into IsoformSwitch
colnames(abundance)[1] <- c("isoform_id")

# design dataframe input for IsoformSwitchAnalyzeR
myDesign <- data.frame(
  sampleID = colnames(counts)[-1],
  condition = c("WT", "TG", "TG", "TG","WT","TG","WT","TG","WT","TG","WT")
)


# failed misannotation 
#misannotated <- c("PB.429.","PB.436.","PB.5107.","PB.13705","PB.430","PB.6677")
misannotated <- c("PB.13787.","PB.11204.","PB.7815.","PB.14001.","PB.2141.")
counts_mod <- counts[!grepl(paste(misannotated,collapse="|"),counts$isoform_id),]
abundance_mod <- abundance[!grepl(paste(misannotated,collapse="|"),abundance$isoform_id),]

# python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence/fq2fa.py all_samples.chained.rep.fq


# using SQANTI generated gtf and fasta, rather than gencode
#isoformExonAnnoation = paste0(input_dir,"all_samples.chained.rep_classification.filtered_lite.gtf"),
aSwitchList <- importRdata(
  isoformCountMatrix   = counts,
  #isoformRepExpression = abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = paste0(input_dir,"all_samples.chained.gff"),
  isoformNtFasta       = paste0(input_dir,"all_samples.chained.rep.fasta"),
)



# fudged aSwitchList to take the gene_id from the class_list rather than generated gene_id from gtf 
#df <- merge(aSwitchList[["isoformFeatures"]], class[,c("isoform","associated_gene")], by.x = "isoform_id", by.y = "isoform")
#df <- df[,c(2,30,1,30,5:29)]
#colnames(df)[2] <- "gene_ref"
#colnames(df)[4] <- "gene_id" 
#df$gene_ref <- as.character(df$gene_ref)
#df$gene_id <- as.character(df$gene_id)
#aSwitchList[["isoformFeatures"]] <- df


# isoformExpressionCutoff =, filtering on isoform expression allows removal of non-used isoforms that only appear in the switchAnalyzeRlist because they were in the isoform/gene annotation used. Furthermore, the expression filtering allows removal of lowly expressed isoforms where the expression levels might be untrustworthy. The filtering on gene expression allows for removal of lowly expressed genes which causes the calculations of the Isoform Fractions (IF) to become untrustworthy
# removeSingleIsoformGenes = TRUE: Removal of single isoform genes is the default setting in preFilter() since these genes, per definition, cannot have changes in isoform usage.
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 5,
  isoformExpressionCutoff = 2,
  removeSingleIsoformGenes = TRUE)


SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)

SwitchListAnalyzed_Sg <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  alpha = 0.05,
  dIFcutoff = 0.1,
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