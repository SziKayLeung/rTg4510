#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## proteogenonomics analysis following paper review
##    process merged (ONT + Iso-Seq) targeted data after collapsing by protein ORF following G.Shenkyman pipeline
##    re-determined reference isoform using most abundant transcript after collapsing by ORF
##    differential transcript analysis at protein level 
## --------------------------------

## ---------- packages -----------------

suppressMessages(library("utils"))


## ---------- config file -----------------

source("/lustre/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Reviews/rTg4510_config.R")


## ---------- process data from pipeline -----------------

## data-wrangle orf_refined.tsv 
  # filter to target genes
  # generate column of the number of transcripts collapsed by delimiting the pb_accs column
  # output: a table of the pb_accs, base_acc (the representative collapsed isoform selected by G.Shenkyman pipeline) and numtxCollapsed
mouseProtein$t2p.collapse <- mouseProtein$t2p.collapse %>% 
  filter(gene%in%TargetGene) %>% 
  mutate(numtxCollapsed = count.fields(textConnection(as.character(pb_accs)), sep = "|"))
char <- strsplit(as.character(mouseProtein$t2p.collapse $pb_accs), '|', fixed = T)
t2p.collapse.dissected <- data.frame(pb_accs=unlist(char), base_acc=rep(mouseProtein$t2p.collapse$base_acc, sapply(char, FUN=length)))
mouseProtein$t2p.collapse <- merge(t2p.collapse.dissected,mouseProtein$t2p.collapse[,c("base_acc","numtxCollapsed")])

## re-determine representative colalsped isoform: using ONT abundance (sum across all samples) rather than arbitrary (G.Shenkyman pipeline)
 # take the ONT_sum read counts from the classification file
 # max = grouping by the base_acc (i.e. the previously selected isoform), select the rows with the maximum ONT FL reads
 # create an index to remap and create a "corrected_acc" column with the corresponding isoform that has the highest number of ONT FL reads
mouseProtein$t2p.collapse <- merge(mouseProtein$t2p.collapse,class.files$targ_filtered[,c("isoform","ONT_sum_FL")],by.x = "pb_accs", by.y = "isoform", all.x = TRUE)
max = mouseProtein$t2p.collapse %>% group_by(base_acc) %>% filter(ONT_sum_FL == max(ONT_sum_FL))
idx <- match(mouseProtein$t2p.collapse$base_acc, max$base_acc)
mouseProtein$t2p.collapse = transform(mouseProtein$t2p.collapse, corrected_acc = ifelse(!is.na(idx), as.character(max$pb_accs[idx]), base_acc))

## include in the original classification file the collapsed PB.ID 
class.files$targ_filtered <- merge(class.files$targ_filtered, mouseProtein$t2p.collapse[,c("pb_accs","numtxCollapsed","base_acc","corrected_acc")], by.x = "isoform", by.y = "pb_accs", all.x = TRUE)

## Statistics
message("Total number of RNA transcripts: ", nrow(class.files$targ_filtered))
message("Number of protein-coding RNA transcripts: ", nrow(class.files$targ_filtered[!is.na(class.files$targ_filtered$base_acc),]))
message("Number of non-protein-coding RNA transcripts: ", nrow(class.files$targ_filtered[is.na(class.files$targ_filtered$base_acc),]))


## ---------- differential expression analysis: proteogenomics -----------------

## 1. expression file
# aggregate sum by same peptide sequence
# datawrangle for input to run_DESeq2()
pFL <- class.files$targ_filtered %>% filter(!is.na(corrected_acc) & associated_gene %in% TargetGene) 
ontpFL <- pFL %>% select(corrected_acc, contains("ONT"), -ONT_sum_FL)
isopFL <- pFL %>% select(corrected_acc, contains("Iso.Seq"), -Iso.Seq_sum_FL)
expressionFiles <- list(
  ontProtein =  aggregate(. ~ corrected_acc, ontpFL, sum) %>% tibble::column_to_rownames(., var = "corrected_acc"),
  isoProtein =  aggregate(. ~ corrected_acc, isopFL, sum) %>% tibble::column_to_rownames(., var = "corrected_acc")
)

## 2. phenotype file 
# matched number of samples across age for genotype analysis
phenotype$targ_matched_ont <- phenotype$targ_ont %>% filter(!sample %in% c("Q20","Q18"))
input <- list(
  ontPhenotype = phenotype$targ_matched_ont %>% mutate(sample = paste0("ONT_",sample)),
  isoPhenotype = phenotype$targ_iso %>% mutate(sample = paste0("Iso.Seq_",sample))
)


## ---------- Creating DESeq2 object and analysis -----------------

ResTran <- list(
  pisoWaldGenotype = run_DESeq2(test="Wald",expressionFiles$isoProtein,input$isoPhenotype,threshold=10,exprowname=NULL,controlname="CONTROL",design="case_control",interaction="On"),
  pontWaldGenotype = run_DESeq2(test="Wald",expressionFiles$ontProtein,input$ontPhenotype,threshold=10,exprowname=NULL,controlname="CONTROL",design="case_control",interaction="On")
)

annoResTran <- list(
  pisoWaldGenotype = anno_DESeq2(ResTran$pisoWaldGenotype,class.files$targ_filtered,input$isoPhenotype,controlname="CONTROL",level="transcript",sig=0.1),
  pontWaldGenotype = anno_DESeq2(ResTran$pontWaldGenotype,class.files$targ_filtered,input$ontPhenotype,controlname="CONTROL",level="transcript",sig=0.1)
)


## ---------- Output -----------------
write.table(class.files$targ_filtered, paste0(dirnames$targ_root, "/2_sqanti3/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered_pCollapsed.txt"),sep="\t",quote = F)
write.table(mouseProtein$t2p.collapse, paste0(dirnames$mprotein,"/all_iso_ont_orf_refined_collapsed.tsv"),sep="\t",quote = F)
saveRDS(annoResTran, file = paste0(dirnames$targ_output, "/DESeq2ProteinLevel.RDS"))
