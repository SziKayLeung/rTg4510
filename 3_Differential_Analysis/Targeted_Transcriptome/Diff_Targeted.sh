#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Diff_Targeted.o
#SBATCH --error=Diff_Targeted.e

# 28/06/2021: Differential Analysis using SQANTI3 on targeted transcriptome
# 27/07/2021: Rerun Kallisto using only the --rf-stranded and included collapsed filtering output of transcript counts
  # run kallisto and run collapse filter on command line with logged files
# 16/08/2021: Rerun Kallisto using the collapsed sqanti files and obtain bam files for visualisation on igv



#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/RNASeq

#cd $DiffAnalysis_WKD; mkdir SQANTI3 SQANTI_TAMA_FILTER RNASeq_SQANTI3 TAPPAS_INPUT COLLAPSE_FILTER Whole_Targeted
#cd $DiffAnalysis_WKD/RNASeq_SQANTI3; mkdir Individual Individual_Collapsed_FSMExp Individual_Collapsed_FSMLength Whole_Targeted
#cd $DiffAnalysis_WKD/TAPPAS_INPUT; mkdir RNASeq_Expression IsoSeq_Expression
#cd $DiffAnalysis_WKD/COLLAPSE_FILTER/; mkdir Expression Length

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
DIFF_FUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis

module load Miniconda2/4.3.21

dataset=AllMouseTargeted

################################################################################################
echo "#************************************* RNASeq & IsoSeq [Function 12]"
# rerun Kallisto on the AllMouseTargeted RNA-Seq files using the only the --rf_stranded
## 12) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto AllMouseTargeted $PostIsoseq3_WKD/TOFU/AllMouseTargeted.collapsed.rep.fa $RNASeq_WKD/MAPPED $DiffAnalysis_WKD/KALLISTO

################################################################################################
echo "#************************************* SQANTI3 [Function 1]"
## 1) run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
# all the samples RNASeq for junction file
SAMPLES_NAMES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)
run_sqanti3 $dataset".collapsed" $dataset".collapsed.gff" $DiffAnalysis_WKD/TOFU $RNASeq_WKD/MAPPED $DiffAnalysis_WKD/KALLISTO/AllMouseTargeted.mod.abundance.tsv $DiffAnalysis_WKD/TOFU/AllMouseTargeted.Demultiplexed_Abundance.txt $DiffAnalysis_WKD/SQANTI3 genome

## 2) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
TAMA_remove_fragments $DiffAnalysis_WKD/SQANTI3/$dataset".collapsed_classification.filtered_lite.gtf" $dataset $DiffAnalysis_WKD/SQANTI_TAMA_FILTER

## 12) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=$dataset".collapsed_classification.filtered_lite"
TAMA_sqanti_filter $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/$dataset".bed" $DiffAnalysis_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $DiffAnalysis_WKD/SQANTI_TAMA_FILTER


##################################################################################################
echo "#************************************* Prepare files for TappAS"
# counts_subset_4tappas <input_class> <output_class> <type_genes>
# regenerate expression matrix for tappas from long reads FL read counts
counts_subset_4tappas $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/$dataset"_sqantitamafiltered.expression.txt" AD

# Script <sqanti_filered_inputfile> <tama_filtered_inputfile> <output_expression_file> <output_finalisoform_file> <type == "Targeted/Whole">
Rscript $DIFF_FUNC/Sqanti_Collapsed_Counts.R $DiffAnalysis_WKD/SQANTI3/$sqname"_classification.txt" $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/AllMouseTargeted_sqantisubset.expression $DiffAnalysis_WKD/TAPPAS_INPUT/Retained_collapsed_pbid Targeted

Rscript $DIFF_FUNC/Sqanti_Collapsed_Counts.R $DiffAnalysis_WKD/SQANTI3/$sqname"_classification.txt" $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/AllMouseTargeted_allsqantisubset.expression.txt $DiffAnalysis_WKD/TAPPAS_INPUT/tama_retained_allcollapsed_pbid.txt All

## 12) collapse_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
# rerun on command line though logged files
# note tama_retained_collapsed_pbid.txt same as Retained_collapsed_pbid_FSMbylength.txt
sqname=$dataset".collapsed_classification.filtered_lite"
collapse_filter $DiffAnalysis_WKD/TAPPAS_INPUT/tama_retained_collapsed_pbid.txt $DiffAnalysis_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $DiffAnalysis_WKD/COLLAPSE_FILTER/Length

collapse_filter $DiffAnalysis_WKD/TAPPAS_INPUT/Retained_collapsed_pbid_FSMbyexpression.txt $DiffAnalysis_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $DiffAnalysis_WKD/COLLAPSE_FILTER/Expression

################################################################################################
echo "#************************************* Whole + Targeted Transcriptome"
WHOLE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/COLLAPSE_FILTER/WholeIsoSeq_sqantisubset.final.classification.gtf
# convertgtf2bed <name> <input_gtf> <working_directory>
convertgtf2bed WholeIsoSeq $WHOLE $DiffAnalysis_WKD/Whole_Targeted
convertgtf2bed TargetedIsoSeq $DiffAnalysis_WKD/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.classification.gtf $DiffAnalysis_WKD/Whole_Targeted

## 19) run_tamamerge <wkd> <wholebed> <targetedbed> <outputname>
run_tamamerge $DiffAnalysis_WKD/Whole_Targeted $DiffAnalysis_WKD/Whole_Targeted/WholeIsoSeq.mod_mod.bed12 $DiffAnalysis_WKD/Whole_Targeted/TargetedIsoSeq.mod_mod.bed12 merged_whole_targeted

# run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/rnaseq/lncrna>
run_sqanti3 merged_whole_targeted merged_whole_targeted.gtf $DiffAnalysis_WKD/Whole_Targeted NA NA NA $DiffAnalysis_WKD/Whole_Targeted basic
# no need for sqanti filtering as already filtered before merging datasets
cd $DiffAnalysis_WKD/Whole_Targeted; rm *filtered*
# grep G30.1 -A 10 $DiffAnalysis_WKD/Whole_Targeted/merged_whole_targeted_corrected.fasta
# grep PB.18.1 $DiffAnalysis_WKD/SQANTI3/AllMouseTargeted.collapsed_corrected.fasta

# reannotate the output from tama merge with the original sqanti classification files to match tama_id with pb_id
# also generate an inclusion list from the tama_merged file of isoforms associated with target gene to be used for tappAS input for RNA-Seq expression (to speed up process)
# note this script is hard-coded with the classification files etc
Rscript $DIFF_FUNC/MergeTargetedWholeTranscriptome.R

###############################################################################################
echo "#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold"
# all the samples RNASeq
SAMPLES_NAMES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)

## 8) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
#mouse_merge_fastq $RNASeq_Filtered $DiffAnalysis_WKD/RNASeq_SQANTI3 AllRNASeq

## 9) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
#run_kallisto AllRNASeq $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/$dataset"_sqantifiltered_tamafiltered_classification.fasta" $DiffAnalysis_WKD/RNASeq_SQANTI3 $DiffAnalysis_WKD/RNASeq_SQANTI3

# individual all the RNASeq samples
# first index the fasta file (note this made in run_kallisto but recreating for ease)
# later align all 59 RNASeq samples separately using Kallisto to IsoSeq scaffold (Diff_TargetedRNASeq.sh)
source activate sqanti2
#cd $DiffAnalysis_WKD/RNASeq_SQANTI3/Individual
#kallisto index -i AllRNASeq_Kallisto.idx $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/$dataset"_sqantifiltered_tamafiltered_classification.fasta" 2> AllRNASeq_Kallisto.index.log

#### Targeted Transcriptome
chrom_set=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/mm10/mm10.chromosome

# rerun on command line though logged files
cd $DiffAnalysis_WKD/RNASeq_SQANTI3/Individual_Collapsed_FSMLength
#kallisto index -i AllRNASeq_Kallisto.idx $DiffAnalysis_WKD/COLLAPSE_FILTER/Length/$dataset"_sqantisubset_classification.fasta" 2> AllRNASeq_Kallisto.index.log
#kallisto quant -i AllRNASeq_Kallisto.idx -b 30 -o AllRNASeq --genomebam --gtf $DiffAnalysis_WKD/COLLAPSE_FILTER/Length/AllMouseTargeted_sqantisubset.final.classification.gtf --chromosomes $chrom_set $DiffAnalysis_WKD/RNASeq_SQANTI3/AllRNASeq_R1.fq $DiffAnalysis_WKD/RNASeq_SQANTI3/AllRNASeq_R2.fq
#cd ALLRNASeq; mv pseudoalignments.bam AllRNASeq_FSMLegth_pseudoalignments.bam
## knight ## scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/RNASeq_SQANTI3/Individual_Collapsed_FSMLength/AllRNASeq/AllRNASeq_FSMLegth_pseudoalignments.bam /mnt/data1/Szi/IGV/

# rerun on command line though logged files
cd $DiffAnalysis_WKD/RNASeq_SQANTI3/Individual_Collapsed_FSMExp
#kallisto index -i AllRNASeq_Kallisto.idx $DiffAnalysis_WKD/COLLAPSE_FILTER/Expression/$dataset"_sqantisubset_classification.fasta" 2> AllRNASeq_Kallisto.index.log
kallisto quant -i AllRNASeq_Kallisto.idx -b 30 -o AllRNASeq --genomebam --gtf $DiffAnalysis_WKD/COLLAPSE_FILTER/Expression/AllMouseTargeted_sqantisubset.final.classification.gtf $DiffAnalysis_WKD/RNASeq_SQANTI3/AllRNASeq_R1.fq $DiffAnalysis_WKD/RNASeq_SQANTI3/AllRNASeq_R2.fq
