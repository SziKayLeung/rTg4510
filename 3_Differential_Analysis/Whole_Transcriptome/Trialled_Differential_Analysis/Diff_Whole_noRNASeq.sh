#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=15:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Diff_WholenoRNASeq.o
#SBATCH --error=Diff_WholenoRNASeq.e

# 28/06/2021: Differential Analysis using SQANTI3 on whole transcriptome
# 05/07/2021: Rerun kallisto as stranded -rf-stranded
# 22/11/2021: SQANTI no RNA-Seq followed by collapse


#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_noRNASEQ
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq

cd $DiffAnalysis_WKD; mkdir SQANTI3 SQANTI_TAMA_FILTER RNASeq_SQANTI3 TAPPAS_INPUT COLLAPSE_FILTER AS
cd $DiffAnalysis_WKD/COLLAPSE_FILTER; mkdir GROUPS
cd $DiffAnalysis_WKD/TAPPAS_INPUT; mkdir IsoSeq_Expression RNASeq_Expression


# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
TAPPASFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis
DIFF_FUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis

module load Miniconda2/4.3.21

dataset=WholeIsoSeq
sqname=$dataset".collapsed_classification.filtered_lite"

################################################################################################
echo "#************************************* SQANTI3 [Function 1]"
## 1) run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
SAMPLES_NAMES=(Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)
run_sqanti3 WholeIsoSeq.collapsed WholeIsoSeq.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/WholeIsoSeq.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt $DiffAnalysis_WKD/SQANTI3 nornaseq

## 2) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
TAMA_remove_fragments $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite.gtf WholeIsoSeq $DiffAnalysis_WKD/SQANTI_TAMA_FILTER

## 12) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
TAMA_sqanti_filter $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/WholeIsoSeq.bed $DiffAnalysis_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" WholeIsoSeq $DiffAnalysis_WKD/SQANTI_TAMA_FILTER

# Script <sqanti_filered_inputfile> <tama_filtered_inputfile> <output_expression_file> <output_finalisoform_file> <type == "Targeted/Whole">
Rscript $DIFF_FUNC/Sqanti_Collapsed_Counts.R $DiffAnalysis_WKD/SQANTI3/$sqname"_classification.txt" $DiffAnalysis_WKD/SQANTI3/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/WholeIsoSeq_sqantisubset.expression.txt $DiffAnalysis_WKD/TAPPAS_INPUT/tama_retained_collapsed_pbid.txt Whole

collapse_filter $DiffAnalysis_WKD/TAPPAS_INPUT/tama_retained_collapsed_pbid.txt_FSMbylength.txt $DiffAnalysis_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $DiffAnalysis_WKD/COLLAPSE_FILTER

################################################################################################
echo "#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold"
# all the samples RNASeq
SAMPLES_NAMES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)

# already merged from previous analysis (Diff_Whole.sh)
## 8) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
#mouse_merge_fastq $RNASeq_Filtered $DiffAnalysis_WKD/RNASeq_SQANTI3 AllRNASeq

## 9) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
#run_kallisto AllRNASeq $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantifiltered_tamafiltered_classification.fasta $DiffAnalysis_WKD/RNASeq_SQANTI3 $DiffAnalysis_WKD/RNASeq_SQANTI3

# individual all the RNASeq samples
# first index the fasta file (note this made in run_kallisto but recreating for ease)
# later align all 59 RNASeq samples separately using Kallisto to IsoSeq scaffold (Diff_WholeRNASeq.sh)
source activate sqanti2
cd $DiffAnalysis_WKD/RNASeq_SQANTI3
kallisto index -i AllRNASeq_Kallisto.idx $DiffAnalysis_WKD/COLLAPSE_FILTER/WholeIsoSeq_sqantisubset_classification.fasta 2> AllRNASeq_Kallisto.index.log

##################################################################################################
echo "#************************************* Alternative Splicing Events for Groups"
# for the individual groups
#cd $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/GENOME/; mkdir GROUP SAMPLE
#cd $PostIsoseq3_WKD/SUPPA/; mkdir GROUP
sqname=WholeIsoSeq_sqantisubset
source activate sqanti2_py3; Rscript $TAPPASFUNC/Counts_Groups.R $DiffAnalysis_WKD/COLLAPSE_FILTER/WholeIsoSeq_sqantisubset.classification.txt $DiffAnalysis_WKD/COLLAPSE_FILTER/GROUPS

for group in WT_2mos TG_2mos WT_8mos TG_8mos; do
  echo "Processing: $group"
  # Rscript .R <retained_id> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_junc_txt> <output_prefix_name> <output_dir>
  #Rscript $GENERALFUNC/sqanti_classgtfsubset.R $DiffAnalysis_WKD/COLLAPSE_FILTER/GROUPS/$group"_counts.txt" $DiffAnalysis_WKD/COLLAPSE_FILTER $sqname".classification.txt" $sqname".classification.gtf" $sqname".junction.txt" $group $DiffAnalysis_WKD/COLLAPSE_FILTER/GROUPS

  # run_suppa2 <input_gtf> <input_class> <output_dir> <output_name>
  run_suppa2 $DiffAnalysis_WKD/COLLAPSE_FILTER/GROUPS/$group"_sqantisubset.classification.gtf" $DiffAnalysis_WKD/COLLAPSE_FILTER/GROUPS/$group"_sqantisubset.classification.txt" $DiffAnalysis_WKD/AS $group
done

for group in WT TG; do Rscript $GENERALFUNC/sqanti_classgtfsubset.R $DiffAnalysis_WKD/COLLAPSE_FILTER/GROUPS/$group"_counts.txt" $DiffAnalysis_WKD/COLLAPSE_FILTER $sqname".classification.txt" $sqname".classification.gtf" $sqname".junction.txt" $group $DiffAnalysis_WKD/COLLAPSE_FILTER/GROUPS; done
