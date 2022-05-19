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
#SBATCH --output=Diff_Targeted_noRNASeq.o
#SBATCH --error=Diff_Targeted_noRNASeq.e

# 22/11/2021: Rerun TargetedIsoSeq but with no RNA-Seq as filter


#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/RNASeq
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

cd $DiffAnalysis_WKD; mkdir -p SQANTI3 SQANTI_TAMA_FILTER RNASeq_SQANTI3 TAPPAS_INPUT COLLAPSE_FILTER WHOLE_TARGETED
cd $DiffAnalysis_WKD/TAPPAS_INPUT; mkdir -p RNASeq_Expression IsoSeq_Expression
#cd $DiffAnalysis_WKD/COLLAPSE_FILTER/; mkdir Expression Length

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
DIFF_FUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis

module load Miniconda2/4.3.21

dataset=AllMouseTargeted
sqname=$dataset".collapsed_classification.filtered_lite"
sqpath=$DiffAnalysis_WKD/SQANTI3/$sqname

################################################################################################
echo "#************************************* RNASeq & IsoSeq [Function 12]"
# rerun Kallisto on the AllMouseTargeted RNA-Seq files using the only the --rf_stranded
## 12) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
#run_kallisto AllMouseTargeted $PostIsoseq3_WKD/TOFU/AllMouseTargeted.collapsed.rep.fa $RNASeq_WKD/MAPPED $DiffAnalysis_WKD/KALLISTO

################################################################################################
echo "#************************************* SQANTI3 [Function 1]"
## 1) run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
# all the samples RNASeq for junction file
SAMPLES_NAMES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)
run_sqanti3 $dataset".collapsed" $dataset".collapsed.gff" $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllMouseTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/AllMouseTargeted.Demultiplexed_Abundance.txt $DiffAnalysis_WKD/SQANTI3 nornaseq

# Rscript script.R <input.classfile> <input.gtf> <output.dir> <prefix>
source activate sqanti2_py3
Rscript $GENERALFUNC/3ISM_remove_classification.R $sqpath"_classification.txt" $sqpath".gtf" $sqpath"_junctions.txt"  $DiffAnalysis_WKD/SQANTI3 $dataset
python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $sqpath".fasta" $DiffAnalysis_WKD/SQANTI3/$dataset"_ISMrem.isoform.txt" $sqpath"_ISMrem.fasta"
