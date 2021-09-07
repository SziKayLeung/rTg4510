#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=0:30:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-58%8 #59 samples (64-5)
#SBATCH --output=Diff_WholeRNASeq-%A_%a.o
#SBATCH --error=Diff_WholeRNASeq-%A_%a.e

# 29/06/2021: Align all 59 RNA-Seq samples to the Iso-Seq scaffold generated in Diff_Whole.sh (SQ3)

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq

#cd $DiffAnalysis_WKD; mkdir SQANTI3 SQANTI_TAMA_FILTER RNASeq_SQANTI3 TAPPAS_INPUT
#cd $DiffAnalysis_WKD/RNASeq_SQANTI3; mkdir Individual

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered

module load Miniconda2/4.3.21
################################################################################################
echo "#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold"
# all the samples RNASeq (for downstream TAPPAS) except samples K22, L19, L20,P22,O17
SAMPLES_NAMES=(K24 L22 M20 O24 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 M18 O22 P20 Q18 S22 T20 K21 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 P23 S23 S17 T23)
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

# run_kallisto_1sample <input_RNASEQ_rawdir> <sample> <input_ref_name_idx> <output_dir>
run_kallisto_1sample $RNASeq_Filtered ${SAMPLE} AllRNASeq_Kallisto.idx $DiffAnalysis_WKD/RNASeq_SQANTI3/Individual
