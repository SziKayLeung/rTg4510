#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=STAR_All-%A_%a.o
#SBATCH --array=0-63%8

# 11/12/2020: Run STAR across all 64 RNASeq Tg4510 data for input into SQANTI as junction coverage (from STAR output SJ.bed files)

#************************************* DEFINE VARIABLES
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/RNASeq
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/STAR_main
RAW_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered/Tg4510_filtered
MAPPED=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples/RNASeq

TG_SAMPLES=(K24 L22 M20 O24 P22 Q20 S24	T22	K22	L20	M18	O22	P20	Q18	S22	T20	K20	L18	M24	O20	P18	Q24	S20	T18	K18	L24	M22	O18	P24	Q22	S18	T24)
WT_SAMPLES=(K17 L21 M19 K23 P21 Q19 M21 T21 K21 L19 M17 O21 P19 Q17 S21 T19 K19 L17 M23 O19 P17 Q23 S19 T17 O23 L23 Q21 O17 P23 S23 S17 T23)
ALL_SAMPLES+=( "${TG_SAMPLES[@]}" "${WT_SAMPLES[@]}" )
sample=${ALL_SAMPLES[${SLURM_ARRAY_TASK_ID}]}


#************************************* Run Star in parallel
source $FUNCTIONS/STAR_Functions.sh
# run_star $sample $Tg4510/J20/input_directory $MAPPED_output_directory  $REFERENCE
run_star ${sample} $RAW_DIR $MAPPED/${sample} $REFERENCE
