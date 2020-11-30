#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=16 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 07/10/2019: run_rsem_prepare using RSEM.sh functions script (only need to be run once)
# 07/10/2019: run RSEM.sh functions script on Tg4510 1-16 samples, and all WT and TG merged, and WT8 merged

#************************************* DEFINE GLOBAL VARIABLES
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/RNASeq
RSEM_Reference=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/RSEM
RSEM=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/RSEM
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
INPUT_RNASeq=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq

#************************************* RUN FUNCTIONS SCRIPT 
module load Miniconda2/4.3.21
source activate sqanti2

SAMPLES_NAMES=(Q21 O18 C21 E18 C20 B21 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)

# Testing out RSEM_Functions.sh in RSEM branch 
cd $FUNCTIONS
git checkout RSEM

source $FUNCTIONS/RSEM_Functions.sh

#cd $RSEM_Reference
#run_rsem_prepare

# RSEM output for individual samples
for i in ${SAMPLES_NAMES[@]}; do

    F_name=$(find $INPUT_RNASeq -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R1" )
    R_name=$(find $INPUT_RNASeq -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R2" )
    # save path directory of files as variable for later mapping
    F_File=$(find $INPUT_RNASeq -name "$F_name")
    R_File=$(find $INPUT_RNASeq -name "$R_name")

    #run_rsem <output_sample_prefix> <path_input_R1_fastq.filtered> <path_input_R2_fastq_filtered> <output_dir>
    run_rsem $i $F_File $R_File $RSEM
done
